/*
 * $Revision: 3837 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 15:19:30 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of GDF format parsing utilities
 *
 * \author Łukasz Hanuszczak
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.txt in the root directory of the OGDF installation for details.
 *
 * \par
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * Version 2 or 3 as published by the Free Software Foundation;
 * see the file LICENSE.txt included in the packaging of this file
 * for details.
 *
 * \par
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * \par
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/

#include <ogdf/fileformats/GdfParser.h>
#include <ogdf/fileformats/Utils.h>

namespace ogdf {

namespace gdf {


Parser::Parser(std::istream &is) : m_istream(is), m_nodeId(NULL)
{
}


/*
 * A method to avoid code duplication: reads either a "nodedef" or
 * "edgedef" GDF file header.
 */
template <typename Attr>
static inline bool readDef(
	const std::string &str,
	Attr toAttribute(const std::string &str), Attr a_unknown,
	std::vector<Attr> &attrs)
{
	std::istringstream ss(str);
	std::string attr;

	/*
	 * Okay, so there is a attribute name first and then optional type.
	 * Therefore we chop this stream by commas, and then we read first
	 * non-whitespace string out of it.
	 */
	while(std::getline(ss, attr, ',')) {
		std::istringstream attrss(attr);
		std::string name;
		attrss >> name;

		Attr attr = toAttribute(name);
		if(attr == a_unknown) {
			std::cerr << "WARNING: attribute \"" << name << "\""
			          << " not supported. Ignoring.\n";
		}
		attrs.push_back(attr);
	}

	return true;

}


bool Parser::readNodeDef(const std::string &str)
{
	return readDef(str, toNodeAttribute, na_unknown, m_nodeAttrs);
}


bool Parser::readEdgeDef(const std::string &str)
{
	return readDef(str, toEdgeAttribute, ea_unknown, m_edgeAttrs);
}


static bool scanQuoted(
	const std::string &str, size_t pos,
	std::string &buff)
{
	for(size_t j = 1; pos + j < str.length(); j++) {
		if(str[pos] == str[pos + j] && str[pos + j - 1] != '\\') {
			return true; // was j before, but j is always >= 1
		}
		buff += str[pos + j];
	}

	return false;
}


static bool split(
	const std::string &str,
	std::vector<std::string> &result)
{
	result.clear();
	std::string buff = "";

	const size_t len = str.length();
	for(size_t i = 0; i < len; i++) {
		if(str[i] == '"' || str[i] == '\'') {
			size_t quoted = scanQuoted(str, i, buff);
			if(quoted) {
				i += quoted;
			} else {
				std::cerr << "ERROR: Unescaped quote.\n";
				return false;
			}
		} else if(str[i] == ',') {
			result.push_back(buff);
			buff = "";
		} else {
			buff += str[i];
		}
	}

	// Last buffer is not inserted during the loop.
	result.push_back(buff);

	return true;
}


bool Parser::readNodeStmt(
	Graph &G, GraphAttributes *GA,
	const std::string &str, size_t line)
{
	std::vector<std::string> values;
	split(str, values);

	if(values.size() != m_nodeAttrs.size()) {
		std::cerr << "ERROR: node definition does not match the header "
		          << "(line " << line << ").\n";
		return false;
	}

	node v = G.newNode();
	for(size_t i = 0; i < values.size(); i++) {
		if(m_nodeAttrs[i] == na_name) {
			m_nodeId[values[i]] = v;
		}
	}

	if(GA && !readAttributes(*GA, v, values)) {
		return false;
	}

	return true;
}


bool Parser::readEdgeStmt(
	Graph &G, GraphAttributes *GA,
	const std::string &str, size_t line)
{
	std::vector<std::string> values;
	split(str, values);

	if(values.size() != m_edgeAttrs.size()) {
		std::cerr << "ERROR: edge definition does not match the header "
		          << "(line " << line << ").\n";
		return false;
	}

	// First, we scan a list for source, target and edge direction.
	bool directed = false;
	node source = NULL, target = NULL;
	for(size_t i = 0; i < values.size(); i++) {
		switch(m_edgeAttrs[i]) {
		case ea_directed:
			if(values[i] == "true") {
				directed = true;
			} else if(values[i] == "false") {
				directed = false;
			} else {
				std::cerr << "ERROR: edge direction must be a boolean "
				          << "(line " << line << ").\n";
			}
			break;
		case ea_source:
			source = m_nodeId[values[i]];
			break;
		case ea_target:
			target = m_nodeId[values[i]];
			break;
		default:
			break;
		}
	}

	// Then, we can create edge(s) and read attributes (if needed).
	if(!source || !target) {
		std::cerr << "ERROR: source or target for edge not found "
		          << "(line " << line << ").\n";
		return false;
	}

	edge st = G.newEdge(source, target);
	edge ts = directed ? NULL : G.newEdge(target, source);

	if(GA && st && !readAttributes(*GA, st, values)) {
		return false;
	}

	if(GA && ts && !readAttributes(*GA, ts, values)) {
		return false;
	}

	return true;
}


static inline Color toColor(const std::string &str)
{
	std::istringstream is(str);
	int r, g, b;
	is >> r >> ',' >> g >> ',' >> b;

	return Color(r, g, b);
}


static bool inline readAttribute(
	GraphAttributes &GA, node v,
	const NodeAttribute &attr, const std::string &value)
{
	const long attrs = GA.attributes();
	switch(attr) {
	case na_name:
		// Not really an attribute, handled elsewhere.
		break;
	case na_label:
		if(attrs & GraphAttributes::nodeLabel) {
			GA.label(v) = value;
		}
		break;
	case na_x:
		if(attrs & GraphAttributes::nodeGraphics) {
			std::istringstream is(value);
			is >> GA.x(v);
		}
		break;
	case na_y:
		if(attrs & GraphAttributes::nodeGraphics) {
			std::istringstream is(value);
			is >> GA.y(v);
		}
		break;
	case na_z:
		if(attrs & GraphAttributes::threeD) {
			std::istringstream is(value);
			is >> GA.z(v);
		}
		break;
	case na_fillColor:
		if(attrs & GraphAttributes::nodeStyle) {
			GA.fillColor(v) = toColor(value);
		}
		break;
	case na_strokeColor:
		if(attrs & GraphAttributes::nodeStyle) {
			GA.strokeColor(v) = toColor(value);
		}
		break;
	case na_shape:
		if(attrs & GraphAttributes::nodeGraphics) {
			GA.shape(v) = toShape(value);
		}
		break;
	case na_width:
		if(attrs & GraphAttributes::nodeGraphics) {
			std::istringstream ss(value);
			ss >> GA.width(v);
		}
		break;
	case na_height:
		if(attrs & GraphAttributes::nodeGraphics) {
			std::istringstream ss(value);
			ss >> GA.height(v);
		}
		break;
	case na_template:
		if(attrs & GraphAttributes::nodeTemplate) {
			GA.templateNode(v) = value;
		}
		break;
	case na_weight:
		if(attrs & GraphAttributes::nodeWeight) {
			std::istringstream ss(value);
			ss >> GA.weight(v);
		}
		break;
	default:
		break;
	}
	return true;
}


static bool inline readAttribute(
	GraphAttributes &GA, edge e,
	const EdgeAttribute &attr, const std::string &value)
{
	const long attrs = GA.attributes();

	switch(attr) {
	case ea_label:
		if(attrs & GraphAttributes::edgeLabel) {
			GA.label(e) = value;
		}
		break;
	case ea_source:
		// Handled elsewhere.
		break;
	case ea_target:
		// Handled elsewhere.
		break;
	case ea_directed:
		// Handled elsewhere.
		break;
	case ea_weight:
		if(attrs & GraphAttributes::edgeDoubleWeight) {
			std::istringstream ss(value);
			ss >> GA.doubleWeight(e);
		} else if(attrs & GraphAttributes::edgeIntWeight) {
			std::istringstream ss(value);
			ss >> GA.intWeight(e);
		}
		break;
	case ea_color:
		if(attrs & GraphAttributes::edgeStyle) {
			GA.strokeColor(e) = toColor(value);
		}
		break;
	case ea_bends:
		if(attrs & GraphAttributes::edgeGraphics) {
			std::istringstream ss(value);
			std::string x, y;

			DPolyline &line = GA.bends(e);
			line.clear();
			while(std::getline(ss, x, ',') && std::getline(ss, y, ',')) {
				std::istringstream conv;
				double dx, dy;

				conv.clear();
				conv.str(x);
				conv >> dx;

				conv.clear();
				conv.str(y);
				conv >> dy;

				line.pushBack(DPoint(dx, dy));
			}
		}
		break;
	default:
		break;
	}

	return true;
}


/*
 * Once again, generic \i readAttributes method to avoid code duplication.
 */
template <typename T, typename A>
static inline bool readAttrs(
	GraphAttributes &GA, T elem,
	const std::vector<A> &attrs,
	const std::vector<std::string> &values)
{
	for(size_t i = 0; i < values.size(); i++) {
		if(!readAttribute(GA, elem, attrs[i], values[i])) {
			return false;
		}
	}

	return true;
}


bool Parser::readAttributes(
	GraphAttributes &GA, node v,
	const std::vector<std::string> &values)
{
	return readAttrs(GA, v, m_nodeAttrs, values);
}


bool Parser::readAttributes(
	GraphAttributes &GA, edge e,
	const std::vector<std::string> &values)
{
	return readAttrs(GA, e, m_edgeAttrs, values);
}


/*
 * Just checks wheter beginning of the string is equal to the pattern.
 * Returns number of matched letters (length of the pattern).
 */
size_t match(const std::string &text, const std::string pattern) {
	const size_t len = pattern.length();
	if(len > text.length()) {
		return 0;
	}

	for(size_t i = 0; i < len; i++) {
		if(pattern[i] != text[i]) {
			return 0;
		}
	}

	return len;
}


bool gdf::Parser::readGraph(
	Graph &G, GraphAttributes *GA)
{
	enum { m_none, m_node, m_edge } mode = m_none;

	size_t line = 0;
	std::string str;
	while(std::getline(m_istream, str)) {
		line += 1;

		/*
		 * We skip empty lines (it is not stated in documentation whether they
		 * are allowed or not, but I like empty lines so since it causes no
		 * charm I think they are fine).
		 */
		if(str.empty()) {
			continue;
		}

		size_t matched = 0;
		if((matched = match(str, "nodedef>"))) {
			if(!readNodeDef(str.substr(matched))) {
				return false;
			}
			mode = m_node;
		} else if((matched = match(str, "edgedef>"))) {
			if(!readEdgeDef(str.substr(matched))) {
				return false;
			}
			mode = m_edge;
		} else if(mode == m_node) {
			if(!readNodeStmt(G, GA, str, line)) {
				return false;
			}
		} else if(mode == m_edge) {
			if(!readEdgeStmt(G, GA, str, line)) {
				return false;
			}
		} else {
			std::cerr << "ERROR: Expected node or edge definition header "
			          << "(line " << line << ").\n";
			return false;
		}
	}

	return true;
}


} // end namespace gdf

} // end namespace ogdf

