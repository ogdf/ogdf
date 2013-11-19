/*
 * $Revision: 3837 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 15:19:30 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of UCINET DL format parser class.
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

#include <ogdf/fileformats/DLParser.h>


namespace ogdf {


DLParser::DLParser(std::istream &is) : m_istream(is)
{
	init();
}


void DLParser::init()
{
	m_initialized = false;
	m_nodeId.resize(1, NULL);

	m_embedded = false;
	m_nodes = -1;
	m_format = fullmatrix;
}


bool DLParser::initGraph(Graph &G)
{
	G.clear();

	if(m_nodes < 0) {
		std::cerr << "ERROR: Node count not specified or incorrect.\n";
		return false;
	}

	for(int i = 0; i < m_nodes; i++) {
		m_nodeId.push_back(G.newNode());
	}
	m_initialized = true;

	return true;
}


// Common for both embedded and non-embedded mode.
static inline bool readMatrixRow(
	std::istream &is,
	Graph &G, GraphAttributes *GA, node v)
{
	const long attrs = GA ? GA->attributes() : 0;
	const bool iweight = (attrs & GraphAttributes::edgeIntWeight) != 0;
	const bool dweight = (attrs & GraphAttributes::edgeDoubleWeight) != 0;

	node u;
	forall_nodes(u, G) {
		double weight;
		if(!(is >> weight)) {
			std::cerr << "ERROR: Expected matrix value.\n";
			return false;
		}

		edge e = NULL;
		if(weight != 0) {
			e = G.newEdge(v, u);
		}

		if(e && iweight) {
			GA->doubleWeight(e) = weight;
		} else if(e && dweight) {
			GA->intWeight(e) = static_cast<int>(weight);
		}
	}

	return true;
}


bool DLParser::readMatrix(Graph &G, GraphAttributes *GA)
{
	node v;
	forall_nodes(v, G) {
		if(!readMatrixRow(m_istream, G, GA, v)) {
			return false;
		}
	}

	return true;
}


bool DLParser::readEmbeddedMatrix(Graph &G, GraphAttributes *GA)
{
	node v;

	// First, top-label line.
	forall_nodes(v, G) {
		std::string label;
		if(!(m_istream >> label)) {
			std::cerr << "ERROR: Expected node embedded label.\n";
			return false;
		}
		toLower(label);

		if(GA && (GA->attributes() & GraphAttributes::nodeLabel)) {
			GA->label(v) = label;
		}
		m_nodeLabel[label] = v;
	}

	// Now, each row have a label and then "normal" row.
	for(int i = 0; i < G.numberOfNodes(); i++) {
		std::string label;
		if(!(m_istream >> label)) {
			std::cerr << "ERROR: Expected node embedded label.\n";
			return false;
		}
		toLower(label);

		if(!(v = m_nodeLabel[label])) {
			std::cerr << "ERROR: Node with given label \""
					  << label << "\" not found.\n";
			return false;
		}

		if(!readMatrixRow(m_istream, G, GA, v)) {
			return false;
		}
	}

	return true;
}


// Common for both embedded and non-emedded mode.
static inline bool readEdgeListRow(
	std::istringstream &is,
	Graph &G, GraphAttributes *GA, node v, node u)
{
	edge e = G.newEdge(v, u);
	double weight;
	if(GA && (is >> weight)) {
		if(GA->attributes() & GraphAttributes::edgeDoubleWeight) {
			GA->doubleWeight(e) = weight;
		} else if(GA->attributes() & GraphAttributes::edgeIntWeight) {
			GA->intWeight(e) = static_cast<int>(weight);
		}
	}

	return true;
}


inline node DLParser::requestLabel(
	GraphAttributes *GA, node &nextFree,
	const std::string &label)
{
	node v = m_nodeLabel[label];

	if(!v) {
		m_nodeLabel[label] = v = nextFree;
		if(GA && (GA->attributes() & GraphAttributes::nodeLabel)) {
			GA->label(v) = label;
		}
		nextFree = nextFree->succ();
	}

	return v;
}


bool DLParser::readEdgeList(Graph &G, GraphAttributes *GA)
{
	std::string buffer;
	for(size_t line = 1; std::getline(m_istream, buffer); line++) {
		// Not necessary I guess, but does not do any harm.
		if(buffer.empty()) {
			continue;
		}

		std::istringstream is(buffer);
		int vid, uid;

		if(!(is >> vid >> uid) || !fineId(vid) || !fineId(uid)) {
			std::cerr << "ERROR: Node id incorrect (data line "
					  << line << "), maximum value is "
					  << m_nodeId.size() - 1 << ".\n";
			return false;
		}

		readEdgeListRow(is, G, GA, m_nodeId[vid], m_nodeId[uid]);
	}

	return true;
}


bool DLParser::readEmbeddedEdgeList(Graph &G, GraphAttributes *GA)
{
	std::string buffer;

	node nextFree = G.firstNode();
	for(size_t line = 1; std::getline(m_istream, buffer); line++) {
		if(buffer.empty()) {
			continue;
		}
		std::istringstream is(buffer);

		std::string vlabel, ulabel;
		if(!(is >> vlabel >> ulabel)) {
			std::cerr << "ERROR: Expected embedded node labels (data line "
					  << line << "), got \"" << is.str() << "\".\n";
			return false;
		}

		node v = requestLabel(GA, nextFree, vlabel);
		node u = requestLabel(GA, nextFree, ulabel);
		readEdgeListRow(is, G, GA, v, u);
	}

	return true;
}


bool DLParser::readNodeList(Graph &G, GraphAttributes *GA)
{
	std::string buffer;
	for(size_t line = 1; std::getline(m_istream, buffer); line++) {
		std::istringstream is(buffer);

		// As always, either ingore incorrect line or throw error.
		int vid;
		if(!(is >> vid)) {
			continue;
		}

		if(!fineId(vid)) {
			std::cerr << "ERROR: Node id incorrect (data line "
					  << line << ".\n";
			return false;
		}
		node v = m_nodeId[vid];

		int uid;
		while((is >> uid)) {
			if(!fineId(uid)) {
				std::cerr << "ERROR: Node id incorrect (data line "
						  << line << ").\n";
				return false;
			}

			G.newEdge(v, m_nodeId[uid]);
		}
	}

	return true;
}


bool DLParser::readEmbeddedNodeList(Graph &G, GraphAttributes *GA)
{
	std::string buffer;

	node nextFree = G.firstNode();
	for(size_t line = 1; std::getline(m_istream, buffer); line++) {
		std::istringstream is(buffer);

		std::string vlabel;
		if(!(is >> vlabel)) {
			continue;
		}

		node v = requestLabel(GA, nextFree, vlabel);

		std::string ulabel;
		while((is >> ulabel)) {
			node u = requestLabel(GA, nextFree, ulabel);
			G.newEdge(v, u);
		}
	}

	return true;
}


bool DLParser::readData(Graph &G, GraphAttributes *GA)
{
	if(m_nodes < 0) {
		std::cerr << "ERROR: Number of nodes not specified or incorret.\n";
		return false;
	}

	if(!m_initialized) {
		initGraph(G);
	}

	// Now, depending on the method choosen we actually read the graph.
	switch(m_format) {
	case fullmatrix:
		return m_embedded ? readEmbeddedMatrix(G, GA) : readMatrix(G, GA);
	case edgelist:
		return m_embedded ? readEmbeddedEdgeList(G, GA) : readEdgeList(G, GA);
	case nodelist:
		return m_embedded ? readEmbeddedNodeList(G, GA) : readNodeList(G, GA);
	}

	return false;
}


/*
 * This function is quite ugly (sphagetti code all the way). That is because
 * it is trying to mirror shitty DL format design. It is terrible, seriously.
 */
bool DLParser::readWithLabels(Graph &G, GraphAttributes *GA)
{

	std::string buffer;

	initGraph(G);
	for(node v = G.firstNode(); v;) {
		if(!(m_istream >> buffer)) {
			std::cerr << "ERROR: Expected node labels.\n";
			return false;
		}
		toLower(buffer); // Labels should be lowercase.

		// We check whether we need to end reading labels.
		if(buffer == "data:") {
			return readData(G, GA);
		} else if(buffer == "labels") {
			// Or we have "labels embedded" information.
			m_istream >> buffer;
			toLower(buffer);
			if(buffer != "embedded:" && buffer != "embedded") {
				std::cerr << "ERROR: Expected embedded keyword, got \""
						  << buffer << "\".\n";
				return false;
			}

			m_embedded = true;
			break;
		}

		// We split input via comma and read labels for succesive nodes.
		std::istringstream is(buffer);
		while(std::getline(is, buffer, ',')) {
			// There is no need parsing labels if GA is not given.
			if(GA && (GA->attributes() & GraphAttributes::nodeLabel)) {
				GA->label(v) = buffer;
			}
			m_nodeLabel[buffer] = v;
			v = v->succ();
		}
	}

	m_istream >> buffer;
	toUpper(buffer);

	if(buffer == "LABELS") {
		m_istream >> buffer;
		toUpper(buffer);
		if(buffer != "EMBEDDED:" && buffer != "EMBEDDED") {
			std::cerr << "ERROR: Expected \"EMBEDDED\" keyword, got \""
					  << buffer << "\".\n";
			return false;
		}

		m_embedded = true;
		m_istream >> buffer;
		toUpper(buffer);
	}

	if(buffer != "DATA:") {
		std::cerr << "ERROR: Expected \"DATA:\" statement, got \""
				  << buffer << "\".\n";
		return false;
	}

	return readData(G, GA);
}


bool DLParser::readAssignment(
	Graph &G,
	const std::string &lhs, const std::string &rhs)
{

	if(lhs == "N") {
		std::istringstream is(rhs);
		if(!(is >> m_nodes)) {
			std::cerr << "ERROR: Incorrect number of nodes.\n";
			return false;
		}
	} else if(lhs == "FORMAT") {
		if(rhs == "FULLMATRIX") {
			m_format = fullmatrix;
		} else if(rhs == "EDGELIST1") {
			m_format = edgelist;
		} else if(rhs == "NODELIST1") {
			m_format = nodelist;
		} else {
			std::cerr << "ERROR: Unknown data format \"" << rhs << "\".\n";
			return false;
		}
	} else {
		std::cerr << "ERROR: Unkown assignment statement: "
				  << "\"" << lhs << "\".\n";
		return false;
	}

	return true;
}


bool DLParser::readStatements(Graph &G, GraphAttributes *GA)
{
	std::string buffer;

	if(!(m_istream >> buffer)) {
		std::cerr << "ERROR: Expected statement.\n";
		return false;
	}
	toUpper(buffer);

	if(buffer == "DATA:") {
		return readData(G, GA);
	}

	if(buffer == "LABELS:") {
		return readWithLabels(G, GA);
	}

	if(buffer == "LABELS") {
		m_istream >> buffer;
		toUpper(buffer);
		if(buffer != "EMBEDDED" && buffer != "EMBEDDED:") {
			std::cerr << "ERROR: Unknown statement "
					  << "\"LABELS " << buffer << "\". "
					  << "Did you mean \"LABELS:\" or \"LABELS EMBEDDED\"?"
					  << "\n";
			return false;
		}

		m_embedded = true;

		// ... and here we go again.
		return readStatements(G, GA);
	}

	// If none of the above, try interpreting this as assignment statement.
	size_t eq = buffer.find('=');
	std::string lhs, rhs;
	if(eq == std::string::npos) {
		// '=' not found inside, therefore buffer has to be left side.
		lhs = buffer;
		char c;
		if(!(m_istream >> c) || c != '=') {
			std::cerr << "ERROR: Expected definition or assignment "
					  << "statement, got: \"" << lhs << "\".\n";
			return false;
		}

		if(!(m_istream >> rhs)) {
			std::cerr << "ERROR: Expected assignment right side.\n";
			return false;
		}
	} else if(eq == buffer.size() - 1) {
		// 'lhs= rhs' case.
		if(!(m_istream >> rhs)) {
			std::cerr << "ERROR: Expected assignment right side.\n";
			return false;
		}
		lhs = buffer.substr(0, eq);
	} else {
		// 'lhs=rhs' case.
		lhs = buffer.substr(0, eq);
		rhs = buffer.substr(eq + 1);
	}
	toUpper(lhs);
	toUpper(rhs);

	if(!readAssignment(G, lhs, rhs)) {
		return false;
	}

	return readStatements(G, GA);
}


bool DLParser::readGraph(Graph &G, GraphAttributes *GA)
{
	init();
	std::string buffer;

	m_istream >> buffer;
	toUpper(buffer);

	if(buffer != "DL") {
		std::cerr << "ERROR: Expected the \"DL\" header, got: \""
				  << buffer << "\".\n";
	}

	return readStatements(G, GA);
}


} // end namespace ogdf

