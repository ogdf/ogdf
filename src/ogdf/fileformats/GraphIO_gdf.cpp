/*
 * $Revision: 3837 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 15:19:30 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Implements GDF write functionality of class GraphIO.
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

#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/fileformats/GDF.h>


namespace ogdf {

namespace gdf {


static inline void writeColor(std::ostream &os, const Color &color)
{
	int r = color.red(), g = color.green(), b = color.blue();
	std::ostringstream ss;
	os << "\"" << r << "," << g << "," << b << "\"";
}


static inline void writeNodeHeader(
	std::ostream &os,
	const GraphAttributes *GA)
{
	os << "nodedef>";
	os << toString(na_name);

	const long attrs = GA ? 0 : GA->attributes();
	if(attrs & GraphAttributes::nodeLabel) {
		os << "," << toString(na_label);
	}
	if(attrs & GraphAttributes::nodeGraphics) {
		os << "," << toString(na_x);
		os << "," << toString(na_y);
		if(attrs & GraphAttributes::threeD) {
			os << "," << toString(na_z);
		}
		os << "," << toString(na_shape);
		os << "," << toString(na_width);
		os << "," << toString(na_height);
	}
	if(attrs & GraphAttributes::nodeStyle) {
		os << "," << toString(na_fillColor);
		os << "," << toString(na_strokeColor);
	}
	if(attrs & GraphAttributes::nodeTemplate) {
		os << "," << toString(na_template);
	}
	if(attrs & GraphAttributes::nodeWeight) {
		os << "," << toString(na_weight);
	}

	os << "\n";
}


static inline void writeNode(
	std::ostream &os,
	const GraphAttributes *GA, node v)
{
	/*
	 * According to official documentation, it is preferred to not name nodes
	 * as number literals and preceed it with some string if needed.
	 */
	os << "n" << v->index();

	const long attrs = GA ? 0 : GA->attributes();
	if(attrs & GraphAttributes::nodeLabel) {
		os << "," << GA->label(v);
	}
	if(attrs & GraphAttributes::nodeGraphics) {
		os << "," << GA->x(v);
		os << "," << GA->y(v);
		if(attrs & GraphAttributes::threeD) {
			os << "," << GA->z(v);
		}
		os << "," << toString(GA->shape(v));
		os << "," << GA->width(v);
		os << "," << GA->height(v);
	}
	if(attrs & GraphAttributes::nodeStyle) {
		os << ",";
		writeColor(os, GA->fillColor(v));
		os << ",";
		writeColor(os, GA->strokeColor(v));
	}
	if(attrs & GraphAttributes::nodeTemplate) {
		os << "," << GA->templateNode(v);
	}
	if(attrs & GraphAttributes::nodeWeight) {
		os << "," << GA->weight(v);
	}

	os << "\n";
}


static inline void writeEdgeHeader(
	std::ostream &os,
	const GraphAttributes *GA)
{
	os << "edgedef>";
	os << toString(ea_source);
	os << "," << toString(ea_target);
	if(GA && GA->directed()) {
		os << "," << toString(ea_directed);
	}

	const long attrs = GA ? 0 : GA->attributes();
	if(attrs & GraphAttributes::edgeLabel) {
		os << "," << toString(ea_label);
	}
	if(attrs & (GraphAttributes::edgeIntWeight |
		        GraphAttributes::edgeDoubleWeight))
	{
		os << "," << toString(ea_weight);
	}
	if(attrs & GraphAttributes::edgeStyle) {
		os << "," << toString(ea_color);
	}
	if(attrs & GraphAttributes::edgeGraphics) {
		os << "," << toString(ea_bends);
	}

	os << "\n";
}


static inline void writeEdge(
	std::ostream &os,
	const GraphAttributes *GA, edge e)
{
	os << "n" << e->source()->index() << ","
	   << "n" << e->target()->index();
	if(GA && GA->directed()) {
		os << "," << "true";
	}

	const long attrs = GA ? 0 : GA->attributes();
	if(attrs & GraphAttributes::edgeLabel) {
		os << "," << GA->label(e);
	}
	if(attrs & GraphAttributes::edgeDoubleWeight) {
		os << "," << GA->doubleWeight(e);
	} else if(attrs & GraphAttributes::edgeIntWeight) {
		os << "," << GA->intWeight(e);
	}
	if(attrs & GraphAttributes::edgeStyle) {
		os << ",";
		writeColor(os, GA->strokeColor(e));
	}
	if(attrs & GraphAttributes::edgeGraphics) {
		os << "," << "\"";

		bool comma = false;
		forall_listiterators(DPoint, it, GA->bends(e)) {
			if(comma) {
				os << ",";
			}
			comma = true;

			const DPoint &p = *it;
			os << p.m_x << "," << p.m_y;
		}

		os << "\"";
	}

	os << "\n";
}


static void writeGraph(
	std::ostream &os,
	const Graph &G, const GraphAttributes *GA)
{
	// Node definition section.
	writeNodeHeader(os, GA);

	node v;
	forall_nodes(v, G) {
		writeNode(os, GA, v);
	}

	// Edge definition section.
	writeEdgeHeader(os, GA);

	edge e;
	forall_edges(e, G) {
		writeEdge(os, GA, e);
	}
}


} // end namespace gdf


bool GraphIO::writeGDF(const Graph &G, std::ostream &os)
{
	gdf::writeGraph(os, G, NULL);
	return true;
}


bool GraphIO::writeGDF(const GraphAttributes &GA, std::ostream &os)
{
	gdf::writeGraph(os, GA.constGraph(), &GA);
	return true;
}


} // end namespace ogdf

