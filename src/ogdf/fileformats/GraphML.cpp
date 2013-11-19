/*
 * $Revision: 3837 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 15:19:30 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of GraphML string conversion functions.
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

#include <ogdf/fileformats/GraphML.h>


namespace ogdf {

namespace graphml {


std::string toString(const Attribute &attr)
{
	switch(attr) {
	// (moved down to default case)
	//case a_unknown: return "unknown";

	case a_nodeLabel: return "label";
	case a_edgeLabel: return "edgelabel";

	case a_x: return "x";
	case a_y: return "y";
	case a_z: return "z";
	case a_width: return "width";
	case a_height:	return "height";
	case a_size: return "size";

	case a_shape: return "shape";

	case a_nodeStroke: return "nodestroke";
	case a_edgeStroke: return "edgestroke";
	case a_clusterStroke: return "clusterstroke";
	case a_nodeFill: return "nodefill";
	case a_r: return "r";
	case a_g: return "g";
	case a_b: return "b";

	case a_nodeWeight: return "nodeweight";
	case a_edgeWeight: return "weight";

	case a_nodeType: return "nodetype";
	case a_edgeType: return "edgetype";

	case a_template: return "template";

	case a_edgeArrow: return "arrow";
	case a_edgeSubGraph: return "avaliable-for";
	case a_edgeBends: return "bends";

	// default case (to avoid compiler warnings)
	default: return "unknown";
	}
}


std::string toString(const Shape &shape)
{
	switch(shape) {
	case shRect: return "rect";
	case shRoundedRect: return "rounded-rect";
	case shEllipse: return "ellipse";
	case shTriangle: return "triangle";
	case shPentagon: return "pentagon";
	case shHexagon: return "hexagon";
	case shOctagon: return "octagon";
	case shRhomb: return "rhomb";
	case shTrapeze: return "trapeze";
	case shParallelogram: return "parallelogram";
	case shInvTriangle: return "inv-triangle";
	case shInvTrapeze: return "inv-trapeze";
	case shInvParallelogram: return "inv-parallelogram";
	case shImage: return "image";
	default: return "rect";
	}
}


std::string toString(const EdgeArrow &arrow)
{
	switch(arrow) {
	case eaNone: return "none";
	case eaLast: return "last";
	case eaFirst: return "first";
	case eaBoth: return "both";
	case eaUndefined: return "undefined";
	default: return "undefined";
	}
}


std::string toString(const Graph::NodeType &type)
{
	switch(type) {
	case Graph::vertex: return "vertex";
	case Graph::dummy: return "dummy";
	case Graph::generalizationMerger: return "generalization-merger";
	case Graph::generalizationExpander: return "generalization-expander";
	case Graph::highDegreeExpander: return "high-degree-expander";
	case Graph::lowDegreeExpander: return "low-degree-expander";
	case Graph::associationClass: return "association-class";
	default: return "vertex";
	}
}


std::string toString(const Graph::EdgeType &type)
{
	switch(type) {
	case Graph::association: return "association";
	case Graph::generalization: return "generalization";
	case Graph::dependency: return "dependency";
	default: return "association";
	}
}


template <typename E>
static inline E toEnum(
	const std::string &str, // A string we want to convert.
	Hashing<std::string, E> *&map, // A map to be lazily evaluated.
	const E first, const E last, const E def) // Enum informations.
{
	if(!map) {
		map = new Hashing<std::string, E>();

		// Iterating over enums is potentially unsafe... (fixable in C++11).
		for(int it = first; it <= last; it++) {
			const E e = static_cast<E>(it);
			map->insert(toString(e), e);
		}
	}

	HashElement<std::string, E> *elem = map->lookup(str);
	return elem ? elem->info() : def;
}


// Map is lazily-evaluated (this could be avoided with C++11 constexpr).
static Hashing<std::string, Attribute> *attrMap = NULL;

Attribute toAttribute(const std::string &str)
{
	return toEnum(
		str, attrMap,
		static_cast<Attribute>(0), a_unknown, a_unknown);
}


// Same as attrMap but with shapes.
static Hashing<std::string, Shape> *shapeMap = NULL;

Shape toShape(const std::string &str)
{
	return toEnum(str, shapeMap, shRect, shImage, shRect);
}


// Same as attrMap but with arrows.
static Hashing<std::string, EdgeArrow> *arrowMap = NULL;

EdgeArrow toArrow(const std::string &str)
{
	return toEnum(str, arrowMap, eaNone, eaUndefined, eaUndefined);
}


// Same as attrMap but with node types.
static Hashing<std::string, Graph::NodeType> *nodeTypeMap = NULL;

Graph::NodeType toNodeType(const std::string &str)
{
	return toEnum(
		str, nodeTypeMap,
		static_cast<Graph::NodeType>(0), Graph::vertex, Graph::vertex);
}

// Same as attrMap but with edge types.
static Hashing<std::string, Graph::EdgeType> *edgeTypeMap = NULL;

Graph::EdgeType toEdgeType(const std::string &str)
{
	return toEnum(
		str, edgeTypeMap,
		static_cast<Graph::EdgeType>(0), Graph::dependency, Graph::association);
}


} // end namespace graphml

} // end namespace ogdf

