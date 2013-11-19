/*
 * $Revision: 3837 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 15:19:30 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of DOT string conversion functions.
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

#include <ogdf/fileformats/DOT.h>
#include <ogdf/fileformats/Utils.h>

namespace ogdf {

namespace dot {


std::string toString(const Attribute &attr)
{
	switch(attr) {
	case a_id: return "id";
	case a_label: return "label";
	case a_template: return "comment";
	case a_width: return "width";
	case a_height: return "height";
	case a_shape: return "shape";
	case a_position: return "pos";
	case a_stroke: return "color";
	case a_fill: return "fillcolor";
	case a_weight: return "weight";
	case a_arrow: return "arrow";
	default: return "comment";
	}
}


std::string toString(const Shape &shape)
{
	switch(shape) {
	case shRect: return "rect";
	case shRoundedRect: return "rect"; // Not supported.
	case shEllipse: return "ellipse";
	case shTriangle: return "triangle";
	case shPentagon: return "pentagon";
	case shHexagon: return "hexagon";
	case shOctagon: return "octagon";
	case shRhomb: return "diamond";
	case shTrapeze: return "trapezium";
	case shParallelogram: return "parallelogram";
	case shInvTriangle: return "invtriangle";
	case shInvTrapeze: return "invtrapezium";
	case shInvParallelogram: return "parallelogram"; // Not supported.
	case shImage: return "box"; // Not supported.
	default: return "rect";
	}
}


std::string toString(const EdgeArrow &arrow)
{
	switch(arrow) {
	case eaNone: return "none";
	case eaLast: return "forward";
	case eaFirst: return "back";
	case eaBoth: return "both";
	case eaUndefined: return "none"; // Not supported.
	default: return "none";
	}
}


std::string toString(const Graph::EdgeType &type)
{
	// Based on IBM UML documentation:
	// http://publib.boulder.ibm.com/infocenter/rsahelp/v7r0m0/index.jsp?topic=
	// /com.ibm.xtools.modeler.doc/topics/crelsme_clssd.html
	switch(type) {
	case Graph::association: return "none";
	case Graph::generalization: return "empty";
	case Graph::dependency: return "open";
	default: return "normal";
	}
}


// Map is lazily-evaluated (this could be avoided with C++11 constexpr).
static Hashing<std::string, Attribute> *attrMap = NULL;

Attribute toAttribute(const std::string &str)
{
	return toEnum(
		str, attrMap, toString,
		static_cast<Attribute>(0), a_unknown, a_unknown);
}


// Same as attrMap but with shapes.
static Hashing<std::string, Shape> *shapeMap = NULL;
Shape toShape(const std::string &str) {
	return toEnum(
		str, shapeMap, toString,
		shRect, shImage, shRect);
}

// Same as attrMap but with arrows.
static Hashing<std::string, EdgeArrow> *arrowMap = NULL;

EdgeArrow toArrow(const std::string &str)
{
	return toEnum(
		str, arrowMap, toString,
		eaNone, eaUndefined, eaUndefined);
}


} // end namespace graphml

} // end namespace ogdf
