/*
 * $Revision: 3837 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 15:19:30 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief String conversions and Hashing for GDF fileformat
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

#include <ogdf/fileformats/GDF.h>
#include <ogdf/fileformats/Utils.h>


namespace ogdf {

namespace gdf {


std::string toString(const NodeAttribute &attr)
{
	switch(attr) {
	case na_name: return "name";
	case na_x: return "x";
	case na_y: return "y";
	case na_z: return "z";
	case na_fillColor: return "color";
	case na_strokeColor: return "strokecolor";
	case na_shape: return "style";
	case na_width: return "width";
	case na_height: return "height";
	case na_label: return "label";
	case na_template: return "template";
	case na_weight: return "weight";
	case na_unknown: return "unknown";
	}

	return "";
}


std::string toString(const EdgeAttribute &attr)
{
	switch(attr) {
	case ea_label: return "label";
	case ea_source: return "node1";
	case ea_target: return "node2";
	case ea_weight: return "weight";
	case ea_directed: return "directed";
	case ea_color: return "color";
	case ea_bends: return "bends";
	case ea_unknown: return "unknown";
	}

	return "";
}


std::string toString(const Shape &shape)
{
	/*
	 * Based on official documentation:
	 * http://guess.wikispot.org/The_GUESS_.gdf_format
	 */
	switch(shape) {
	case shRect: return "1";
	case shEllipse: return "2";
	case shRoundedRect: return "3";
	case shImage: return "7";
	default: return "1";
	}
}


static Hashing<std::string, NodeAttribute> *nodeAttrMap = NULL;

NodeAttribute toNodeAttribute(const std::string &str)
{
	return toEnum(
		str, nodeAttrMap, toString,
		static_cast<NodeAttribute>(0), na_unknown, na_unknown);
}


static Hashing<std::string, EdgeAttribute> *edgeAttrMap = NULL;

EdgeAttribute toEdgeAttribute(const std::string &str)
{
	return toEnum(
		str, edgeAttrMap, toString,
		static_cast<EdgeAttribute>(0), ea_unknown, ea_unknown);
}


static Hashing<std::string, Shape> *shapeMap = NULL;

Shape toShape(const std::string &str)
{
	return toEnum(
		str, shapeMap, toString,
		shRect, shImage, shRect);
}



} // end namespace gdf

} // end namespace ogdf

