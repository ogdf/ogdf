/*
 * $Revision: 3837 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 15:19:30 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Declarations for GDF file format
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

#ifdef _MSC_VER
#pragma once
#endif


#ifndef OGDF_GDF_H
#define OGDF_GDF_H


#include <ogdf/basic/graphics.h>
#include <ogdf/basic/Hashing.h>
#include <ogdf/basic/HashArray.h>

#include <string>


namespace ogdf {

namespace gdf {


enum NodeAttribute {
	// GDF standard
	na_name = 0,
	na_label,
	na_x, na_y, na_z,
	na_fillColor, na_strokeColor,
	na_shape,
	na_width, na_height,
	// OGDF specific
	na_template,
	na_weight,
	na_unknown
};


enum EdgeAttribute {
	// GDF standard
	ea_label = 0,
	ea_source, ea_target,
	ea_weight,
	ea_directed,
	ea_color,
	// OGDF specific
	ea_bends,
	ea_unknown
};


std::string toString(const NodeAttribute &attr);
std::string toString(const EdgeAttribute &attr);
std::string toString(const Shape &shape);

NodeAttribute toNodeAttribute(const std::string &str);
EdgeAttribute toEdgeAttribute(const std::string &str);
Shape toShape(const std::string &str);


} // end namespace gdf

} // end namespace ogdf


#endif
