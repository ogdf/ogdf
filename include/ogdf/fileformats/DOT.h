/*
 * $Revision: 3837 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 15:19:30 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief DOT related enums and string conversion functions.
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


#ifndef OGDF_DOT_H
#define OGDF_DOT_H


#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/HashArray.h>
#include <string>


namespace ogdf {

namespace dot {

	enum Attribute {
		a_id = 0,
		a_label,
		a_template,
		a_stroke,
		a_fill,
		a_width,
		a_height,
		a_shape,
		a_weight,
		a_position,
		a_arrow,
		a_unknown
	};

	std::string toString(const Attribute &attr);
	std::string toString(const Shape &shape);
	std::string toString(const EdgeArrow &arrow);

	Attribute toAttribute(const std::string &str);
	Shape toShape(const std::string &str);
	EdgeArrow toArrow(const std::string &str);

} // end namespace dot

} // end namespace ogdf


#endif
