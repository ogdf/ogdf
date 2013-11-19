/*
 * $Revision: 3837 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 15:19:30 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of GEXF string conversion functions.
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

#include <ogdf/fileformats/GEXF.h>


namespace ogdf {

namespace gexf {


std::string toString(const Shape &shape)
{
	switch(shape) {
	case shRect: return "square";
	case shRoundedRect: return "rect"; // Not supported.
	case shEllipse: return "disc";
	case shTriangle: return "triangle";
	case shRhomb: return "diamond";
	case shImage: return "image";
	default: return "disc";
	}
}


Shape toShape(const std::string &str)
{
	if(str == "square") {
		return shRect;
	} else if(str == "disc") {
		return shEllipse;
	} else if(str == "triangle") {
		return shTriangle;
	} else if(str == "diamond") {
		return shRhomb;
	} else if(str == "image") {
		return shImage;
	} else {
		return shRect;
	}
}


} // end namespace gexf

} // end namespace ogdf

