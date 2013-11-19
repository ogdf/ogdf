/*
 * $Revision: 3837 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 15:19:30 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of TLP string conversion functions.
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

#include <ogdf/fileformats/Tlp.h>


namespace ogdf {

namespace tlp {


std::string toString(const Attribute &attr)
{
	switch(attr) {
	case a_label: return "viewLabel";
	case a_color: return "viewColor";
	case a_position: return "viewLayout";
	case a_size: return "viewSize";
	case a_shape: return "viewShape";
	default: return "unknown";
	}
}


Attribute toAttribute(const std::string &str)
{
	if(str == "viewLabel") {
		return a_label;
	}
	if(str == "viewColor") {
		return a_color;
	}
	if(str == "viewLayout") {
		return a_position;
	}
	if(str == "viewSize") {
		return a_size;
	}
	if(str == "viewShape") {
		return a_shape;
	}
	return a_unknown;
}


} // end namespace tlp

} // end namespace ogdf

