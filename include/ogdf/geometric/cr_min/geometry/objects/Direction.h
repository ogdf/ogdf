/** \file
 *
 * \author Marcel Radermacher
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.md in the OGDF root directory for details.
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
 * License along with this program; if not, see
 * http://www.gnu.org/copyleft/gpl.html
 */

#pragma once

#ifdef OGDF_INCLUDE_CGAL

#	include <CGAL/Aff_transformation_2.h>
#	include <CGAL/Direction_2.h>
#	include <CGAL/aff_transformation_tags.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace geometry {
template<typename kernel>
using Direction_t = CGAL::Direction_2<kernel>;

template<typename kernel>
inline Direction_t<kernel> rotate(const Direction_t<kernel>& v, const double angle) {
	CGAL::Aff_transformation_2<kernel> rotatation(CGAL::ROTATION, std::sin(angle), std::cos(angle));
	return std::move(rotatation(v));
}
}
}
}
}

#endif
