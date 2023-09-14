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

#	include <ogdf/geometric/cr_min/geometry/objects/Point.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Vector.h>

#	include <CGAL/Line_2.h>
#	include <math.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace geometry {
using namespace tools;
template<typename kernel>
using Line_t = CGAL::Line_2<kernel>;

template<typename kernel>
inline bool turn(const Line_t<kernel>& l1, const Line_t<kernel>& l2) {
	return turn(l1.to_vector(), l2.to_vector());
}

template<typename kernel>
inline bool left_turn(const Line_t<kernel>& l1, const Line_t<kernel>& l2) {
	return left_turn(l1.to_vector(), l2.to_vector());
}

template<typename kernel>
inline bool right_turn(const Line_t<kernel>& l1, const Line_t<kernel>& l2) {
	return right_turn(l1.to_vector(), l2.to_vector());
}

template<typename kernel>
inline bool turn(const Line_t<kernel>& l, const Point_t<kernel>& p) {
	return turn(l.to_vector(), p - l.point(0));
}

template<typename kernel>
inline bool left_turn(const Line_t<kernel>& l, const Point_t<kernel>& p) {
	return left_turn(l.to_vector(), p - l.point(0));
}

template<typename kernel>
inline bool right_turn(const Line_t<kernel>& l, const Point_t<kernel>& p) {
	return right_turn(l.to_vector(), p - l.point(0));
}

template<typename kernel>
inline Point_t<kernel> intersect(const Line_t<kernel>& l1, const Line_t<kernel>& l2) {
	OGDF_ASSERT(!l1.is_degenerate());
	OGDF_ASSERT(!l2.is_degenerate());
	auto result = CGAL::intersection(l1, l2);
	if (!result) {
		return Point_t<kernel>(std::numeric_limits<unsigned int>::max(),
				std::numeric_limits<unsigned int>::max());
	}

	Point_t<kernel> intersection;
	if (boost::get<Line_t<kernel>>(&*result)) {
		auto s = boost::get<Line_t<kernel>>(&*result);
		intersection = s->point(0);
	} else {
		intersection = *boost::get<Point_t<kernel>>(&*result);
	}
	return intersection;
}

} /* namespace geometry */
}
}
}

#endif
