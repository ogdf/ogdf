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

#	include <ogdf/geometric/cr_min/geometry/objects/Line.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Point.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Vector.h>

#	include <CGAL/IO/io.h>
#	include <CGAL/Ray_2.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace geometry {
using namespace tools;

template<typename kernel>
using Ray_t = CGAL::Ray_2<kernel>;

template<typename kernel, typename type>
inline Point_t<kernel> rel_point_at(const Ray_t<kernel>& r, const type l) {
	return r.source() + r.to_vector() * l;
}

template<typename kernel>
inline typename kernel::FT squared_convex_parameter_directed(const Ray_t<kernel>& r,
		const Point_t<kernel>& p) {
	OGDF_ASSERT(r.has_on(p));
	return (r.source() - p).squared_length() / r.to_vector().squared_length();
}

template<typename kernel>
inline Point_t<kernel> intersect(const Line_t<kernel>& l, const Ray_t<kernel>& r) {
	OGDF_ASSERT(!l.is_degenerate());
	OGDF_ASSERT(!r.is_degenerate());
	auto result = CGAL::intersection(l, r);
	if (!result) {
		return Point_t<kernel>(std::numeric_limits<unsigned int>::max(),
				std::numeric_limits<unsigned int>::max());
	}
	Point_t<kernel> intersection;
	if (boost::get<Ray_t<kernel>>(&*result)) {
		auto s = boost::get<Ray_t<kernel>>(&*result);
		intersection = s->point(0);
	} else {
		intersection = *boost::get<Point_t<kernel>>(&*result);
	}
	return intersection;
}

} //namespace
}
}
}

#endif
