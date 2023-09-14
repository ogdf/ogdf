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

#	include <ogdf/geometric/cr_min/geometry/objects/Vector.h>
#	include <ogdf/geometric/cr_min/tools/math.h>

#	include <limits>
#	include <sstream>

#	include <CGAL/Aff_transformation_2.h>
#	include <CGAL/Point_2.h>
#	include <CGAL/aff_transformation_tags.h>
#	include <CGAL/squared_distance_2.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace geometry {

template<typename kernel>
using Point_t = CGAL::Point_2<kernel>;

template<typename kernel>
inline bool is_the_same(const Point_t<kernel>& p1, const Point_t<kernel>& p2) {
	return p1.x() == p2.x() && p1.y() == p2.y();
}

template<typename kernel>
inline typename kernel::FT distance(const Point_t<kernel>& p1, const Point_t<kernel>& p2) {
	return CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(p1, p2)));
}

template<typename kernel, typename precision>
inline Point_t<kernel> rotate(const Point_t<kernel>& p, precision angle) {
	CGAL::Aff_transformation_2<kernel> rotate(CGAL::ROTATION, sin(angle), cos(angle));
	return std::move(rotate(p));
}

template<typename kernel>
inline bool is_valid(const Point_t<kernel>& p) {
	//return !isinf(p.x()) && !isinf(p.y()) && !isnan(p.x()) && !isnan(p.y());
	return CGAL::is_finite(p.x()) && CGAL::is_finite(p.y()) && CGAL::is_valid(p.x())
			&& CGAL::is_valid(p.y());
}

template<typename kernel>
void serialize(const Point_t<kernel>& p, std::ostream& os) {
	os.write(reinterpret_cast<const char*>(&p.x()), sizeof(p.x()));
	os.write(reinterpret_cast<const char*>(&p.y()), sizeof(p.y()));
}

template<typename kernel>
void deserialize(Point_t<kernel>& p, std::istream& in) {
	in.read(reinterpret_cast<char*>(&p.x()), sizeof(p.x()));
	in.read(reinterpret_cast<char*>(&p.y()), sizeof(p.y()));
}

template<typename Kernel>
struct ExactPointComparison {
	using Point = Point_t<Kernel>;

	static bool equal(const Point& p1, const Point& p2) { return is_the_same(p1, p2); }

	static bool less(const Point& p1, const Point& p2) {
		return p1.x() < p2.x() || (p1.x() == p2.x() && p1.y() < p2.y());
	}

	static bool less_equal(const Point& p1, const Point& p2) {
		return p1.x() < p2.x() || (p1.x() == p2.x() && p1.y() <= p2.y());
	}

	static bool greater(const Point& p1, const Point& p2) {
		return p1.x() > p2.x() || (p1.x() == p2.x() && p1.y() > p2.y());
	}

	static bool greater_equal(const Point& p1, const Point& p2) {
		return p1.x() > p2.x() || (p1.x() == p2.x() && p1.y() >= p2.y());
	}

	bool operator()(const Point& p1, const Point& p2) const { return less(p1, p2); }
};

template<typename Kernel, typename T>
inline geometry::Point_t<Kernel> cast(const geometry::Point_t<T>& p) {
	return {(typename Kernel::FT)(p.x()), (typename Kernel::FT)(p.y())};
}

template<typename t>
inline Point_t<t> operator+(const Point_t<t>& p1, const Point_t<t>& p2) {
	return {p1.x() + p2.x(), p1.y() + p2.y()};
}

template<typename kernel>
inline Point_t<kernel> operator*(const Point_t<kernel>& p, const typename kernel::FT v) {
	CGAL::Aff_transformation_2<kernel> scale(CGAL::SCALING, v);
	return std::move(scale(p));
}

template<typename kernel>
inline Point_t<kernel> operator*(const typename kernel::FT v, const Point_t<kernel>& p) {
	return std::move(p * v);
}


} //namespace
}
}
}

#endif
