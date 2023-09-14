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

#	ifndef OGDF_GEOMETRIC_INEXACT_NUMBER_TYPE
#		define OGDF_GEOMETRIC_INEXACT_NUMBER_TYPE false
#	endif

#	include <ogdf/geometric/cr_min/geometry/objects/Line.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Point.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Ray.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Vector.h>
#	include <ogdf/geometric/cr_min/tools/math.h>

#	include <iomanip>
#	include <limits>

#	include <CGAL/Segment_2.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace geometry {
using namespace tools;

template<typename kernel>
using LineSegment_t = CGAL::Segment_2<kernel>;

/* computes angle between this segment and @ls
	* @ls line segment
	* @return angle in [0, \Pi]
	* */
template<typename kernel>
inline typename kernel::FT angle(const LineSegment_t<kernel>& l1, const LineSegment_t<kernel>& l2) {
	OGDF_ASSERT(l1.to_vector().squared_length() > 0);
	OGDF_ASSERT(l2.to_vector().squared_length() > 0);
	return geometry::angle(l1.to_vector(), l2.to_vector());
}

/* compute the full angle between this segment and @ls  counterclockwise order
	* @ls line segment
	* @return angle in [0, 2 * \Pi]
	* */
template<typename kernel>
inline typename kernel::FT full_angle(const LineSegment_t<kernel>& l1,
		const LineSegment_t<kernel>& l2) {
	OGDF_ASSERT(l1.to_vector().squared_length() > 0);
	OGDF_ASSERT(l2.to_vector().squared_length() > 0);
	return full_angle(l1.to_vector(), l2.to_vector());
}

template<typename kernel>
inline bool turn(const LineSegment_t<kernel>& ls, const Line_t<kernel>& l) {
	return turn(ls.to_vector(), l.to_vector());
}

template<typename kernel>
inline bool left_turn(const LineSegment_t<kernel>& ls, const LineSegment_t<kernel>& l) {
	return left_turn(ls.to_vector(), l.to_vector());
}

template<typename kernel>
inline bool right_turn(const LineSegment_t<kernel>& ls, const LineSegment_t<kernel>& l) {
	return right_turn(ls.to_vector(), l.to_vector());
}

template<typename kernel>
inline bool turn(const LineSegment_t<kernel>& ls, const Point_t<kernel>& p) {
	return turn(ls.supporting_line(), p);
}

template<typename kernel>
inline bool left_turn(const LineSegment_t<kernel>& ls, const Point_t<kernel>& p) {
	return left_turn(ls.supporting_line(), p);
}

template<typename kernel>
inline bool right_turn(const LineSegment_t<kernel>& ls, const Point_t<kernel>& p) {
	return right_turn(ls.supporting_line(), p);
}

template<typename kernel>
inline typename kernel::FT squared_convex_parameter_directed(const LineSegment_t<kernel>& ls,
		const Point_t<kernel>& p) {
	OGDF_ASSERT(ls.has_on(p));
	return (ls.source() - p).squared_length() / ls.to_vector().squared_length();
}

template<typename kernel>
inline typename kernel::FT squared_convex_parameter(const LineSegment_t<kernel>& ls,
		const Point_t<kernel>& point) {
	OGDF_ASSERT(ls.has_on(point));
	return std::max(squared_convex_parameter_directed(ls, point),
			squared_convex_parameter_directed(ls.opposite(), point));
}

template<typename kernel>
inline Point_t<kernel> rel_point_at(const LineSegment_t<kernel>& ls, const typename kernel::FT l) {
	OGDF_ASSERT(0 <= l && l <= 1.0);
	return ls.source() + (ls.to_vector() * l);
}

template<typename kernel>
inline typename kernel::FT value_at(const LineSegment_t<kernel>& ls, const typename kernel::FT x) {
	if (ls.supporting_line().is_vertical()) {
		return std::numeric_limits<typename kernel::FT>::infinity();
	} else {
		return ls.supporting_line().y_at_x(x);
	}
}

template<typename kernel>
inline typename kernel::FT value_at(const LineSegment_t<kernel>& ls, const Point_t<kernel>& p) {
	if (ls.supporting_line().is_vertical()) {
		return p.y();
	} else {
		return ls.supporting_line().y_at_x(p.x());
	}
}

template<bool closed, typename kernel>
inline bool is_on(const LineSegment_t<kernel>& ls, const Point_t<kernel>& point) {
	if (!closed) {
		return ls.has_on(point) && ls.target() != point;
	} else {
		return ls.has_on(point);
	}
}

template<typename kernel>
inline bool is_on(const LineSegment_t<kernel>& ls, const Point_t<kernel>& point) {
	return is_on<false>(ls, point);
}

template<typename kernel>
inline bool overlapping(const LineSegment_t<kernel>& l1, const LineSegment_t<kernel>& l2) {
	const bool parallel = CGAL::parallel(l1.supporting_line(), l2.supporting_line());
	const bool on_l1 = is_on<true>(l1, l2.source()) || is_on<true>(l1, l2.target());
	const bool on_l2 = is_on<true>(l2, l1.source()) || is_on<true>(l2, l1.target());
	return (on_l1 || on_l2) && parallel;
}

template<typename kernel>
inline bool overlapping(const LineSegment_t<kernel>& ls, const Line_t<kernel>& line) {
	const bool parallel = parallel(ls.supporting_line(), line);
	// if both segments are parallel, then it would be sufficient to check wheather base or end is on the line
	// to be a little more robust, check if one of both is on the line. therefore small error are overseen
	const bool on_line = line.has_on(ls.source()) || line.has_on(ls.target());
	return on_line && parallel;
}

template<typename kernel>
inline bool do_intersect_wo_targets(const LineSegment_t<kernel>& l1, const LineSegment_t<kernel>& l2) {
	return CGAL::do_intersect(l1, l2) && !l1.has_on(l2.target()) && !l2.has_on(l1.target());
}

template<typename kernel>
inline bool do_intersect_open(const LineSegment_t<kernel>& l1, const LineSegment_t<kernel>& l2) {
	return CGAL::do_intersect(l1, l2) && !l1.has_on(l2.target()) && !l2.has_on(l1.target())
			&& !l1.has_on(l2.source()) && !l2.has_on(l1.source());
}

template<typename kernel>
inline bool do_intersect_wo_target(const LineSegment_t<kernel>& l1, const Line_t<kernel>& l) {
	return CGAL::do_intersect(l1, l) && !l.has_on(l1.target());
}

template<typename kernel>
inline bool do_intersect_wo_target(const Line_t<kernel>& l, const LineSegment_t<kernel>& l1) {
	return geometry::do_intersect_wo_target(l1, l);
}

template<typename kernel>
inline bool do_intersect_wo_target(const LineSegment_t<kernel>& l1, const Ray_t<kernel>& r) {
	return CGAL::do_intersect(l1, r) && !r.has_on(l1.target());
}

template<typename kernel>
inline bool do_intersect_wo_target(const Ray_t<kernel>& r, const LineSegment_t<kernel>& l1) {
	return geometry::do_intersect_wo_target(l1, r);
}

template<typename kernel>
inline bool have_common_endpoints(const LineSegment_t<kernel>& l1, const LineSegment_t<kernel>& l2) {
	return l1.source() == l2.source() || l1.source() == l2.target() || l1.target() == l1.source()
			|| l1.target() == l2.target();
}

template<typename kernel>
inline Point_t<kernel> get_common_endpoint(const LineSegment_t<kernel>& l1,
		const LineSegment_t<kernel>& l2) {
	OGDF_ASSERT(have_common_endpoints(l1, l2));

	if (l1.source() == l2.source() || l1.source() == l2.target()) {
		return l1.source();
	} else {
		return l1.target();
	}
}

template<typename kernel>
inline Point_t<kernel> intersect(const LineSegment_t<kernel>& l1, const LineSegment_t<kernel>& l2) {
	OGDF_ASSERT(!l1.is_degenerate());
	OGDF_ASSERT(!l2.is_degenerate());

	auto result = CGAL::intersection(l1, l2);
	if (!result) {
		return Point_t<kernel>(std::numeric_limits<double>::max(),
				std::numeric_limits<double>::max());
	}

	Point_t<kernel> intersection;
	if (boost::get<LineSegment_t<kernel>>(&*result)) {
		auto s = boost::get<LineSegment_t<kernel>>(&*result);
		intersection = s->source();
	} else if (boost::get<Point_t<kernel>>(&*result)) {
		intersection = *boost::get<Point_t<kernel>>(&*result);
	}
	//TODO
	//OGDF_ASSERT(is_on<true>(l1, intersection));
	//OGDF_ASSERT(is_on<true>(l2, intersection));
	return intersection;
}

template<typename kernel>
inline Point_t<kernel> intersect(const LineSegment_t<kernel>& l1, const Line_t<kernel>& l2) {
	OGDF_ASSERT(!l1.is_degenerate());
	OGDF_ASSERT(!l2.is_degenerate());
	auto result = CGAL::intersection(l1, l2);
	if (!result) {
		return Point_t<kernel>(std::numeric_limits<unsigned int>::max(),
				std::numeric_limits<unsigned int>::max());
	}

	Point_t<kernel> intersection;
	if (boost::get<LineSegment_t<kernel>>(&*result)) {
		auto s = boost::get<LineSegment_t<kernel>>(&*result);
		intersection = s->source();
	} else {
		intersection = *boost::get<Point_t<kernel>>(&*result);
	}
	OGDF_ASSERT(is_on<true>(l1, intersection));
	OGDF_ASSERT(l2.has_on(intersection));
	return intersection;
}

template<typename kernel>
inline Point_t<kernel> intersect(const Line_t<kernel>& l1, const LineSegment_t<kernel>& l2) {
	return intersect(l2, l1);
}

template<typename kernel>
inline Point_t<kernel> intersect(const LineSegment_t<kernel>& l1, const Ray_t<kernel>& r) {
	OGDF_ASSERT(!l1.is_degenerate());
	OGDF_ASSERT(!r.is_degenerate());

	auto result = CGAL::intersection(l1, r);
	if (!result) {
		return Point_t<kernel>(std::numeric_limits<unsigned int>::max(),
				std::numeric_limits<unsigned int>::max());
	}


	Point_t<kernel> intersection;
	if (boost::get<LineSegment_t<kernel>>(&*result)) {
		auto s = boost::get<LineSegment_t<kernel>>(&*result);
		intersection = s->source();
	} else {
		intersection = *boost::get<Point_t<kernel>>(&*result);
	}
	OGDF_ASSERT(OGDF_GEOMETRIC_INEXACT_NUMBER_TYPE || is_on<true>(l1, intersection));
	OGDF_ASSERT(OGDF_GEOMETRIC_INEXACT_NUMBER_TYPE || r.has_on(intersection));
	return intersection;
}

template<typename kernel>
inline Point_t<kernel> intersect(const Ray_t<kernel>& l1, const LineSegment_t<kernel>& l2) {
	return intersect(l2, l1);
}

template<typename kernel>
inline typename kernel::FT order_along_seg(const LineSegment_t<kernel>& reference

		,
		const LineSegment_t<kernel>& l_1, const LineSegment_t<kernel>& l_2) {
	auto cross = [&](const Vector_t<kernel>& v_1, const Vector_t<kernel>& v_2) {
		return v_1.x() * v_2.y() - v_1.y() * v_2.x();
	};

	auto s = reference.target() - reference.source();
	auto r_1 = l_1.target() - l_1.source();
	auto r_2 = l_2.target() - l_2.source();

	auto c_1 = cross(reference.source() - l_1.source(), s);
	auto c_2 = cross(reference.source() - l_2.source(), s);

	auto d_1 = cross(r_1, s);
	auto d_2 = cross(r_2, s);

	return c_1 * d_2 - c_2 * d_1;
}

template<typename kernel>
void serialize(const LineSegment_t<kernel>& ls, std::ostream& os) {
	serialize(ls.source(), os);
	serialize(ls.target(), os);
}

template<typename kernel>
inline bool operator<(const LineSegment_t<kernel>& l1, const LineSegment_t<kernel>& l2) {
	return l1.source() < l2.source() || (l1.source() == l2.source() && l1.target() < l2.target());
}

template<typename Kernel, typename T>
LineSegment_t<Kernel> cast(const LineSegment_t<T>& segment) {
	return {cast<Kernel>(segment.source()), cast<Kernel>(segment.target())};
}

} //namespace
}
}
}

#endif
