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

#	include <ogdf/geometric/cr_min/geometry/objects/LineSegment.h>
#	include <ogdf/geometric/cr_min/tools/math.h>

#	include <functional>
#	include <iomanip>
#	include <vector>

#	include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#	include <CGAL/Boolean_set_operations_2.h>
#	include <CGAL/Circle_2.h>
#	include <CGAL/Circular_kernel_2.h>
#	include <CGAL/General_polygon_set_2.h>
#	include <CGAL/Gps_circle_segment_traits_2.h>

namespace ogdf {
namespace internal {
namespace gcm {

namespace geometry {
using namespace tools;

template<typename kernel>
using Circle_t = CGAL::Circle_2<kernel>;

template<typename kernel>
typename CGAL::Circular_kernel_2<kernel, CGAL::Algebraic_kernel_for_circles_2_2<typename kernel::FT>>::Point_2
convert(const typename kernel::Point_2& p) {
	return {p.x(), p.y()};
}

template<typename kernel>
typename CGAL::Circular_kernel_2<kernel, CGAL::Algebraic_kernel_for_circles_2_2<typename kernel::FT>>::Circle_2
convert(const Circle_t<kernel>& circle) {
	return {{circle.center().x(), circle.center().y()}, circle.squared_radius()};
}

template<typename kernel>
typename CGAL::Circular_kernel_2<kernel, CGAL::Algebraic_kernel_for_circles_2_2<typename kernel::FT>>::Line_2
convert(const geometry::Line_t<kernel>& line) {
	return {line.a(), line.b(), line.c()};
}

template<typename kernel>
typename CGAL::Circular_kernel_2<kernel, CGAL::Algebraic_kernel_for_circles_2_2<typename kernel::FT>>::Segment_2
convert(const geometry::LineSegment_t<kernel>& ls) {
	return {convert<kernel>(ls.source()), convert<kernel>(ls.target())};
}

template<typename kernel>
bool do_intersect(const Circle_t<kernel>& circle, const LineSegment_t<kernel>& ls) {
	if (circle.has_on_boundary(ls.source())) {
		return true;
	} else if (circle.has_on_unbounded_side(ls.source()) && circle.has_on_bounded_side(ls.target())) {
		return true;
	} else if (circle.has_on_unbounded_side(ls.target()) && circle.has_on_bounded_side(ls.source())) {
		return true;
	} else {
		return false;
	}
}

template<typename kernel>
void intersect(const Circle_t<kernel>& circle, const LineSegment_t<kernel>& ls,
		Point_t<kernel>& intersection_1, Point_t<kernel>& intersection_2) {
	using AKernel = CGAL::Algebraic_kernel_for_circles_2_2<typename kernel::FT>;
	using CKernel = CGAL::Circular_kernel_2<kernel, AKernel>;
	using CCircle = typename CKernel::Circle_2;
	using CLine = typename CKernel::Line_2;


	auto c = convert(circle);
	auto l = convert(ls);
	typename CKernel::Line_arc_2 line_arc(l);

	OGDF_ASSERT(geometry::do_intersect(c, l));
	using IS_Result = typename CGAL::CK2_Intersection_traits<CKernel, CCircle, CLine>::type;

	std::vector<IS_Result> is;
	CGAL::intersection(c, line_arc, std::back_inserter(is));

	using R = std::pair<CGAL::Circular_arc_point_2<CKernel>, unsigned int>;

	if (is.size() > 0) {
		auto is_1 = boost::get<R>(&is[0]);
		intersection_1 = {CGAL::to_double(is_1->first.x()), CGAL::to_double(is_1->first.y())};
	}
	if (is.size() > 1) {
		auto is_2 = boost::get<R>(&is[1]);
		intersection_2 = {CGAL::to_double(is_2->first.x()), CGAL::to_double(is_2->first.y())};
	}
}

template<typename kernel>
void intersect(const Circle_t<kernel>& c1, const Circle_t<kernel>& c2,
		Point_t<kernel>& intersection_1, Point_t<kernel>& intersection_2) {
	using AKernel = CGAL::Algebraic_kernel_for_circles_2_2<typename kernel::FT>;
	using CKernel = CGAL::Circular_kernel_2<kernel, AKernel>;
	using CCircle = typename CKernel::Circle_2;

	auto cc1 = convert(c1);
	auto cc2 = convert(c2);
	OGDF_ASSERT(CGAL::do_intersect(c1, c2));

	using IS_Result = typename CGAL::CK2_Intersection_traits<CKernel, CCircle, CCircle>::type;
	using R = std::pair<CGAL::Circular_arc_point_2<CKernel>, unsigned int>;

	std::vector<IS_Result> is;
	CGAL::intersection(cc1, cc2, std::back_inserter(is));

	if (is.size() > 0) {
		auto is_1 = boost::get<R>(&is[0]);
		if (is_1) {
			intersection_1 = {CGAL::to_double(is_1->first.x()), CGAL::to_double(is_1->first.y())};
		}
	}

	if (is.size() > 1) {
		auto is_2 = boost::get<R>(&is[1]);
		if (is_2) {
			intersection_2 = {CGAL::to_double(is_2->first.x()), CGAL::to_double(is_2->first.y())};
		}
	}
}

template<typename kernel>
bool contains(const Circle_t<kernel>& c, const Point_t<kernel>& p) {
	return c.has_on_boundary(p) || c.has_on_bounded_side(p);
}

template<typename kernel>
class CircleOperation_t {
private:
	using Point = Point_t<kernel>;
	using LineSegment = LineSegment_t<kernel>;
	using Circle = Circle_t<kernel>;

protected:
	template<typename F>
	void intersect_pair(const Circle& c_1, const Circle& c_2, F&& f) const {
		if (CGAL::do_intersect(c_1, c_2)) {
			Point is_1, is_2;
			geometry::intersect(c_1, c_2, is_1, is_2);
			f(is_1);
			if (is_1 != is_2) {
				f(is_2);
			}
		}
	}

public:
	const Circle& circle_1;
	const Circle& circle_2;

	CircleOperation_t(const Circle& circle1, const Circle& circle2)
		: circle_1(circle1), circle_2(circle2) {
		/*nothing to do*/
	}

	virtual ~CircleOperation_t() { /*nothing to do*/
	}

	virtual void intersect(const LineSegment&, std::function<void(const Point)>&&) const = 0;
	virtual void intersect(const CircleOperation_t<kernel>&,
			std::function<void(const Point)>&&) const = 0;
	virtual bool do_intersect(const LineSegment_t<kernel>&) const = 0;
	virtual bool contains(const Point_t<kernel>&) const = 0;
	virtual bool contains(const Circle& other, const Point& p) const = 0;

	friend std::ostream& operator<<(std::ostream& os, const CircleOperation_t<kernel>& c) {
		os << c.circle_1 << " " << c.circle_2;
		return os;
	}
};

template<typename kernel>
class IntersectingCircles_t : public CircleOperation_t<kernel> {
private:
	using CircleOperation = CircleOperation_t<kernel>;
	using parent = CircleOperation;
	using ICOp = IntersectingCircles_t<kernel>;
	using Circle = Circle_t<kernel>;
	using LineSegment = LineSegment_t<kernel>;
	using Point = Point_t<kernel>;

	template<typename F>
	void intersect_circle(const LineSegment& ls, const Circle& c, const Circle& other, F&& f) const {
		if (geometry::do_intersect(c, ls)) {
			Point is_1, is_2;
			geometry::intersect(c, ls, is_1, is_2);
			//intersection is definetyl in c
			if (contains(other, is_1)) {
				f(is_1);
			}
			if (is_1 != is_2 && contains(other, is_2) && is_on(ls, is_2)) {
				f(is_2);
			}
		}
	}

public:
	bool contains(const Circle& other, const Point& p) const {
		return geometry::contains(other, p);
	}

	IntersectingCircles_t(const Circle& circle_1, const Circle& circle_2)
		: parent(circle_1, circle_2) {
		/*nothing to do*/
	}

	void intersect(const CircleOperation& op, std::function<void(const Point)>&& f) const {
		parent::intersect_pair(parent::circle_1, op.circle_1, [&](const Point& p) {
			if (contains(parent::circle_2, p) && op.contains(op.circle_2, p)) {
				f(p);
			}
		});

		parent::intersect_pair(parent::circle_1, op.circle_2, [&](const Point& p) {
			if (contains(parent::circle_2, p) && op.contains(op.circle_1, p)) {
				f(p);
			}
		});

		parent::intersect_pair(parent::circle_2, op.circle_1, [&](const Point& p) {
			if (contains(parent::circle_1, p) && op.contains(op.circle_2, p)) {
				f(p);
			}
		});

		parent::intersect_pair(parent::circle_2, op.circle_2, [&](const Point& p) {
			if (contains(parent::circle_1, p) && op.contains(op.circle_1, p)) {
				f(p);
			}
		});

		parent::intersect_pair(parent::circle_1, parent::circle_2, [&](const Point& p) {
			if (op.contains(p)) {
				f(p);
			}
		});

		parent::intersect_pair(op.circle_1, op.circle_2, [&](const Point& p) {
			if (contains(p)) {
				f(p);
			}
		});

		if (op.contains(parent::circle_1.center()) && contains(parent::circle_1.center())) {
			f(parent::circle_1.center());
		}
		if (op.contains(op.circle_1.center()) && contains(op.circle_1.center())) {
			f(op.circle_1.center());
		}
	}

	void intersect(const LineSegment& ls, std::function<void(const Point)>&& f) const {
		intersect_circle(ls, parent::circle_1, parent::circle_2, f);
		intersect_circle(ls, parent::circle_2, parent::circle_1, f);
	}

	bool do_intersect(const LineSegment& ls) const {
		return geometry::do_intersect(parent::circle_1, ls)
				&& geometry::do_intersect(parent::circle_2, ls);
	}

	bool contains(const Point& p) const {
		return geometry::contains(parent::circle_1, p) && geometry::contains(parent::circle_2, p);
	}
};

template<typename kernel>
class CombinedCircles_t : public CircleOperation_t<kernel> {
private:
	using CircleOperation = CircleOperation_t<kernel>;
	using parent = CircleOperation;
	using CCOp = CombinedCircles_t<kernel>;
	using Circle = Circle_t<kernel>;
	using LineSegment = LineSegment_t<kernel>;
	using Point = Point_t<kernel>;

	template<typename F>
	void intersect_circle(const LineSegment& ls, const Circle& c, F&& f) const {
		if (geometry::do_intersect(c, ls)) {
			Point is_1, is_2;
			geometry::intersect(c, ls, is_1, is_2);
			f(is_1);
			if (is_1 != is_2 && is_on(ls, is_2)) {
				f(is_2);
			}
		}
	}

public:
	CombinedCircles_t(const Circle& circle_1, const Circle& circle_2) : parent(circle_1, circle_2) {
		/*nothing to do*/
	}

	void intersect(const LineSegment& ls, std::function<void(const Point)>&& f) const {
		intersect_circle(ls, parent::circle_1, f);
		intersect_circle(ls, parent::circle_2, f);
	}

	// TODO assumes circle operation is of same type....
	void intersect(const CircleOperation& op, std::function<void(const Point)>&& f) const {
		parent::intersect_pair(parent::circle_1, op.circle_1, f);
		parent::intersect_pair(parent::circle_1, op.circle_2, f);
		parent::intersect_pair(parent::circle_2, op.circle_1, f);
		parent::intersect_pair(parent::circle_2, op.circle_2, f);

		parent::intersect_pair(parent::circle_1, parent::circle_2, [&](const Point& p) {
			if (op.contains(p)) {
				f(p);
			}
		});
		parent::intersect_pair(op.circle_1, op.circle_2, [&](const Point& p) {
			if (contains(p)) {
				f(p);
			}
		});

		if (op.contains(parent::circle_1.center())) {
			f(parent::circle_1.center());
		}

		if (op.contains(parent::circle_2.center())) {
			f(parent::circle_1.center());
		}

		if (contains(op.circle_1.center())) {
			f(op.circle_1.center());
		}

		if (contains(op.circle_2.center())) {
			f(op.circle_2.center());
		}
	}

	bool do_intersect(const LineSegment& ls) const {
		return geometry::do_intersect(parent::circle_1, ls)
				|| geometry::do_intersect(parent::circle_2, ls);
	}

	bool contains(const Circle&, const Point&) const { return true; }

	bool contains(const Point& p) const {
		return geometry::contains(parent::circle_1, p) || geometry::contains(parent::circle_2, p);
	}
};

template<typename Kernel>
typename CGAL::Gps_circle_segment_traits_2<Kernel>::General_polygon_2 construct_polygon_from_circle(
		const Circle_t<Kernel>& circle) {
	using Traits_2 = CGAL::Gps_circle_segment_traits_2<Kernel>;
	using Polygon_2 = typename Traits_2::General_polygon_2;
	using Curve_2 = typename Traits_2::Curve_2;
	using X_monotone_curve_2 = typename Traits_2::X_monotone_curve_2;

	//Subdivide the circle into two x-monotone arcs.
	Traits_2 traits;
	Curve_2 curve(circle);
	std::list<CGAL::Object> objects;
	traits.make_x_monotone_2_object()(curve, std::back_inserter(objects));
	OGDF_ASSERT(objects.size() == 2);
	// Construct the polygon.
	Polygon_2 pgn;
	X_monotone_curve_2 arc;
	std::list<CGAL::Object>::iterator iter;
	for (iter = objects.begin(); iter != objects.end(); ++iter) {
		CGAL::assign(arc, *iter);
		pgn.push_back(arc);
	}
	return pgn;
}

} //namespace
}
}
}

#endif
