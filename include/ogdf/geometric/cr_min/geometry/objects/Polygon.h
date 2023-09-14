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
#	include <ogdf/geometric/cr_min/geometry/objects/Point.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Ray.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Rectangle.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Vector.h>
#	include <ogdf/geometric/cr_min/graph/GeometricDrawing.h>
#	include <ogdf/geometric/cr_min/graph/Path.h>
#	include <ogdf/geometric/cr_min/graph/PolylineDrawing.h>

#	include <limits>
#	include <list>
#	include <vector>

#	include <CGAL/Polygon_2.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace geometry {

using namespace tools;

template<typename kernel>
class Polygon_t : public CGAL::Polygon_2<kernel> {
private:
	using poly = CGAL::Polygon_2<kernel>;

public:
	using poly::poly; //inherit constructors

	Polygon_t(const poly& p) : poly(p) { }

	Polygon_t(const poly&& p) : poly(std::move(p)) { }

	Polygon_t() { }

	Polygon_t(const Polygon_t<kernel>& p) : poly(p) { }

	inline typename poly::Edge_const_iterator begin() const { return poly::edges_begin(); }

	inline typename poly::Edge_const_iterator end() const { return poly::edges_end(); };

	void push_back(const Point_t<kernel>& x) {
		if (poly::is_empty() || (poly::container().back() != x)) {
			poly::push_back(x);
		}
	}

	Polygon_t<kernel> operator=(const Polygon_t<kernel>& p) {
		(poly)(*this) = p;
		return *this;
	}
};

template<typename kernel, typename Graph>
inline Polygon_t<kernel> get_polygon(const graph::GeometricDrawing<kernel, Graph>& drawing,
		const graph::Path& path) { //TODO
	Polygon_t<kernel> polygon;
	for (auto& v : path.nodes()) {
		polygon.push_back(drawing.get_point(v));
	}
	return polygon;
}

template<typename kernel, typename Graph>
inline Polygon_t<kernel> get_polygon(const graph::PolylineDrawing<kernel, Graph>& drawing,
		const graph::Path& path) {
	Polygon_t<kernel> polygon;

	for (unsigned int i = 0; i < path.edges().size(); ++i) {
		auto e = path.edges()[i];

		if (!path.is_reversed(i)) {
			for (auto& p : drawing.get_polyline(e)) {
				polygon.push_back(p);
			}
		} else {
			auto& p = drawing.get_polyline(e);
			for (auto itr = p.rbegin(); itr != p.rend(); ++itr) {
				polygon.push_back(*itr);
			}
		}
	}
	return polygon;
}

template<typename kernel>
inline Polygon_t<kernel> get_polygon(const Rectangle_t<kernel>& rect) {
	Polygon_t<kernel> p;
	for (unsigned int i = 0; i < 4; ++i) {
		p.push_back(rect[i]);
	}
	return p;
}

template<typename kernel>
inline Polygon_t<kernel> get_polygon(const CGAL::Bbox_2& bb) {
	return get_polygon(Rectangle_t<kernel>(bb));
}

//TODO move somewhere usefull ;)
template<typename kernel>
LineSegment_t<kernel> clip(const Rectangle_t<kernel>& rect, const Line_t<kernel>& line) {
	//assume line intersects rect twice
	std::vector<Point_t<kernel>> is;
	Polygon_t<kernel> p = get_polygon(rect);
	for (auto s : p) {
		if (geometry::do_intersect_wo_target(line, s)) {
			is.push_back(geometry::intersect(line, s));
		}
	}

	return LineSegment_t<kernel>(is[0], is[1]);
}

template<typename kernel, typename Graph>
inline std::vector<typename Graph::Node> graph_from_polygon(const Polygon_t<kernel>& polygon,
		graph::GeometricDrawing<kernel, Graph>& drawing) {
	std::vector<typename Graph::Node> node_map;
	for (unsigned int i = 0; i < polygon.size(); ++i) {
		typename Graph::Node u = drawing.get_graph().add_node();
		drawing.set_point(u, polygon[i]);
		node_map.push_back(u);
	}

	for (unsigned int i = 0; i < polygon.size(); ++i) {
		drawing.get_graph().add_edge(node_map[i], node_map[(i + 1) % polygon.size()]);
	}

	return node_map;
}

template<typename kernel>
inline unsigned int next(const Polygon_t<kernel>& p, unsigned int i) {
	return (i + 1) % p.size();
}

template<typename kernel>
inline unsigned int prev(const Polygon_t<kernel>& p, unsigned int i) {
	return std::min((size_t)i - 1, p.size() - 1);
}

template<typename kernel>
inline Polygon_t<kernel> reverse(const Polygon_t<kernel>& polygon) {
	return std::move(Polygon_t<kernel>(polygon.container().rbegin(), polygon.container().rend()));
}

template<typename kernel>
inline bool contains(const Polygon_t<kernel>& polygon, const geometry::Point_t<kernel>& p) {
	OGDF_ASSERT(!is_clockwise(polygon));
	if (polygon.size() == 0) {
		return false;
	}
	if (polygon.size() == 1) {
		return polygon[0] == p;
	}
	if (CGAL::is_zero(CGAL::abs(polygon.area()))) {
		// all segments are nearly colinear, check if p is on one of them
		for (const auto& s : polygon) {
			if (is_on(s, p)) {
				return true;
			}
		}
		return false;
	}

	unsigned int segment = -1;
	geometry::LineSegment_t<kernel> l;
	do {
		++segment;
		l = {p, polygon[segment] + polygon.edge(segment).to_vector() * 0.5};
	} while (overlapping(l, polygon.edge(segment)));

	OGDF_ASSERT((size_t)segment < polygon.size());

	geometry::Ray_t<kernel> r(p, l.to_vector());
	unsigned int nof_intersections = 0;
	bool is_on_border = false;

	for (unsigned int i = 0; i < polygon.size(); ++i) {
		auto s = polygon.edge(i);

		if (geometry::do_intersect_wo_target(s, r)) {
			auto is = geometry::intersect(s, r);
			if (is == s.source()) {
				auto prev_seg = polygon.edge(prev(polygon, i));
				if (!geometry::right_turn(prev_seg, s)) { // not a concave corner
					++nof_intersections;
				}
			} else {
				++nof_intersections;
			}
		}
		is_on_border = is_on_border || is_on(s, p);
	}
	return nof_intersections % 2 == 1 || is_on_border;
}

//returns the distance of the point to the polygon. The Value is negative, if the point is not contained in the interior of the polygon.
template<typename Kernel>
typename Kernel::FT squared_distance(const Polygon_t<Kernel>& polygon, const Point_t<Kernel>& point) {
	typename Kernel::FT min_dist = CGAL::squared_distance(*polygon.begin(), point);
	for (const auto& e : polygon) {
		typename Kernel::FT sq = CGAL::squared_distance(e, point);
		if (sq < min_dist) {
			min_dist = sq;
		}
	}
	if (min_dist > 1e-12 && !contains(polygon, point)) {
		min_dist = -min_dist;
	}

	return min_dist;
}

//returns the distance of the point to the polygon. The Value is negative, if the point is not contained in the interior of the polygon.
template<typename Kernel>
Point_t<Kernel> centroid(const Polygon_t<Kernel>& polygon) {
	Point_t<Kernel> p(0, 0);
	for (const Point_t<Kernel>& q : polygon.container()) {
		p = p + q;
	}
	p = p * (1.0 / polygon.size());
	return p;
}

template<typename kernel>
inline Point_t<kernel> intersect(const Polygon_t<kernel>& polygon, const Ray_t<kernel>& ray,
		unsigned int& idx) {
	// FIXME implement me
	OGDF_ASSERT(false);
	(void)polygon;
	(void)ray;
	(void)idx;
}

template<typename kernel>
inline Point_t<kernel> intersect(const Polygon_t<kernel>& polygon, const Ray_t<kernel>& ray) {
	OGDF_ASSERT(contains(polygon, ray.source()));
	Point_t<kernel> p;
	for (auto s : polygon) {
		if (geometry::do_intersect_wo_target(s, ray)) {
			return geometry::intersect(s, ray);
		}
	}
	return {-1, -1};
}

template<typename kernel>
inline bool is_clockwise(const Polygon_t<kernel>& polygon) {
	return polygon.area() < 0;
}

template<typename kernel>
std::vector<unsigned int> find_duplicates(const Polygon_t<kernel>& polygon) {
	std::vector<unsigned int> segment_ids(polygon.size(), 0);
	for (unsigned int i = 0; i < polygon.size(); ++i) {
		segment_ids[i] = i;
	}
	std::sort(segment_ids.begin(), segment_ids.end(), [&](unsigned int i, unsigned int j) {
		return (polygon.edge(i).min() < polygon.edge(j).min())
				|| (polygon.edge(i).min() == polygon.edge(j).min()
						&& polygon.edge(i).max() < polygon.edge(j).max());
	});

	std::vector<unsigned int> is_duplicate(polygon.size(), -1);
	for (unsigned int i = 0; i + 1 < polygon.size(); ++i) {
		unsigned int a = segment_ids[i];
		unsigned int b = segment_ids[i + 1];
		if (polygon.edge(a).min() == polygon.edge(b).min()
				&& polygon.edge(a).max() == polygon.edge(b).max()) {
			is_duplicate[a] = b;
			is_duplicate[b] = a;
		}
	}
	return is_duplicate;
}

template<typename kernel>
std::string ggb(const Polygon_t<kernel>& polygon) {
	std::stringstream os;
	os << "polygon[";
	for (unsigned int i = 0; i < polygon.size(); ++i) {
		os << polygon[i];
		if (i + 1 < polygon.size()) {
			os << ",";
		}
	}
	os << "]";
	return os.str();
}

template<typename kernel>
std::ostream& operator<<(std::ostream& os, const Polygon_t<kernel>& p) {
	os << ggb(p);
	return os;
}


} //namespace
}
}
}

#endif
