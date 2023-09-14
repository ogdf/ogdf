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

#	include <ogdf/geometric/cr_min/geometry/objects/Direction.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Line.h>
#	include <ogdf/geometric/cr_min/geometry/objects/LineSegment.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Point.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Polygon.h>

#	include <boost/graph/breadth_first_search.hpp>
#	include <boost/graph/visitors.hpp>
#	include <tuple>

#	include <CGAL/Arr_extended_dcel.h>
#	include <CGAL/Arr_face_index_map.h>
#	include <CGAL/Arr_linear_traits_2.h>
#	include <CGAL/Arr_segment_traits_2.h>
#	include <CGAL/Arr_walk_along_line_point_location.h>
#	include <CGAL/Arrangement_2.h>
#	include <CGAL/Arrangement_with_history_2.h>
#	include <CGAL/graph_traits_dual_arrangement_2.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace geometry {

/*! Compute a planar subdivision of a given set of segments
 */
template<typename kernel>
class CGALPlanarSubdivision {
private:
	using LineSegment = geometry::LineSegment_t<kernel>;
	using Point = geometry::Point_t<kernel>;
	using Polygon = geometry::Polygon_t<kernel>;
	using Line = geometry::Line_t<kernel>;
	using Ray = geometry::Ray_t<kernel>;
	using Direction = geometry::Direction_t<kernel>;
	using Traits = CGAL::Arr_segment_traits_2<kernel>;

	using LinearTraits = CGAL::Arr_linear_traits_2<kernel>;

	using Dcel = CGAL::Arr_extended_dcel<Traits, unsigned int, unsigned int, unsigned int>;
	using ArrWH = CGAL::Arrangement_with_history_2<Traits, Dcel>;

	using Arr = CGAL::Arrangement_2<Traits>;
	using LinearArr = CGAL::Arrangement_2<LinearTraits>;

	using DualArr = CGAL::Dual<Arr>;
	using Face_index_map = CGAL::Arr_face_index_map<Arr>;

	template<typename _Arr>
	Polygon get_face(CGAL::Arr_walk_along_line_point_location<_Arr>& pl, const Point& q) {
		auto obj = pl.locate(q);

		using Vertex_const_handle = typename _Arr::Vertex_const_handle;
		using Halfedge_const_handle = typename _Arr::Halfedge_const_handle;
		using Face_const_handle = typename _Arr::Face_const_handle;

		const Vertex_const_handle* v;
		const Halfedge_const_handle* e;
		const Face_const_handle* f;

		Polygon poly;
		if ((f = boost::get<Face_const_handle>(&obj))) { // located inside a face
			OGDF_ASSERT(!(*f)->is_unbounded());
			auto c = (*f)->outer_ccb();
			do {
				poly.push_back(c->source()->point());
			} while (++c != (*f)->outer_ccb());
		} else if ((e = boost::get<Halfedge_const_handle>(&obj))) { // located on an edge

			poly.push_back((*e)->source()->point());
			poly.push_back((*e)->target()->point());
		} else if ((v = boost::get<Vertex_const_handle>(&obj))) { // located on a vertex
			poly.push_back((*v)->point());
		} else {
			// undefined behaviour
			OGDF_ASSERT(false);
		}

		return poly;
	}

public:
	using LCGEdge = std::tuple<unsigned int, unsigned int, unsigned int>;

	Polygon extract_face(const std::vector<Line>& lines, const Point& q) {
		LinearArr arr;
		CGAL::Arr_walk_along_line_point_location<LinearArr> pl(arr);
		CGAL::insert(arr, lines.begin(), lines.end());
		return get_face(pl, q);
	}

	Polygon extract_face(const std::vector<LineSegment>& segments, const Point& q) {
		Arr arr;
		CGAL::Arr_walk_along_line_point_location<Arr> pl(arr);
		CGAL::insert(arr, segments.begin(), segments.end());
		return get_face(pl, q);
	}

	/*! Compute planar subdivision of segments
	 * \param segments segments to intersect
	 * \param lines lines to insert between segments
	 * \param node_to_point results node to point mapping of node id to geometrical point
	 * \param edge_list resulting list of edges
	 * \param separator where the split happens
	 */
	void process(std::vector<LineSegment>& segments, std::vector<Line>& lines,
			std::vector<Point>& node_to_point, std::vector<LCGEdge>& edge_list, int separator = -1) {
		ArrWH arr;
		if (separator < 0) {
			separator = segments.size();
		}
		CGAL::insert(arr, segments.begin(), segments.begin() + separator);
		CGAL::insert(arr, lines.begin(), lines.end());
		//polygon edges
		CGAL::insert(arr, segments.begin() + separator, segments.end());

		unsigned int node = 1;
		for (auto v_it = arr.vertices_begin(); v_it != arr.vertices_end(); ++v_it) {
			v_it->set_data(node++);
			node_to_point.push_back(v_it->point());
		}

		auto f_make_edge = [&](decltype(arr.induced_edges_begin(arr.curves_begin()))& e_itr,
								   const Direction& ref, unsigned int flag) {
			auto e = *e_itr;
			LineSegment t(e->source()->point(), e->target()->point());

			if (t.direction() == ref) {
				return LCGEdge(e->source()->data() - 1, e->target()->data() - 1, flag);
			} else {
				return LCGEdge(e->target()->data() - 1, e->source()->data() - 1, flag);
			}
		};
		auto itr = arr.curves_begin();
		for (unsigned int i = 0; i < segments.size(); ++i, ++itr) {
			LineSegment s(itr->source(), itr->target());
			for (auto e_itr = arr.induced_edges_begin(itr); e_itr != arr.induced_edges_end(itr);
					++e_itr) {
				edge_list.push_back(f_make_edge(e_itr, s.direction(), i));
			}
		}

		for (unsigned int i = 0; i < lines.size(); ++i, ++itr) {
			OGDF_ASSERT(itr != arr.curves_end());
			for (auto e_itr = arr.induced_edges_begin(itr); e_itr != arr.induced_edges_end(itr);
					++e_itr) {
				auto e = *e_itr;
				if (e->source()->data() > 0 && e->target()->data() > 0) {
					edge_list.push_back(LCGEdge(e->source()->data() - 1, e->target()->data() - 1, i));
				}
			}
		}
	}

	/*! Compute planar subdivision of segments
	 * \param segments segments to intersect
	 * \param node_to_point results node to point mapping of node id to geometrical point
	 * \param edge_list resulting list of edges
	 * \param separator where the split happens
	 */

	void process(std::vector<LineSegment>& segments, std::vector<Point>& node_to_point,
			std::vector<LCGEdge>& edge_list, unsigned int separator = 0) {
		std::vector<Line> r;
		process(segments, r, node_to_point, edge_list, separator);
	}
};

}
}
}
}

#endif
