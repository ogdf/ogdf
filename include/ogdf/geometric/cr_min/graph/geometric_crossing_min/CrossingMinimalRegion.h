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

#	include <ogdf/basic/Logger.h>
#	include <ogdf/basic/extended_graph_alg.h>
#	include <ogdf/geometric/cr_min/geometry/algorithm/BloatedDual.h>
#	include <ogdf/geometric/cr_min/geometry/algorithm/CGALPlanarSubdivision.h>
#	include <ogdf/geometric/cr_min/geometry/algorithm/ExtractCellFromBloatedDual.h>
#	include <ogdf/geometric/cr_min/geometry/algorithm/PlanarSubdivision.h>
#	include <ogdf/geometric/cr_min/graph/algorithms/Dijkstra.h>
#	include <ogdf/planarity/BoothLueker.h>

#	include <limits>

#	include <CGAL/Min_circle_2.h>
#	include <CGAL/Min_circle_2_traits_2.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace graph {

template<typename Kernel, typename Graph>
class CrossingMinimalRegion {
private:
	using Point = geometry::Point_t<Kernel>;
	using Segment = geometry::LineSegment_t<Kernel>;
	using Polygon = geometry::Polygon_t<Kernel>;
	using Drawing = graph::GeometricDrawing<Kernel, Graph>;

	using Node = typename Graph::Node;
	using Edge = typename Graph::Edge;

	using BD = BloatedDualGraph<Kernel>;

	static Segment get_segment(const Drawing& d, const Node& w, const Node& u,
			const geometry::Rectangle_t<Kernel>& rect) {
		OGDF_ASSERT(u != w);
		using FT = typename Kernel::FT;
		FT delta_x = (d.get_point(u).x() - d.get_point(w).x());
		FT delta_y = (d.get_point(u).y() - d.get_point(w).y());

		geometry::Ray_t<Kernel> r(d.get_point(u), geometry::Vector_t<Kernel>(delta_x, delta_y));
		double t = 0;
		for (unsigned int i = 0; i < 4; ++i) {
			t = std::max(t, CGAL::to_double(geometry::distance(rect.vertex(i), d.get_point(u))));
		}
		OGDF_ASSERT(t > 0);
		Segment s = {d.get_point(u), d.get_point(u) + geometry::normalize(r.to_vector()) * t * 1.1};

		return s;
	}

	static std::vector<int> generate_segments(const Drawing& d, const Node& v,
			const std::vector<Node>& neighbors, const std::vector<Edge>& edges, BD& bd,
			geometry::Rectangle_t<Kernel>& rect_box //limiting the area
	) {
		const auto& g = d.get_graph();

		bd.segments.reserve(g.number_of_nodes() * v->degree() + edges.size());
		std::vector<int> left_to_right_cost;

		std::set<Node> nodes;
		std::set<Edge> sampled_edges;
		for (auto e : edges) {
			auto s = d.get_segment(e);
			if (s.target() < s.source()) {
				s = {s.target(), s.source()};
			}

			OGDF_ASSERT(s.squared_length() > 0);
			bd.segments.push_back(s);
			nodes.insert(e->target());
			nodes.insert(e->source());
			sampled_edges.insert(e);

			left_to_right_cost.push_back(0);
			for (auto w : neighbors) {
				if (e->isIncident(w)) {
					continue;
				}
				if (geometry::left_turn(s, d.get_point(w))) {
					left_to_right_cost.back()++;
				} else {
					left_to_right_cost.back()--;
				}
			}
		}
		for (auto u : nodes) {
			for (auto w : neighbors) {
				if (u != w && u != v) {
					auto s = get_segment(d, w, u, rect_box);
					if (s.target() < s.source()) {
						s = {s.target(), s.source()};
					}
					OGDF_ASSERT(s.squared_length() > 0);
					bd.segments.push_back(s);

					left_to_right_cost.push_back(0);
					for (auto e : g.edges(u)) {
						if (sampled_edges.find(e) == sampled_edges.end()) {
							continue;
						}
						auto x = e->opposite(u);
						if (e->isIncident(w) || e->isIncident(v)) {
							continue;
						}
						if (geometry::left_turn(s, d.get_point(x))) {
							left_to_right_cost.back()--;
						} else {
							left_to_right_cost.back()++;
						}
					}
				}
			}
		}
		std::vector<Point> p_rect;
		for (unsigned int i = 0; i < 4; ++i) {
			double x = ogdf::randomNumber(0, INT_MAX) % 1000 / 1000.0;
			double y = ogdf::randomNumber(0, INT_MAX) % 1000 / 1000.0;
			x = y = 0;
			p_rect.push_back(rect_box.vertex(i) + geometry::Vector_t<Kernel>(x, y));
		}

		for (unsigned int i = 0; i < 4; ++i) {
			Segment s(p_rect[i], p_rect[(i + 1) % 4]);

			if (s.target() < s.source()) {
				s = {s.target(), s.source()};
			}

			OGDF_ASSERT(s.squared_length() > 0);
			bd.segments.push_back(s);
			if (geometry::left_turn(s, d.get_point(v))) {
				left_to_right_cost.push_back(
						d.get_graph().number_of_edges() * d.get_graph().number_of_edges());
			} else {
				left_to_right_cost.push_back(
						-d.get_graph().number_of_edges() * d.get_graph().number_of_edges());
			}
		}

		return left_to_right_cost;
	}

	static bool check_segments(std::vector<Segment>& segments) {
		bool valid = true;
		auto f_check_point = [&](const Segment& s, const Point& p) {
			if (s.has_on(p) && s.source() != p && s.target() != p) {
				// p is an interior point on s
				return false;
			} else {
				return true;
			}
		};
		for (unsigned int i = 0; i < segments.size(); ++i) {
			for (unsigned int j = i + 1; j < segments.size(); ++j) {
				if (!f_check_point(segments[i], segments[j].source())) {
					Logger::slout() << "[math] " << i << " " << j << "; " << segments[i] << " "
									<< segments[j] << std::endl;
					valid = false;
				}
				if (!f_check_point(segments[i], segments[j].target())) {
					Logger::slout() << "[math] " << i << " " << j << "; " << segments[i] << " "
									<< segments[j] << std::endl;
					valid = false;
				}
			}
		}
		return valid;
	}

	static typename BD::Node get_min_vertex(const BD& bd, const Point& p,
			const std::vector<int>& left_to_right_cost) {
		typename BD::Node min_node = bd.segm_to_node_range[left_to_right_cost.size()
				- 4]; // minimal node id of the nodes that represent the canvas / rect


		auto rep = geometry::find_representative_node(bd, p);

		typename BD::Node min_vertex = rep;
		double min_weight = 0; // this value becomes negative, if we can improve...

		std::vector<typename BD::Node> parent(bd.number_of_nodes());
		std::vector<typename BD::Edge> parent_edge(bd.number_of_nodes());

		auto f_weight = [&](typename BD::Node v, typename BD::Edge e) {
			if (bd.is_face_edge(e) || bd.degree(v) == 2) {
				return 0;
			} else if (bd.is_left(v)) {
				return left_to_right_cost[bd.node_to_segment[v]];
			} else {
				return -left_to_right_cost[bd.node_to_segment[v]];
			}
		};

		auto expand = [&](typename BD::Node w, typename BD::Edge e) {
			bool bad_edge = w >= min_node && bd.edges[e] >= min_node; //TODO correct??
			return !bad_edge;
		};

		parent[rep] = rep;

		std::vector<typename BD::Node> queue;
		std::vector<int> weight(bd.number_of_nodes(), 0);
		std::vector<bool> visited(bd.number_of_nodes(), false); //TODO use flags?

		queue.push_back(rep);
		weight[rep] = 0;
		visited[rep] = true;

		auto handle_edge = [&](typename BD::Node v, typename BD::Edge e) {
			typename BD::Node w = bd.edges[e];
			OGDF_ASSERT(!visited[w] || v == w || weight[w] == weight[v] + f_weight(v, e));
			if (!visited[w] && expand(v, e) && v != w) {
				weight[w] = weight[v] + f_weight(v, e);
				visited[w] = true;
				queue.push_back(w);
				parent[w] = v;
				parent_edge[w] = e;
			}
		};

		while (!queue.empty()) {
			typename BD::Node current = queue.back();
			queue.pop_back();
			if (weight[current] < min_weight) {
				min_vertex = current;
				min_weight = weight[current];
			}


			handle_edge(current, 3 * current);
			handle_edge(current, 3 * current + 1);
			handle_edge(current, 3 * current + 2);
		}
		return min_vertex;
	}

public:
	BD bd;
	std::vector<int> left_to_right_cost;

	void build_bd(const Drawing& d, const Node& v, const std::vector<Node>& neighbors,
			const std::vector<Edge>& edges, geometry::Rectangle_t<Kernel>& rect_box,
			unsigned int& arr_n_segs, unsigned int& nof_intersections_in_arr) {
		bd.clear();
		left_to_right_cost.clear();

		left_to_right_cost = generate_segments(d, v, neighbors, edges, bd, rect_box);
		OGDF_ASSERT(check_segments(bd.segments));

#	ifdef OGDF_GEOMETRIC_CR_MIN_DEBUG
		std::cout << "intersect..." << std::flush;
#	endif
		geometry::PlanarSubdivision<Kernel, false> ps;
		bd.intersecting_segments = std::move(ps.subdivision(bd.segments));
		arr_n_segs = bd.segments.size();
		nof_intersections_in_arr = bd.intersecting_segments.size();
#	ifdef OGDF_GEOMETRIC_CR_MIN_DEBUG
		std::cout << "done" << std::endl;

		std::cout << "nof segments in bd: " << arr_n_segs << std::endl;
		std::cout << "nof intersectsion in arr: " << nof_intersections_in_arr << std::endl;
		std::cout << "construct bd..." << std::flush;
#	endif
		static geometry::BloatedDual<Kernel> bdc;
		bdc.construct(bd);
#	ifdef OGDF_GEOMETRIC_CR_MIN_DEBUG
		std::cout << "done" << std::endl;
#	endif
	}

	Polygon compute(const Drawing& d, const Node& v, const std::vector<Node>& neighbors,
			const std::vector<Edge>& edges, geometry::Rectangle_t<Kernel>& rect_box,
			unsigned int& arr_n_segs, unsigned int& nof_intersections_in_arr) {
		build_bd(d, v, neighbors, edges, rect_box, arr_n_segs, nof_intersections_in_arr);

		auto dr = geometry::BloatedDual<Kernel>::compute_drawing(bd);
		auto min_vertex = get_min_vertex(bd, d.get_point(v), left_to_right_cost);
		auto seq = geometry::extract_cell(bd, min_vertex);
		auto opt_region = geometry::extract_polygon(bd.segments, seq);
		return opt_region;
	}
};

}
}
}
}

#endif
