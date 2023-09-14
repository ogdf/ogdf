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

#	include <ogdf/geometric/cr_min/datastructure/UnionFind.h>
#	include <ogdf/geometric/cr_min/geometry/algorithm/ExtractCellFromBloatedDual.h>
#	include <ogdf/geometric/cr_min/geometry/algorithm/RandomPointInPolygon.h>
#	include <ogdf/geometric/cr_min/graph/BloatedDual.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace graph {

template<typename Kernel>
class RandomPointsInCell {
private:
	using BD = BloatedDualGraph<Kernel>;

	using Point = geometry::Point_t<Kernel>;
	using Segment = geometry::LineSegment_t<Kernel>;
	using Polygon = geometry::Polygon_t<Kernel>;

	static std::vector<int> vertex_weights(const BD& bd, const Point& p,
			const std::vector<int>& left_to_right_cost, std::vector<bool>& visited,
			datastructure::UnionFind& uf) {
		typename BD::Node min_node = bd.segm_to_node_range[left_to_right_cost.size() - 4];

		auto rep = geometry::find_representative_node(bd, p);

		int min_weight = 0; // this value becomes negative, if we can improve...
		int max_weight = 0; // this value becomes negative, if we can improve...


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
				if (bd.is_face_edge(e)) {
					uf.merge(v, w);
				}
			}
		};

		while (!queue.empty()) {
			typename BD::Node current = queue.back();
			queue.pop_back();
			min_weight = std::min(min_weight, weight[current]);
			max_weight = std::max(max_weight, weight[current]);

			handle_edge(current, 3 * current);
			handle_edge(current, 3 * current + 1);
			handle_edge(current, 3 * current + 2);
		}

		for (auto& w : weight) {
			w = max_weight - w + 1;
		}

		return weight;
	}

public:
	static std::vector<Point> compute(const BD& bd, const std::vector<int>& left_to_right_cost,
			const Point& reference, unsigned int nof_points, std::mt19937_64& gen) {
		datastructure::UnionFind uf(bd.number_of_nodes());
		uf.all_to_singletons();

		std::vector<bool> visited(bd.number_of_nodes(), false);
		std::vector<int> weights = vertex_weights(bd, reference, left_to_right_cost, visited, uf);

		std::vector<unsigned int> repr;
		for (unsigned int u = 0; u < bd.number_of_nodes(); ++u) {
			if (uf.find(u) == u && visited[u]) {
				repr.push_back(u);
			}
		}


		std::vector<double> repr_weights(repr.size(), 0);
		std::vector<unsigned int> interval;

		for (unsigned int i = 0; i < repr.size(); ++i) {
			repr_weights[i] = std::pow(2, weights[repr[i]]);

			interval.push_back(i);
		}

		interval.push_back(repr.size());

		std::piecewise_constant_distribution<double> dist(interval.begin(), interval.end(),
				repr_weights.begin());

		std::vector<Point> point_set;
		if (!repr.empty()) {
			for (unsigned int r = 0; r < nof_points; ++r) {
				//select a region randomly inverse proportional to its number of crossings
				unsigned int v_id = std::floor(dist(gen));
				OGDF_ASSERT(v_id < repr.size());
				unsigned int v = repr[v_id];

				OGDF_ASSERT(v < bd.number_of_nodes());
				auto seq = geometry::extract_cell(bd, v);
				auto region = geometry::extract_polygon(bd.segments, seq);
				if (CGAL::abs(region.area()) > 1e-5) {
					auto p = geometry::random_point_in_polygon(region, gen);
					point_set.push_back(p);
				}
			}
		}
		return point_set;
	}
};


}
}
}
}

#endif
