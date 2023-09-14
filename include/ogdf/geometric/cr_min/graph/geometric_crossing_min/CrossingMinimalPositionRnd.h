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

#	include <ogdf/geometric/cr_min/geometry/algorithm/CollinearTriple.h>
#	include <ogdf/geometric/cr_min/geometry/algorithm/RandomPointInPolygon.h>
#	include <ogdf/geometric/cr_min/graph/geometric_crossing_min/CrossingMinimalRegion.h>
#	include <ogdf/geometric/cr_min/graph/geometric_crossing_min/RandomPoint.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace graph {

template<typename Kernel, typename Graph, typename CRegion = CrossingMinimalRegion<Kernel, Graph>>
class CrossingMinimalPositionRnd {
private:
	using Point = geometry::Point_t<Kernel>;
	using Segment = geometry::LineSegment_t<Kernel>;
	using Polygon = geometry::Polygon_t<Kernel>;
	using Drawing = graph::GeometricDrawing<Kernel, Graph>;
	using Node = typename Graph::Node;
	using Edge = typename Graph::Edge;

public:
	static Point compute(const Drawing& d, const Node& v, const std::vector<Edge>& sample,
			const unsigned int random_draws, const unsigned int neighbor_threshold,
			geometry::Rectangle_t<Kernel>& rect_box,
			bool use_min_region //otherwise select with probability to cr. number
			,
			std::mt19937_64& rnd) {
		unsigned int dump_a, dump_b;

		std::vector<Node> neighbors;
		for (Node u : d.get_graph().neighbors(v)) {
			neighbors.push_back(u);
		}

		std::shuffle(neighbors.begin(), neighbors.end(), rnd);

		std::vector<Node> current;


		Point p = d.get_point(v);
		auto count = [&](Point p_v) {
			unsigned int n_cr = 0;
			for (auto edge : d.get_graph().edges()) {
				if (!edge->isIncident(v)) {
					Segment s = d.get_segment(edge);
					for (auto u : d.get_graph().neighbors(v)) {
						if (!edge->isIncident(u)) {
							Segment s_f(p_v, d.get_point(u));
							n_cr += CGAL::do_intersect(s, s_f);
						}
					}
				}
			}
			return n_cr;
		};

		unsigned int min_cr = count(p);

		for (size_t i = 0; i < neighbors.size(); i = i + neighbor_threshold) {
			current.clear();
			std::copy(neighbors.begin() + i,
					neighbors.begin() + std::min(i + neighbor_threshold, neighbors.size()),
					std::back_inserter(current));
			if (use_min_region) {
				CRegion cr;

				auto region = cr.compute(d, v, current, sample, rect_box, dump_a, dump_b);

				OGDF_ASSERT(!CGAL::is_zero(CGAL::abs(region.area())));

				if (geometry::is_clockwise(region)) {
					region = geometry::reverse(region);
				}

				if (region.size() > 2) {
					for (unsigned int x = 0; x < random_draws; ++x) {
						Point p_v = geometry::random_point_in_polygon(region, rnd);
						unsigned int n_cr = count(p_v);
						if (n_cr < min_cr) {
							min_cr = n_cr;
							p = p_v;
						}
					}
				}
			} else {
				CRegion cr;

				cr.build_bd(d, v, current, sample, rect_box, dump_a, dump_b);

				auto rnd_points = RandomPointsInCell<Kernel>::compute(cr.bd, cr.left_to_right_cost,
						d.get_point(v), random_draws, rnd);

				for (auto p_v : rnd_points) {
					unsigned int n_cr = count(p_v);
					if (n_cr < min_cr) {
						min_cr = n_cr;
						p = p_v;
					}
				}
			}
		}

		//round
		OGDF_ASSERT(rect_box.has_on_bounded_side(p));
		p = Point(CGAL::to_double(p.x()), CGAL::to_double(p.y()));
		OGDF_ASSERT(rect_box.has_on_bounded_side(p));
		return p;
	}
};

}
}
}
}

#endif
