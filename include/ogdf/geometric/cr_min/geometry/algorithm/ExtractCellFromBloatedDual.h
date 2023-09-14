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

#	include <ogdf/basic/basic.h>
#	include <ogdf/geometric/cr_min/graph/BloatedDual.h>

#	include <vector>

namespace ogdf {
namespace internal {
namespace gcm {
namespace geometry {

namespace details {


template<typename Kernel>
typename graph::BloatedDualGraph<Kernel>::Edge forward(
		const typename graph::BloatedDualGraph<Kernel>::Node prev,
		const typename graph::BloatedDualGraph<Kernel>::Node cur,
		const typename graph::BloatedDualGraph<Kernel>& bd) {
	OGDF_ASSERT(bd.degree(cur) >= 2);

	auto is_good = [&](const unsigned int e) {
		//return bd.edge_to_segments[e].seg_a != bd.edge_to_segments[e].seg_b
		return bd.is_face_edge(e) && bd.edges[e] != prev && bd.edges[e] != cur; // not a self-loop
	};

	using Edge = typename graph::BloatedDualGraph<Kernel>::Edge;
	const Edge e_1 = 3 * cur;
	const Edge e_2 = 3 * cur + 1;
	const Edge e_3 = 3 * cur + 2;

	if (is_good(e_1)) {
		return e_1;
	} else if (is_good(e_2)) {
		return e_2;
	} else {
		return e_3;
	}
}

template<typename Kernel>
typename graph::BloatedDualGraph<Kernel>::Edge find_start(
		const typename graph::BloatedDualGraph<Kernel>::Node cur,
		const graph::BloatedDualGraph<Kernel>& bd) {
	using Edge = typename graph::BloatedDualGraph<Kernel>::Edge;
	const Edge e_1 = 3 * cur;
	const Edge e_2 = 3 * cur + 1;
	const Edge e_3 = 3 * cur + 2;

	if (bd.is_face_edge(e_1) || bd.degree(cur) == 1) {
		return e_1;
	} else if (bd.is_face_edge(e_2)) {
		return e_2;
	} else {
		return e_3;
	}
}

}

/**
 * Given a vertex @p v of a bloated dual @p dual, this procedure returns an ordered list of segments
 *  such that the intersection of the respective half planes represent v
 * @warning does not handle regions with disconnected boundaries properly
 **/
//use clockwise and counter clockwise to decode which side??
template<typename Kernel>
std::vector<unsigned int> extract_cell(const graph::BloatedDualGraph<Kernel>& bd,
		const typename graph::BloatedDualGraph<Kernel>::Node v) {
	using Graph = graph::BloatedDualGraph<Kernel>;
	using Node = typename Graph::Node;

	std::vector<unsigned int> segments;

	auto cur_edge = details::find_start(v, bd);

	//in case that bd.degree(v) == 1, the arrangment is a single segment
	if (bd.degree(v) > 1) {
		Node prev = v;
		segments.push_back(bd.node_to_segment[v]);
		Node current = bd.edges[cur_edge];
		while (current != v) {
			cur_edge = details::forward(prev, current, bd);

			//segment_pairs.push_back(bd.edge_to_segments[cur_edge]);
			segments.push_back(bd.node_to_segment[current]);
			prev = current;
			current = bd.edges[cur_edge];
		}
	}
	return segments;
}

/*
 *@caution does not handle open regions properly
 */
template<typename Kernel>
Polygon_t<Kernel> extract_polygon(const std::vector<LineSegment_t<Kernel>>& segments,
		const std::vector<unsigned int>& seq) {
	Polygon_t<Kernel> poly;

	//for (auto x : seq) {
	for (unsigned int i = 0; i < seq.size(); ++i) {
		Point_t<Kernel> p;
		bool common_endpoint = false;
		unsigned int seg_a = seq[i];

		unsigned int seg_b = -1;
		if (i + 1 < seq.size()) {
			seg_b = seq[i + 1];
		} else {
			seg_b = seq[0];
		}

		for (unsigned int l = 0; l <= 1; ++l) {
			for (unsigned int j = 0; j <= 1; ++j) {
				OGDF_ASSERT(seg_a < segments.size());
				OGDF_ASSERT(seg_b < segments.size());
				if (segments[seg_a].vertex(l) == segments[seg_b].vertex(j)) {
					p = segments[seg_a].vertex(l);
					common_endpoint = true;
				}
			}
		}

		if (!common_endpoint) {
			p = geometry::intersect(segments[seg_a], segments[seg_b]);
		}

		if (poly.is_empty() || *(--poly.vertices_end()) != p) {
			poly.push_back(p);
		}
	}

	if (*poly.vertices_begin() == *(--poly.vertices_end())) {
		OGDF_ASSERT(poly.size() > 0);
		poly.erase(--poly.vertices_end());
	}

	return poly;
}

template<typename Kernel>
typename graph::BloatedDualGraph<Kernel>::Node find_representative_node(
		const graph::BloatedDualGraph<Kernel>& bd, const Point_t<Kernel>& p) {
	auto& segments = bd.segments;
	unsigned int opt = -1; // segment with the smallest distance to p
	Point_t<Kernel> is_with_r;
	while (opt > segments.size()) {
		Ray_t<Kernel> r(p,
				Vector_t<Kernel>(ogdf::randomNumber(0, INT_MAX) % 1000 - 500,
						ogdf::randomNumber(0, INT_MAX) % 1000 - 500));
		for (unsigned i = 0; i < segments.size(); ++i) {
			if (CGAL::do_intersect(segments[i], r)) {
				if (opt > segments.size()) {
					opt = i;
					is_with_r = geometry::intersect(segments[i], r);
				} else {
					auto o = geometry::intersect(segments[i], r);
					if (squared_distance(o, p) < squared_distance(is_with_r, p)) {
						opt = i;
						is_with_r = o;
					}
				}
			}
		}
	}

	unsigned int opt_is = 0;
	typename Kernel::FT distance = CGAL::squared_distance(segments[opt].source(), is_with_r);
	// find subsegment on segments[opt]
	for (unsigned int i = 0; i < bd.segment_to_intersections[opt].size(); ++i) {
		auto is = bd.get_intersection_point(bd.segment_to_intersections[opt][i]);
		if (CGAL::squared_distance(is, is_with_r) < distance) {
			distance = CGAL::squared_distance(is, is_with_r);
			opt_is = i;
		}
	}

	typename graph::BloatedDualGraph<Kernel>::Node v;
	bool left = geometry::left_turn(segments[opt], p);

	// endpoint should be an intersection
	OGDF_ASSERT(!bd.segment_to_intersections[opt].empty());

	if (CGAL::squared_distance(segments[opt].source(), is_with_r)
			< CGAL::squared_distance(segments[opt].source(),
					bd.get_intersection_point(bd.segment_to_intersections[opt][opt_is]))) {
		if (left) {
			v = bd.first_left(opt, bd.segment_to_intersections[opt][opt_is]);
		} else {
			v = bd.first_right(opt, bd.segment_to_intersections[opt][opt_is]);
		}

	} else {
		if (left) {
			v = bd.second_left(opt, bd.segment_to_intersections[opt][opt_is]);
		} else {
			v = bd.second_right(opt, bd.segment_to_intersections[opt][opt_is]);
		}
	}

	return v;
}

}
}
}
}

#endif
