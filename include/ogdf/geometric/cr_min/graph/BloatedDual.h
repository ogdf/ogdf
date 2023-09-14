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
#	include <ogdf/geometric/cr_min/geometry/algorithm/PlanarSubdivision.h>
#	include <ogdf/geometric/cr_min/geometry/objects/LineSegment.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Point.h>

#	include <vector>

#	include <omp.h>

namespace ogdf {
namespace internal {
namespace gcm {

namespace graph {
template<typename Kernel, bool parallel = false>
class BloatedDualGraph {
private:
	using Segment = geometry::LineSegment_t<Kernel>;
	using Point = geometry::Point_t<Kernel>;

	unsigned int m_number_of_nodes = 0;

	std::vector<omp_lock_t> locks;

public:
	//intersection from the perspective of a @reference_segment
	struct LocalIntersection {
		unsigned int left_intersection_ids[2] = {(unsigned int)-1, (unsigned int)-1};
		unsigned int right_intersection_ids[2] = {(unsigned int)-1, (unsigned int)-1};
		bool is_end_point = false;

		inline bool left_is_valid() const { return left_intersection_ids[0] != (unsigned int)-1; }

		inline bool right_is_valid() const { return right_intersection_ids[0] != (unsigned int)-1; }

		unsigned int intersection_id() const {
			if (left_is_valid()) {
				return left_intersection_ids[0];
			} else {
				return right_intersection_ids[0];
			}
		}
	};

	using Node = unsigned int;
	using Edge = unsigned int;

	std::vector<Segment> segments;
	std::vector<geometry::Intersection> intersecting_segments;
	std::vector<Point> intersections;

	std::vector<unsigned int> segm_to_node_range; //TODO actually edge range, now
	std::vector<Node> edges;
	std::vector<unsigned int> node_to_segment;
	std::vector<std::vector<LocalIntersection>> segment_to_intersections;

	void clear() {
		m_number_of_nodes = 0;
		segments.clear();
		intersecting_segments.clear();
		intersections.clear();
		segm_to_node_range.clear();
		edges.clear();
		node_to_segment.clear();
		segment_to_intersections.clear();
		locks.clear();
	}

	BloatedDualGraph() { }

	unsigned int number_of_nodes() const { return m_number_of_nodes; }

	unsigned int degree(const Node v) const {
		if (edges[3 * v + 2] != v) {
			return 3;
		} else if (edges[3 * v + 1] != v) {
			return 2;
		} else if (edges[3 * v] != v) {
			return 1;
		} else {
			return 0;
		}
	}

	Node add_node(unsigned int segment_id) {
		edges.push_back(m_number_of_nodes); // initial every edge is a self loop
		edges.push_back(m_number_of_nodes);
		edges.push_back(m_number_of_nodes);

		node_to_segment.push_back(segment_id);

		if (parallel) {
			locks.push_back(omp_lock_t());
		}

		return m_number_of_nodes++;
	}

	void add_edge(const Node u, const Node v, bool face_edge) {
		OGDF_ASSERT(u < number_of_nodes());
		OGDF_ASSERT(v < number_of_nodes());

		auto add_edge = [&](Node _u, Node _v) {
			if (parallel) {
				omp_set_lock(&locks[_u]);
			}

			unsigned int i = _u * 3;
			OGDF_ASSERT(i < edges.size());
			if (edges[i] != _v && edges[i + 1] != _v && edges[i + 2] != _v) {
				OGDF_ASSERT(edges[i + 2] == _u);

				unsigned int off = 0;
				if (face_edge) {
					off = 1;
					if (edges[i + 1] != _u) {
						off = 2; //self loop serves as sentinel
					}
				}

				edges[i + off] = _v;
			}
			if (parallel) {
				omp_unset_lock(&locks[_u]);
			}
		};

		add_edge(u, v);
		add_edge(v, u);
	}

	bool is_left(const Node v) const {
		OGDF_ASSERT(v < number_of_nodes());
		return v % 2 == 0;
	}

	bool is_right(const Node v) const {
		OGDF_ASSERT(v < number_of_nodes());
		return v % 2 == 1;
	}

	Point get_intersection_point(const LocalIntersection& li) const {
		unsigned int intersection_id = -1;
		if (li.left_is_valid()) {
			intersection_id = li.left_intersection_ids[0];
		} else {
			intersection_id = li.right_intersection_ids[1];
		}
		return intersections[intersection_id];
	}

	unsigned int get_first_id(const unsigned int segment, const unsigned int intersection_id) const {
		OGDF_ASSERT(intersection_id < intersecting_segments.size());

		const geometry::Intersection& is = intersecting_segments[intersection_id];
		OGDF_ASSERT(is.is_incident(segment));
		unsigned int pos = 2 * (is.pos(segment) - 1); //2 accounts for the left and right vertices...
		OGDF_ASSERT(segment < segment_to_intersections.size());
		OGDF_ASSERT(segm_to_node_range[segment] + pos < segm_to_node_range[segment + 1]);
		return segm_to_node_range[segment] + pos;
	}

	unsigned int get_second_id(const unsigned int segment, const unsigned int intersection_id) const {
		OGDF_ASSERT(intersection_id < intersecting_segments.size());

		const geometry::Intersection& is = intersecting_segments[intersection_id];
		OGDF_ASSERT(is.is_incident(segment));

		unsigned int pos = 2 * is.pos(segment);

		OGDF_ASSERT(segment < segment_to_intersections.size());
		OGDF_ASSERT(segm_to_node_range[segment] + pos < segm_to_node_range[segment + 1]);

		return segm_to_node_range[segment] + pos;
	}

	unsigned int get_first_id(const unsigned int segment, const LocalIntersection& li) const {
		if (li.left_is_valid()) {
			return get_first_id(segment, li.left_intersection_ids[0]);
		} else {
			return get_first_id(segment, li.right_intersection_ids[0]);
		}
	}

	unsigned int get_second_id(const unsigned int segment, const LocalIntersection& li) const {
		if (li.left_is_valid()) {
			return get_second_id(segment, li.left_intersection_ids[0]);
		} else {
			return get_second_id(segment, li.right_intersection_ids[0]);
		}
	}

	Node first_left(const unsigned int segment, const unsigned int intersection_id) const {
		OGDF_ASSERT(get_first_id(segment, intersection_id) < number_of_nodes());
		return get_first_id(segment, intersection_id);
	}

	Node first_right(const unsigned int segment, const unsigned int intersection_id) const {
		OGDF_ASSERT(get_first_id(segment, intersection_id) + 1 < number_of_nodes());
		return get_first_id(segment, intersection_id) + 1;
	}

	Node second_left(const unsigned int segment, const unsigned int intersection_id) const {
		OGDF_ASSERT(get_second_id(segment, intersection_id) < number_of_nodes());
		return get_second_id(segment, intersection_id);
	}

	Node second_right(const unsigned int segment, const unsigned int intersection_id) const {
		OGDF_ASSERT(get_second_id(segment, intersection_id) + 1 < number_of_nodes());
		return get_second_id(segment, intersection_id) + 1;
	}

	Node first_left(const unsigned int segment, const LocalIntersection& li) const {
		OGDF_ASSERT(get_first_id(segment, li) < number_of_nodes());
		return get_first_id(segment, li);
	}

	Node first_right(const unsigned int segment, const LocalIntersection& li) const {
		OGDF_ASSERT(get_first_id(segment, li) + 1 < number_of_nodes());
		return get_first_id(segment, li) + 1;
	}

	Node second_left(const unsigned int segment, const LocalIntersection& li) const {
		OGDF_ASSERT(get_second_id(segment, li) < number_of_nodes());
		return get_second_id(segment, li);
	}

	Node second_right(const unsigned int segment, const LocalIntersection& li) const {
		OGDF_ASSERT(get_second_id(segment, li) + 1 < number_of_nodes());
		return get_second_id(segment, li) + 1;
	}

	bool is_face_edge(const Edge e) const {
		OGDF_ASSERT(e < 3 * number_of_nodes());
		return e % 3 != 0;
	}
};

}
}
}
}

#endif
