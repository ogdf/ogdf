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

#	include <ogdf/geometric/cr_min/graph/BloatedDual.h>

#	ifdef OGDF_GEOMETRIC_CR_MIN_DEBUG
#		include <ogdf/basic/Logger.h>
#		define OGDF_GEOMETRIC_BD_ECHO(x) Logger::slout() << "[BD] " << x << std::endl;
#	else
#		define OGDF_GEOMETRIC_BD_ECHO(x)
#	endif

namespace ogdf {
namespace internal {
namespace gcm {
namespace geometry {

template<typename Kernel, bool parallel = false>
class BloatedDual {
private:
	using Point = geometry::Point_t<Kernel>;
	using Segment = LineSegment_t<Kernel>;
	using ExactKernel = CGAL::Simple_cartesian<CGAL::Gmpq>;

	using BD = graph::BloatedDualGraph<Kernel>;
	using Node = typename BD::Node;
	using Edge = typename BD::Edge;
	using LocalIntersection = typename BD::LocalIntersection;

	inline bool consider_equal(const Point& p, const Point& q) const { return p == q; }

	bool has_point_on_segment(const Segment& reference, const Point& p) {
		if (reference.has_on(p)) {
			return true;
		} else if (consider_equal(p, reference.supporting_line().projection(p))) {
			auto ex_p = geometry::cast<ExactKernel>(p);
			auto ex_segment = geometry::cast<ExactKernel>(reference);
			return ex_segment.has_on(ex_p);
		}

		return false;
	}

	inline bool has_endpoint_on_segment(const Segment& reference, const Segment& opposite) {
		return has_point_on_segment(reference, opposite.source())
				|| has_point_on_segment(reference, opposite.target());
	}

	std::set<std::pair<unsigned int, unsigned int>> seen;

	void assign_intersections_to_segment(const std::vector<Segment>& segments,
			const std::vector<Intersection>& intersecting_segments,
			std::vector<std::vector<unsigned int>>& _segment_to_intersections,
			std::vector<Point>& intersections) {
		//TODO set really required?
		seen.clear();
		for (unsigned int i = 0; i < intersecting_segments.size(); ++i) {
			auto& pair = intersecting_segments[i];
			if (seen.find({pair.seg_a, pair.seg_b}) == seen.end()
					&& seen.find({pair.seg_b, pair.seg_a}) == seen.end()) {
				_segment_to_intersections[pair.seg_a].push_back(i);
				if (pair.seg_a != pair.seg_b) {
					_segment_to_intersections[pair.seg_b].push_back(i);
				}

				auto& seg_a = segments[pair.seg_a];
				auto& seg_b = segments[pair.seg_b];
				if (pair.is_source) {
					OGDF_ASSERT(pair.seg_a == pair.seg_b);
					intersections.push_back(seg_a.source());
				} else if (pair.is_target) {
					OGDF_ASSERT(pair.seg_a == pair.seg_b);
					intersections.push_back(seg_a.target());
				} else if (geometry::have_common_endpoints(seg_a, seg_b)) {
					intersections.push_back(geometry::get_common_endpoint(seg_a, seg_b));
				} else {
					Point p_a = geometry::intersect(seg_a, seg_b);
					intersections.push_back(p_a);
				}
			}
		}
	}

	void sort_intersections_along_segment(const std::vector<Segment>& segments,
			const std::vector<Intersection>& intersecting_segments,
			std::vector<std::vector<unsigned int>>& _segment_to_intersections,
			std::vector<Point>& intersections) {
		for (unsigned int i = 0; i < _segment_to_intersections.size(); ++i) {
			auto& is = _segment_to_intersections[i];

			auto compare = [&](const unsigned int a, const unsigned int b) {
				const Point& p_a = intersections[a];

				const Point& p_b = intersections[b];

				auto d_a = CGAL::squared_distance(segments[i].source(), p_a);
				auto d_b = CGAL::squared_distance(segments[i].source(),
						p_b); //order a long line with respect to source

				bool comp = d_a < d_b;
				return comp;
			};

			if (parallel) {
				//boost::sort::block_indirect_sort(is.begin(), is.end(), compare);
				std::sort(is.begin(), is.end(), compare);
			} else {
				std::sort(is.begin(), is.end(), compare);
			}
		}
	}

	std::vector<unsigned int> left_intersections;
	std::vector<unsigned int> right_intersections;

	/* If several segments intersect in on point on @segment only two segments
	 * are necessary to construct the bloated dual. This function filters
	 * all non-relevant segments.
	 */
	void clean_up(unsigned int segment, const std::vector<Segment>& segments,
			const std::vector<Intersection>& intersecting_segments //maps intersection id to pair of segments
			,
			const std::vector<Point>& intersections // maps intersection id to intersections
			,
			const std::vector<unsigned int>& intersections_along_segment //maps to intersection_id
			,
			std::vector<LocalIntersection>& clean_intersecting_segments // filtered intersection_id's
	) {
		auto is_proper = [&](const Segment& _segment, const Point& p) {
			return !consider_equal(_segment.source(), p) && !consider_equal(_segment.target(), p)
					&& !has_point_on_segment(_segment, p);
		};

		auto is_proper_left_turn = [&](const Segment& _segment, const Point& p) {
			return is_proper(_segment, p) && geometry::left_turn(_segment, p);
		};


		auto is_proper_right_turn = [&](const Segment& _segment, const Point& p) {
			return is_proper(_segment, p) && geometry::right_turn(_segment, p);
		};


		for (unsigned int j = 0; j
				< intersections_along_segment.size();) { //incrementation is done in the next while loop
			left_intersections.clear();
			right_intersections.clear();

			const unsigned int ref_intersection = j;

			int self_intersection_id = -1;

			OGDF_GEOMETRIC_BD_ECHO("reference: " << segment);

			while (j < intersections_along_segment.size()
					&& consider_equal(intersections[intersections_along_segment[ref_intersection]],
							intersections[intersections_along_segment[j]])) { //TODO make robust?

				unsigned int intersection_id = intersections_along_segment[j];
				unsigned int opp = intersecting_segments[intersection_id].opposite(segment);

				bool is_left_turn = is_proper_left_turn(segments[segment], segments[opp].source())
						|| is_proper_left_turn(segments[segment], segments[opp].target());
				bool is_right_turn = is_proper_right_turn(segments[segment], segments[opp].source())
						|| is_proper_right_turn(segments[segment], segments[opp].target());


				if (segment == opp) {
					self_intersection_id = intersection_id;
				} else if (!is_left_turn && !is_right_turn) { //is collinear

					left_intersections.push_back(intersection_id);
					right_intersections.push_back(intersection_id);
				} else {
					if (is_left_turn) {
						left_intersections.push_back(intersection_id);
					}

					if (is_right_turn) {
						right_intersections.push_back(intersection_id);
					}
				}

				++j;
			}

			if (left_intersections.empty() && right_intersections.empty()) {
				left_intersections.push_back(self_intersection_id);
			}

			if (left_intersections.size() > 1) {
				std::sort(left_intersections.begin(), left_intersections.end(),
						[&](const unsigned int a, const unsigned int b) {
							auto& seg_a = segments[intersecting_segments[a].opposite(segment)];
							auto& seg_b = segments[intersecting_segments[b].opposite(segment)];

							int sign_a = 1;
							int sign_b = 1;
							if (is_proper_left_turn(segments[segment], seg_a.source())
									|| segments[segment].source() == seg_a.target()) {
								sign_a = -1;
							}

							if (is_proper_left_turn(segments[segment], seg_b.source())
									|| segments[segment].source() == seg_b.target()) {
								sign_b = -1;
							}

							return geometry::right_turn(seg_a.to_vector() * sign_a,
									seg_b.to_vector() * sign_b);
						});
			}


			if (right_intersections.size() > 1) {
				std::sort(right_intersections.begin(), right_intersections.end(),
						[&](const unsigned int a, const unsigned int b) {
							auto& seg_a = segments[intersecting_segments[a].opposite(segment)];
							auto& seg_b = segments[intersecting_segments[b].opposite(segment)];

							int sign_a = 1;
							int sign_b = 1;

							if (is_proper_right_turn(segments[segment], seg_a.source())
									|| segments[segment].source() == seg_a.target()) {
								sign_a = -1;
							}

							if (is_proper_right_turn(segments[segment], seg_b.source())
									|| segments[segment].source() == seg_b.target()) {
								sign_b = -1;
							}

							return geometry::left_turn(seg_a.to_vector() * sign_a,
									seg_b.to_vector() * sign_b);
						});
			}

			OGDF_ASSERT(!left_intersections.empty() || !right_intersections.empty());
			LocalIntersection li;
			li.is_end_point = self_intersection_id >= 0;
			if (!left_intersections.empty()) {
				li.left_intersection_ids[0] = left_intersections.front();
				li.left_intersection_ids[1] = left_intersections.back();
			}

			if (!right_intersections.empty()) {
				li.right_intersection_ids[0] = right_intersections.front();
				li.right_intersection_ids[1] = right_intersections.back();
			}

			clean_intersecting_segments.push_back(li);
		}
	}

	void handle_non_proper_intersection(BD& bd, unsigned int ref_seg, LocalIntersection& li) {
		if (li.left_is_valid() //self intersection is stored on the left side
				&& (bd.intersecting_segments[li.left_intersection_ids[0]].is_source
						|| bd.intersecting_segments[li.left_intersection_ids[0]].is_target)) {
			return;
		}

		OGDF_GEOMETRIC_BD_ECHO("---------------------------");
		OGDF_GEOMETRIC_BD_ECHO("reference segment " << ref_seg);

		// assign segments to left / right, depending on
		// whether they have an endpoint to left / right of the reference segments;

		bool intersection_is_source =
				consider_equal(bd.segments[ref_seg].source(), bd.get_intersection_point(li));
		bool intersection_is_target =
				consider_equal(bd.segments[ref_seg].target(), bd.get_intersection_point(li));
		bool proper = !intersection_is_source && !intersection_is_target;

		bool source_is_left[2] = {false, false};
		bool source_is_right[2] = {false, false};

		if (li.left_is_valid()) {
			unsigned int left_segment_id[2];
			OGDF_GEOMETRIC_BD_ECHO("ref: " << ref_seg << "; " << bd.segments[ref_seg]);
			OGDF_GEOMETRIC_BD_ECHO("intersection is source: " << intersection_is_source)
			OGDF_GEOMETRIC_BD_ECHO("intersection is target: " << intersection_is_target)
			OGDF_GEOMETRIC_BD_ECHO("proper: " << proper)
			for (unsigned int i = 0; i < 2; ++i) {
				left_segment_id[i] =
						bd.intersecting_segments[li.left_intersection_ids[i]].opposite(ref_seg);

				source_is_left[i] = !consider_equal(bd.segments[left_segment_id[i]].source(),
											bd.get_intersection_point(li))
						&& geometry::left_turn(bd.segments[ref_seg],
								bd.segments[left_segment_id[i]].source());

				OGDF_GEOMETRIC_BD_ECHO("\t l" << left_segment_id[i] << "; "
											  << bd.segments[left_segment_id[i]] << ": "
											  << source_is_left[i]);
			}

			if (proper || intersection_is_target) {
				Node u1 = bd.first_left(ref_seg, li.left_intersection_ids[0]);
				Node v1 = -1;

				if (source_is_left[0]
						|| consider_equal(bd.segments[left_segment_id[0]].target(),
								bd.get_intersection_point(li))) {
					v1 = bd.first_right(left_segment_id[0], li.left_intersection_ids[0]);
				} else {
					v1 = bd.second_left(left_segment_id[0], li.left_intersection_ids[0]);
				}

				bd.add_edge(u1, v1, true);
			}

			if (proper || intersection_is_source) {
				Node u2 = bd.second_left(ref_seg, li.left_intersection_ids[0]);
				Node v2 = -1;

				if (source_is_left[1]
						|| consider_equal(bd.segments[left_segment_id[1]].target(),
								bd.get_intersection_point(li))) {
					v2 = bd.first_left(left_segment_id[1], li.left_intersection_ids[1]);
				} else {
					v2 = bd.second_right(left_segment_id[1], li.left_intersection_ids[1]);
				}

				bd.add_edge(u2, v2, true);
			}

		} else if (proper) {
			Node u1 = bd.first_left(ref_seg, li.right_intersection_ids[0]);
			Node v1 = bd.second_left(ref_seg, li.right_intersection_ids[0]);
			bd.add_edge(u1, v1, true);
		} else if (intersection_is_source) {
			auto id = li.right_intersection_ids[0];

			Node u1 = bd.second_left(ref_seg, id);
			Node v1 = -1;
			unsigned int opp = bd.intersecting_segments[id].opposite(ref_seg);
			if (consider_equal(bd.segments[opp].target(), bd.get_intersection_point(li))) {
				v1 = bd.first_left(opp, id);
			} else {
				v1 = bd.second_right(opp, id);
			}

			bd.add_edge(u1, v1, true);

		} else if (intersection_is_target) {
			auto id = li.right_intersection_ids[1];

			Node u1 = bd.first_left(ref_seg, id);
			Node v1 = -1;
			unsigned int opp = bd.intersecting_segments[id].opposite(ref_seg);
			if (consider_equal(bd.segments[opp].target(), bd.get_intersection_point(li))) {
				v1 = bd.first_right(opp, id);
			} else {
				v1 = bd.second_left(opp, id);
			}

			bd.add_edge(u1, v1, true);
		}

		if (li.right_is_valid()) {
			unsigned int right_segment_id[2];
			OGDF_GEOMETRIC_BD_ECHO("ref: " << ref_seg << "; " << bd.segments[ref_seg]);
			for (unsigned int i = 0; i < 2; ++i) {
				right_segment_id[i] =
						bd.intersecting_segments[li.right_intersection_ids[i]].opposite(ref_seg);

				source_is_right[i] = !consider_equal(bd.segments[right_segment_id[i]].source(),
											 bd.get_intersection_point(li))
						&& geometry::right_turn(bd.segments[ref_seg],
								bd.segments[right_segment_id[i]].source());

				OGDF_GEOMETRIC_BD_ECHO("\t r" << right_segment_id[i] << "; "
											  << bd.segments[right_segment_id[i]] << ": "
											  << source_is_right[i]);
			}

			if (proper || intersection_is_target) {
				Node u1 = bd.first_right(ref_seg, li.right_intersection_ids[0]);
				Node v1 = -1;
				if (source_is_right[0]
						|| consider_equal(bd.segments[right_segment_id[0]].target(),
								bd.get_intersection_point(li))) {
					v1 = bd.first_left(right_segment_id[0], li.right_intersection_ids[0]);
				} else {
					v1 = bd.second_right(right_segment_id[0], li.right_intersection_ids[0]);
				}

				bd.add_edge(u1, v1, true);
			}

			if (proper || intersection_is_source) {
				Node u2 = bd.second_right(ref_seg, li.right_intersection_ids[0]);
				Node v2 = -1;

				if (source_is_right[1]
						|| consider_equal(bd.segments[right_segment_id[1]].target(),
								bd.get_intersection_point(li))) {
					v2 = bd.first_right(right_segment_id[1], li.right_intersection_ids[1]);
				} else {
					v2 = bd.second_left(right_segment_id[1], li.right_intersection_ids[1]);
				}

				bd.add_edge(u2, v2, true);
			}
		} else if (proper) {
			Node u1 = bd.first_right(ref_seg, li.left_intersection_ids[0]);
			Node v1 = bd.second_right(ref_seg, li.left_intersection_ids[0]);
			bd.add_edge(u1, v1, true);

		} else if (intersection_is_source) {
			auto id = li.left_intersection_ids[0];

			Node u1 = bd.second_right(ref_seg, id);
			Node v1 = -1;
			unsigned int opp = bd.intersecting_segments[id].opposite(ref_seg);
			if (consider_equal(bd.segments[opp].target(), bd.get_intersection_point(li))) {
				v1 = bd.first_right(opp, id);
			} else {
				v1 = bd.second_left(opp, id);
			}

			bd.add_edge(u1, v1, true);
		} else if (intersection_is_target) {
			auto id = li.left_intersection_ids[1];

			Node u1 = bd.first_right(ref_seg, id);
			Node v1 = -1;
			unsigned int opp = bd.intersecting_segments[id].opposite(ref_seg);
			if (consider_equal(bd.segments[opp].target(), bd.get_intersection_point(li))) {
				v1 = bd.first_left(opp, id);
			} else {
				v1 = bd.second_right(opp, id);
			}

			bd.add_edge(u1, v1, true);
		}
	}

	inline void construct_graph(BD& bd) {
		int nof_threads = 1;
		(void)nof_threads;
		if (parallel) {
			nof_threads = omp_get_max_threads();
		}
#	pragma omp parallel for num_threads(nof_threads)
		for (unsigned int i = 0; i < bd.segments.size(); ++i) {
			for (unsigned int j = 0; j < bd.segment_to_intersections[i].size(); ++j) {
				handle_non_proper_intersection(bd, i, bd.segment_to_intersections[i][j]);
			}
		}
	}

public:
	// a bloated dual has a number of vertices in each face corresponding to the number of edges on the face
	// vertices within a face are connected via circle corresponding to order of the edges of the face
	// functions assumes that endpoints of semgents are in lexicographical order
	// @ASSUMPTION: no overlapping segments
	// @ASSUMPTION: no segment contains end endpoint of another segment in its interior....
	// @return a vector containg a mapping an edge to pair of segment ids. if the edge crosses a segment,
	// then the entries coincide.


	std::vector<std::vector<unsigned int>> segment_to_intersections;

	void construct(BD& bd) {
		//sort intersections along segment
		segment_to_intersections.clear();
		segment_to_intersections.resize(bd.segments.size());

		for (unsigned int i = 0; i < bd.segments.size(); ++i) {
			Intersection s(i, i);

			s.is_source = true;
			bd.intersecting_segments.push_back(s);
			s.is_source = false;

			s.is_target = true;
			bd.intersecting_segments.push_back(s);
		}

		assign_intersections_to_segment(bd.segments, bd.intersecting_segments,
				segment_to_intersections, bd.intersections);

		sort_intersections_along_segment(bd.segments, bd.intersecting_segments,
				segment_to_intersections, bd.intersections);

		bd.segment_to_intersections.resize(bd.segments.size());
		for (unsigned int i = 0; i < bd.segments.size(); ++i) {
			clean_up(i, bd.segments, bd.intersecting_segments, bd.intersections,
					segment_to_intersections[i], bd.segment_to_intersections[i]);
			OGDF_ASSERT(segment_to_intersections[i].size() <= 1
					|| bd.segment_to_intersections[i].size() >= 1);
		}

		for (unsigned int i = 0; i < bd.segment_to_intersections.size(); ++i) { // i: i-th segment
			auto& is = bd.segment_to_intersections[i];
			//assign positions
			for (unsigned int j = 0; j < is.size(); ++j) { // j: j-th intersection on segment i
				LocalIntersection& x = is[j];

				auto set_pos = [&](unsigned int intersection_id) {
					if (intersection_id < bd.intersecting_segments.size()) {
						auto& y = bd.intersecting_segments[intersection_id];
						if (y.seg_a == i) {
							y.pos_on_a = j;
						}
						if (y.seg_b == i) {
							y.pos_on_b = j;
						}
					}
				};

				set_pos(x.left_intersection_ids[0]);
				set_pos(x.left_intersection_ids[1]);
				set_pos(x.right_intersection_ids[0]);
				set_pos(x.right_intersection_ids[1]);
			}
		}

		bd.segm_to_node_range.reserve(bd.intersections.size());
		bd.edges.reserve(3 * bd.intersections.size());
		bd.node_to_segment.reserve(bd.intersections.size());
		;

		//add all vertices
		for (unsigned int i = 0; i < bd.segments.size(); ++i) {
			bd.segm_to_node_range.push_back(bd.number_of_nodes());
			//add #is + 1 nodes
			auto& is = bd.segment_to_intersections[i];

			for (unsigned int j = 0; j + 1 < is.size(); ++j) {
				Node left = bd.add_node(i);
				Node right = bd.add_node(i);
				bd.add_edge(left, right, false);
			}
		}
		bd.segm_to_node_range.push_back(bd.number_of_nodes());

		construct_graph(bd);
	}

	static std::vector<Segment> compute_drawing(graph::BloatedDualGraph<Kernel>& bd) {
		std::vector<Point> node_to_point(bd.number_of_nodes());
		std::vector<Segment> edge_segments;

		for (unsigned int i = 0; i < bd.segments.size(); ++i) {
			auto left_dir = geometry::normalize(bd.segments[i].to_vector())
									.perpendicular(CGAL::COUNTERCLOCKWISE);
			auto right_dir =
					geometry::normalize(bd.segments[i].to_vector()).perpendicular(CGAL::CLOCKWISE);
			if (i >= bd.segment_to_intersections.size()) {
				continue;
			}

			auto& is = bd.segment_to_intersections[i];
			OGDF_GEOMETRIC_BD_ECHO(bd.segments[i]);

			// source and target should be intersections
			OGDF_ASSERT(is.size() >= 2);
			Point last = bd.intersections[is[0].intersection_id()];

			for (unsigned int j = 0; j + 1 < is.size(); j = j + 1) {
				Point current = bd.intersections[is[j + 1].intersection_id()];

				auto mid = last + (current - last) * 0.5;

				auto left = bd.segm_to_node_range[i] + 2 * j;
				auto right = bd.segm_to_node_range[i] + 2 * j + 1;

				double c = std::max(
						CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(current, mid))) / 10, 0.0);
				node_to_point[left] = mid + left_dir * c;
				node_to_point[right] = mid + right_dir * c;
				last = current;
			}
		}

		auto add_segment = [&](Node u, Node v) {
			if (u != v) {
				edge_segments.push_back({node_to_point[u], node_to_point[v]});
			}
		};

		for (unsigned int v = 0; v < bd.number_of_nodes(); ++v) {
			add_segment(v, bd.edges[3 * v]);
			add_segment(v, bd.edges[3 * v + 1]);
			add_segment(v, bd.edges[3 * v + 2]);
		}

		return edge_segments;
	}
};

}
}
}
}

#endif
