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
#	include <ogdf/geometric/cr_min/geometry/objects/Direction.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Line.h>
#	include <ogdf/geometric/cr_min/geometry/objects/LineSegment.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Polygon.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Ray.h>

#	include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#	ifndef OGDF_GEOMETRIC_INEXACT_NUMBER_TYPE
#		define OGDF_GEOMETRIC_INEXACT_NUMBER_TYPE false
#	endif

namespace ogdf {
namespace internal {
namespace gcm {
namespace geometry {


struct Intersection {
	Intersection(unsigned int _seg_a, unsigned int _seg_b) : seg_a(_seg_a), seg_b(_seg_b) {
		//nothing to
	}

	unsigned int seg_a = -1; // first segment id
	unsigned int seg_b = -1; // second segment id
	unsigned int pos_on_a = -1; // ordered position on segment a
	unsigned int pos_on_b = -1; // ordered position on segment b
	bool proper = true; // is a proper intersection , i.e.,
			// intersection point is an interior point of both segments;
	bool is_source = false;
	bool is_target = false;

	bool is_incident(const unsigned int segment) const {
		return seg_a == segment || seg_b == segment;
	}

	unsigned int opposite(const unsigned int segment) const {
		if (seg_a == segment) {
			return seg_b;
		} else {
			return seg_a;
		}
	}

	unsigned int pos(const unsigned int segment) const {
		if (seg_a == segment) {
			return pos_on_a;
		} else {
			return pos_on_b;
		}
	}
};

/*! Compute a planar subdivision of a given set of segments
 */
template<typename Kernel, bool OPEN = true>
class PlanarSubdivision {
private:
	using LineSegment = geometry::LineSegment_t<Kernel>;
	using Point = geometry::Point_t<Kernel>;
	using Polygon = geometry::Polygon_t<Kernel>;
	using Line = geometry::Line_t<Kernel>;
	using Ray = geometry::Ray_t<Kernel>;
	using Direction = geometry::Direction_t<Kernel>;

	using ExactKernel = CGAL::Exact_predicates_exact_constructions_kernel;
	using ExactLineSegment = geometry::LineSegment_t<ExactKernel>;
	using ExactPoint = geometry::Point_t<ExactKernel>;


	enum class EventType { Start, End };

	struct Event {
		Point p;
		unsigned int segment_id;
		EventType event;

		Event(Point _p, unsigned int _id, EventType _ev) : p(_p), segment_id(_id), event(_ev) {
			//nothing to do
		}
	};


public:
	using LCGEdge = std::tuple<unsigned int, unsigned int, unsigned int>;

	/*! Compute sweep
	 * \param segments segments to intersect
	 * \param intersection_event callback
	 */
	template<typename IntersectionEvent>
	void sweep(std::vector<LineSegment>& segments, IntersectionEvent&& intersection_event) {
		std::vector<unsigned int> segment_to_active(segments.size(), -1);
		std::vector<unsigned int> active;

		auto add_to_active = [&](unsigned int segment_id) {
			segment_to_active[segment_id] = active.size();
			active.push_back(segment_id);
		};

		auto remove_from_active = [&](unsigned int segment_id) {
			unsigned int active_id = segment_to_active[segment_id];
			std::swap(active[active_id], active.back());
			segment_to_active[active[active_id]] = active_id;
			active.pop_back();
			segment_to_active[segment_id] = -1;
		};
		std::vector<Event> events;
		events.reserve(2 * segments.size());


		for (unsigned int i = 0; i < segments.size(); ++i) {
			Event ev_start(std::min(segments[i].source(), segments[i].target()), i, EventType::Start);
			Event ev_end(std::max(segments[i].source(), segments[i].target()), i, EventType::End);

			events.push_back(ev_start);
			events.push_back(ev_end);
		}


		std::sort(events.begin(), events.end(), [&](const Event& a, const Event& b) -> bool {
			return a.p < b.p || (a.p == b.p && a.event == EventType::Start);
		});


		for (const Event& e : events) {
			auto& new_segment = segments[e.segment_id];
			if (e.event == EventType::Start) {
				for (unsigned int active_id : active) {
					OGDF_ASSERT(active_id < segments.size());
					auto& active_segment = segments[active_id];

					if (OPEN && geometry::do_intersect_open(new_segment, active_segment)) {
						intersection_event(e.segment_id, active_id);
					}

					if (!OPEN && CGAL::do_intersect(new_segment, active_segment)) {
						intersection_event(e.segment_id, active_id);
					}
				}
				add_to_active(e.segment_id);
			} else {
				remove_from_active(e.segment_id);
			}
		}
	}

	//robustness depends on kernel
	void subdivision(std::vector<LineSegment>& segment,
			std::vector<Point>& intersections // node_to_point
			,
			std::vector<LCGEdge>& edge_list) {
		//compute intersections
		std::vector<std::vector<unsigned int>> intersections_of_segment(segment.size());

		//TODO handle start and end points;

		auto intersection_event = [&](unsigned int seg_1, unsigned int seg_2) {
			//TODO handle collinearity?

			unsigned int v = intersections.size();

			intersections_of_segment[seg_1].push_back(v);
			intersections_of_segment[seg_2].push_back(v);

			Point is = geometry::intersect(segment[seg_1], segment[seg_2]);
			intersections.push_back(is);

			OGDF_ASSERT(CGAL::squared_distance(segment[seg_1], is) < 1);
			OGDF_ASSERT(CGAL::squared_distance(segment[seg_2], is) < 1);

			OGDF_ASSERT(OGDF_GEOMETRIC_INEXACT_NUMBER_TYPE || segment[seg_1].has_on(is));
			OGDF_ASSERT(OGDF_GEOMETRIC_INEXACT_NUMBER_TYPE || segment[seg_2].has_on(is));
		};

		for (unsigned int i = 0; i < segment.size(); ++i) {
			auto& s = segment[i];
			intersections_of_segment[i].push_back(intersections.size());
			intersections.push_back(s.source());
			intersections_of_segment[i].push_back(intersections.size());
			intersections.push_back(s.target());
		}
		sweep(segment, intersection_event);

		//identify and remove duplication of intersections

		datastructure::UnionFind node_mapping(intersections.size());
		node_mapping.all_to_singletons();

		std::vector<unsigned int> same(intersections.size(), 0);
		for (unsigned int i = 0; i < same.size(); ++i) {
			same[i] = i;
		}

		std::sort(same.begin(), same.end(),
				[&](unsigned int a, unsigned int b) { return intersections[a] < intersections[b]; });

		for (unsigned int i = 0; i + 1 < same.size(); ++i) {
			unsigned int u = same[i];
			unsigned int v = same[i + 1];
			if (intersections[u] == intersections[v]) {
				node_mapping.merge(u, v);
			}
		}

		// build planar arrangement
		for (unsigned int i = 0; i < segment.size(); ++i) {
			auto& is = intersections_of_segment[i];
			//intersections += is.size();
			std::sort(is.begin(), is.end(), [&](const unsigned int a, unsigned int b) {
				const Point& p_a = intersections[a];
				const Point& p_b = intersections[b];

				auto d_a = CGAL::squared_distance(segment[i].source(), p_a);
				auto d_b = CGAL::squared_distance(segment[i].source(),
						p_b); //order a long line with respect to source
				return d_a < d_b;
			});
			for (unsigned int j = 0; j + 1 < is.size(); ++j) {
				unsigned int u = node_mapping[is[j]];
				unsigned int v = node_mapping[is[j + 1]];

				if (u != v) {
					OGDF_ASSERT(OGDF_GEOMETRIC_INEXACT_NUMBER_TYPE
							|| segment[i].has_on(intersections[u]));
					OGDF_ASSERT(OGDF_GEOMETRIC_INEXACT_NUMBER_TYPE
							|| segment[i].has_on(intersections[v]));

					edge_list.push_back(LCGEdge(u, v, i));
				}
			}
		}
	}

	std::vector<Intersection> subdivision(std::vector<LineSegment>& segment) {
		//compute intersections
		std::vector<Intersection> intersecting_segments;
		//TODO handle start and end points;

		auto intersection_event = [&](unsigned int seg_1, unsigned int seg_2) {
			intersecting_segments.push_back({seg_1, seg_2});
		};

		sweep(segment, intersection_event);
		return intersecting_segments;
	}
};


}
}
}
}

#endif
