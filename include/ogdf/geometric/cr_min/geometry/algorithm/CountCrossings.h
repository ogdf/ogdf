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

#	include <ogdf/geometric/cr_min/datastructure/OGDFVector.h>
#	include <ogdf/geometric/cr_min/geometry/objects/LineSegment.h>
#	include <ogdf/geometric/cr_min/graph/GeometricDrawing.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace geometry {
template<typename Kernel>
int count_crossings(std::vector<LineSegment_t<Kernel>>& segments) {
	unsigned int crossings = 0;
	for (unsigned int i = 0; i + 1 < segments.size(); ++i) {
		for (unsigned int j = i + 1; j < segments.size(); ++j) {
			crossings += geometry::do_intersect_open(segments[i], segments[j]);
		}
	}

	return crossings;
}

template<typename Kernel, typename Graph>
int count_crossings(graph::GeometricDrawing<Kernel, Graph>& drawing) {
	std::vector<LineSegment_t<Kernel>> segments;

	for (auto e : drawing.get_graph().edges()) {
		segments.push_back(drawing.get_segment(e));
	}
	return count_crossings(segments);
}

template<typename Kernel, typename Graph>
int count_crossings(graph::GeometricDrawing<Kernel, Graph>& d, typename Graph::Node& v) {
	int cr = 0;
	for (auto f : d.get_graph().edges()) {
		for (auto e : d.get_graph().edges(v)) {
			cr += (f != e) && geometry::do_intersect_open(d.get_segment(e), d.get_segment(f));
		}
	}
	return cr;
}

template<typename Kernel, typename Graph>
std::vector<int> count_crossings_vec(const graph::GeometricDrawing<Kernel, Graph>& d,
		const typename Graph::Node& v) {
	std::vector<int> cr(v->degree(), 0);

	for (auto f : d.get_graph().edges()) {
		unsigned int i = 0;
		for (auto e : d.get_graph().edges(v)) {
			cr[i] += (f != e) && geometry::do_intersect_open(d.get_segment(e), d.get_segment(f));
			++i;
		}
	}
	return cr;
}

template<typename Kernel, typename Graph>
int count_crossing_edges(graph::GeometricDrawing<Kernel, Graph>& d, typename Graph::Node& v) {
	int cr = 0;
	auto& g = d.get_graph();

	for (auto f : g.edges()) {
		bool counted = false;
		for (auto e : g.edges(v)) {
			if ((f != e) && !counted
					&& geometry::do_intersect_open(d.get_segment(e), d.get_segment(f))) {
				++cr;
				counted = true;
			}
		}
	}
	return cr;
}

template<typename Kernel, typename Graph>
int count_crossings(graph::GeometricDrawing<Kernel, Graph>& d, typename Graph::Edge& e) {
	int cr = 0;
	for (auto f : d.get_graph().edges()) {
		cr += (f != e) && (f->commonNode(e) == nullptr)
				&& geometry::do_intersect_open(d.get_segment(e), d.get_segment(f));
	}
	return cr;
}

}
}
}
}

#endif
