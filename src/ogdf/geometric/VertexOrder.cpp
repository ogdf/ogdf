/** \file
 * \brief Computes an vertex order based on the number of crossings in a given (straight-line) drawing.
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

#include <ogdf/geometric/VertexOrder.h>

#include <queue>

#ifdef OGDF_INCLUDE_CGAL

#	include <ogdf/geometric/cr_min/geometry/algorithm/CountCrossings.h>

#	include <CGAL/Gmpq.h>
#	include <CGAL/Simple_cartesian.h>

namespace ogdf {


void CrossingVertexOrder::sort() {
	auto asc = [&](const QElement& a, const QElement& b) { return a.second < b.second; };

	auto desc = [&](const QElement& a, const QElement& b) { return a.second > b.second; };

	if (o == OrderEnum::asc) {
		std::sort(vertex_order.begin(), vertex_order.end(), asc);
	} else {
		std::sort(vertex_order.begin(), vertex_order.end(), desc);
	}
}

double CrossingVertexOrder::crossings(int c) {
	if (m == MeasureEnum::log) {
		return std::log(c + 1);
	} else if (m == MeasureEnum::squared) {
		return c * c;
	} else {
		return c;
	}
}

void CrossingVertexOrder::init() {
	using _Graph = internal::gcm::graph::OGDFGraphWrapper;

	using _Kernel = CGAL::Simple_cartesian<double>; //todo
	using Drawing = internal::gcm::graph::GeometricDrawing<_Kernel, _Graph>;

	_Graph g(ga.constGraph());
	Drawing d(g);

	internal::gcm::graph::ogdf_attributes_to_geometric_drawing(ga, d);

	vertex_order.clear();
	if (o == OrderEnum::rnd) {
		for (auto v : d.get_graph().nodes()) {
			vertex_order.push_back({v, 1});
		}
		std::mt19937_64 rd(ogdf::randomSeed());
		std::shuffle(vertex_order.begin(), vertex_order.end(), rd);
	} else {
		for (auto v : d.get_graph().nodes()) {
			std::vector<int> cr = internal::gcm::geometry::count_crossings_vec(d, v);
			int x = 0;
			for (auto c : cr) {
				x += crossings(c);
			}
			vertex_order.push_back(std::make_pair(v, x));
		}
		sort();
	}
}

void CrossingVertexOrder::init_cr(ogdf::edge e) {
	using _Graph = internal::gcm::graph::OGDFGraphWrapper;

	using _Kernel = CGAL::Simple_cartesian<double>; //todo
	using Drawing = internal::gcm::graph::GeometricDrawing<_Kernel, _Graph>;

	_Graph g(ga.constGraph());
	Drawing d(g);

	internal::gcm::graph::ogdf_attributes_to_geometric_drawing(ga, d);

	vertex_order.clear();

	if (m != MeasureEnum::zero) {
		for (auto w : d.get_graph().nodes()) {
			unsigned int cr = 0;
			for (auto f : d.get_graph().edges(w)) {
				cr += CGAL::do_intersect(d.get_segment(e), d.get_segment(f));
			}
			if (cr > 0) {
				vertex_order.push_back({w, crossings(cr)});
			}
		}
		sort();
	}
}


}
#else

namespace ogdf {

void CrossingVertexOrder::sort() { }

double CrossingVertexOrder::crossings(int) { return 0; }

void CrossingVertexOrder::init() {
	OGDF_THROW_PARAM(LibraryNotSupportedException, LibraryNotSupportedCode::Cgal);
}

void CrossingVertexOrder::init_cr(ogdf::edge e) {
	OGDF_THROW_PARAM(LibraryNotSupportedException, LibraryNotSupportedCode::Cgal);
}


}

#endif
