/** \file
 * \brief In the Geometric Vertex Insertion Approach we iteratively insert a given set of vertices into a drawing of a graph.  Each vertex is placed
 * at its optimal position. The optimal position can be for example a position that minimizes the number of crossings for this vertex.
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

#include <ogdf/geometric/GeometricVertexInsertion.h>

#ifdef OGDF_INCLUDE_CGAL

#	include <ogdf/geometric/cr_min/geometry/algorithm/CountCrossings.h>

#	include <iostream>

#	include <CGAL/Simple_cartesian.h>

namespace ogdf {

void GeometricVertexInsertion::call(GraphAttributes& GA) {
	using _Graph = internal::gcm::graph::OGDFGraphWrapper;
	using _Kernel = CGAL::Simple_cartesian<double>; //todo
	using Drawing = internal::gcm::graph::GeometricDrawing<_Kernel, _Graph>;

	_Graph _g(g);
	Drawing d(_g);

	internal::gcm::graph::ogdf_attributes_to_geometric_drawing(GA, d);

	std::mt19937_64 gen(ogdf::randomSeed());

	std::uniform_real_distribution<double> dist(-1, 1);

	std::vector<std::unique_ptr<Graph::HiddenEdgeSet>> hidden_edges;
	std::vector<node> hidden_vertex;

	for (auto pp : *m_vertex_order) {
		if (internal::gcm::geometry::count_crossings(d) > 0) {
			std::unique_ptr<Graph::HiddenEdgeSet> u(new Graph::HiddenEdgeSet(g));

			hidden_edges.push_back(std::move(u));
			hidden_vertex.push_back(pp);
			for (auto e : _g.edges(pp)) {
				hidden_edges.back()->hide(e);
			}
		}
	}

	if (m_initial_layout_module) {
		(*m_initial_layout_module)(GA);
	}

	for (int i = hidden_vertex.size() - 1; i >= 0; --i) {
		node v = hidden_vertex[i];
		hidden_edges[i]->restore();

		auto p = (*m_pos)(GA, v);
		GA.x(v) = p.m_x;
		GA.y(v) = p.m_y;
	}
}

}
#else
namespace ogdf {

void GeometricVertexInsertion::call(GraphAttributes&) {
	OGDF_THROW_PARAM(LibraryNotSupportedException, LibraryNotSupportedCode::Cgal);
}

}

#endif
