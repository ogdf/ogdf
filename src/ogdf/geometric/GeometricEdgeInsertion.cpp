/** \file
 * \brief In the Geometric Edge Insertiont Approach we iteratively insert a given set of edges into a drawing of a graph.
 * In each insertion the end vertices are moved to their optimal position.
 * The optimal position can be for example a position that minimizes the number of crossings for this vertex.
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

#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/fileformats/SvgPrinter.h>
#include <ogdf/geometric/GeometricEdgeInsertion.h>

namespace ogdf {


void GeometricEdgeInsertion::call(GraphAttributes& GA) {
	Graph::HiddenEdgeSet hidden_edges(g);
	for (auto e : *m_edge_set) {
		hidden_edges.hide(e);
	}

	if (m_initial_layout_module) {
		(*m_initial_layout_module)(GA);
	}

	auto move = [&](node v) {
		auto p = (*m_pos)(GA, v);
		GA.x(v) = p.m_x;
		GA.y(v) = p.m_y;
	};

	for (auto e : *m_edge_set) {
		hidden_edges.restore(e);
		move(e->source());
		move(e->target());
	}
}


}
