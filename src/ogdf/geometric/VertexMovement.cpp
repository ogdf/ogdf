/** \file
 * \brief In the VertexMovement Approach the vertices are moved one by one to its optimal position.
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
#include <ogdf/geometric/VertexMovement.h>

namespace ogdf {

VertexMovement::VertexMovement() {
	//nothing to do here
}

VertexMovement::~VertexMovement() {
	//nothing to do here
}

void VertexMovement::call(GraphAttributes& GA) {
	for (auto v : *m_vertex_order) {
		auto p = (*m_pos)(GA, v);

		GA.x(v) = p.m_x;
		GA.y(v) = p.m_y;
	}
}


}
