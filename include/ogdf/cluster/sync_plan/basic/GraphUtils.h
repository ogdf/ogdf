/** \file
 * \brief TODO Document
 *
 * \author Simon D. Fink <ogdf@niko.fink.bayern>
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

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/basic.h>

#include <functional>

namespace ogdf {
template<class E>
class List;
} // namespace ogdf

namespace ogdf {

OGDF_EXPORT void moveEnd(Graph& G, edge e, node keep_end, node new_end);

OGDF_EXPORT void moveEnd(Graph& G, adjEntry keep_adj, adjEntry new_adj,
		Direction dir = Direction::after);

OGDF_EXPORT edge splitEdge(Graph& G, edge old_edge, node new_adj_to_source, node new_adj_to_target,
		edge new_edge = nullptr);

OGDF_EXPORT adjEntry splitEdge(Graph& G, adjEntry adj, node new_adj_to_node, node new_adj_to_twin,
		edge new_edge = nullptr);

OGDF_EXPORT bool joinEdge(Graph& G, edge u_e, edge v_e, node u, node v);

OGDF_EXPORT bool joinEdge(Graph& G, adjEntry u_adj, adjEntry v_adj, node u, node v);

OGDF_EXPORT bool joinEdge(Graph& G, edge u_e, edge v_e, node u, node v,
		const std::function<void(edge)>& deleteEdge);

OGDF_EXPORT bool joinEdge(Graph& G, adjEntry u_adj, adjEntry v_adj, node u, node v,
		const std::function<void(edge)>& deleteEdge);

OGDF_EXPORT void assertStarCentreAndRay(node centre, node ray);

OGDF_EXPORT node getCentreOfStar(node g_n);

enum class OrderComp { SAME, REVERSED, DIFFERENT };

OGDF_EXPORT OrderComp compareCyclicOrder(node n, List<adjEntry>& o, bool full_check = false);

OGDF_EXPORT void moveAdjToFront(Graph& G, adjEntry f);

OGDF_EXPORT void moveAdjToBack(Graph& G, adjEntry b);

}
