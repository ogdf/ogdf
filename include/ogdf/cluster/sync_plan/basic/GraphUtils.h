/** \file
 * \brief Some missing utilities for working with Graphs, their embeddings, and cyclic orders. TODO should be moved to a central location (some maybe part of the Graph class?).
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

//! Change one endpoint of an edge, no matter its direction.
OGDF_EXPORT void moveEnd(Graph& G, edge e, node keep_end, node new_end);

//! Change one endpoint of an edge observing a certain embedding, no matter its direction.
OGDF_EXPORT void moveEnd(Graph& G, adjEntry keep_adj, adjEntry new_adj,
		Direction dir = Direction::after);

//! Split an edge, moving the two new middle endpoints to some other vertices.
OGDF_EXPORT edge splitEdge(Graph& G, edge old_edge, node new_adj_to_source, node new_adj_to_target,
		edge new_edge = nullptr);

//! Split an edge, moving the two new middle endpoints to some other vertices observing a certain embedding.
OGDF_EXPORT adjEntry splitEdge(Graph& G, adjEntry adj, node new_adj_to_node, node new_adj_to_twin,
		edge new_edge = nullptr);

//! Join two edges into one, keeping the two given endpoints.
OGDF_EXPORT bool joinEdge(Graph& G, edge u_e, edge v_e, node u, node v);

//! Join two edges into one, keeping the two endpoints corresponding to the adjEntries.
OGDF_EXPORT bool joinEdge(Graph& G, adjEntry u_adj, adjEntry v_adj, node u, node v);

//! Join two edges into one, keeping the two given endpoints and using a custom function for deleting the old edge.
OGDF_EXPORT bool joinEdge(Graph& G, edge u_e, edge v_e, node u, node v,
		const std::function<void(edge)>& deleteEdge);

//! Join two edges into one, keeping the two endpoints corresponding to the adjEntries and using a custom function for deleting the old edge.
OGDF_EXPORT bool joinEdge(Graph& G, adjEntry u_adj, adjEntry v_adj, node u, node v,
		const std::function<void(edge)>& deleteEdge);

//! Check that one vertex is the centre of star while the other is one of its rays.
OGDF_EXPORT void assertStarCentreAndRay(node centre, node ray);

//! Given a vertex that is either the centre or ray of a star, return the centre of the respective star.
OGDF_EXPORT node getCentreOfStar(node g_n);

enum class OrderComp { SAME, REVERSED, DIFFERENT };

//! Cyclically compare the rotation of a node with a given cyclic order.
OGDF_EXPORT OrderComp compareCyclicOrder(node n, List<adjEntry>& o, bool full_check = false);

//! Rotate a node to move a given adjEntry to the front of its list.
OGDF_EXPORT void moveAdjToFront(Graph& G, adjEntry f);

//! Rotate a node to move a given adjEntry to the back of its list.
OGDF_EXPORT void moveAdjToBack(Graph& G, adjEntry b);

}
