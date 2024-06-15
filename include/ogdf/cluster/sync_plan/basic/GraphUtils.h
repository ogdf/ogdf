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
#include <ogdf/cluster/sync_plan/basic/Iterators.h>
#include <ogdf/cluster/sync_plan/utils/Bijection.h>

using namespace ogdf;

void moveEnd(Graph& G, edge e, node keep_end, node new_end);

void moveEnd(Graph& G, adjEntry keep_adj, adjEntry new_adj, Direction dir = Direction::after);

edge splitEdge(Graph& G, edge old_edge, node new_adj_to_source, node new_adj_to_target,
		int new_edge_idx = -1);

adjEntry splitEdge(Graph& G, adjEntry adj, node new_adj_to_node, node new_adj_to_twin,
		int new_edge_idx = -1);

bool joinEdge(Graph& G, edge u_e, edge v_e, node u, node v);

bool joinEdge(Graph& G, adjEntry u_adj, adjEntry v_adj, node u, node v);

std::pair<node, node> split(Graph& G, PipeBij& bij, const EdgeArray<int>* split_idcs = nullptr,
		const EdgeArray<bool>* split_reverse = nullptr, int src_idx = -1, int tgt_idx = -1);

void join(Graph& G, node u, node v, PipeBij& bij, List<bool>* reverse_v = nullptr);

void assertStarCentreAndRay(node centre, node ray);

node getCentreOfStar(node g_n);

enum class OrderComp { SAME, REVERSED, DIFFERENT };

OrderComp compareCyclicOrder(node n, List<adjEntry>& o, bool full_check = false);

void moveAdjToFront(Graph& G, adjEntry f);

void moveAdjToBack(Graph& G, adjEntry b);
