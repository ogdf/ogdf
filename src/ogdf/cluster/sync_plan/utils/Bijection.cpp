/** \file
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
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/basic.h>
#include <ogdf/cluster/sync_plan/basic/GraphUtils.h>
#include <ogdf/cluster/sync_plan/utils/Bijection.h>

#include <functional>
#include <utility>

namespace ogdf::sync_plan {

PipeBijRange getPipeBijection(node u, node v) {
	OGDF_ASSERT(u->degree() == v->degree());
	return PipeBijRange {PipeBijIterator {u->adjEntries.begin(), v->adjEntries.rbegin()},
			PipeBijIterator {u->adjEntries.end(), v->adjEntries.rend()}};
}

void getPipeBijection(node u, node v, PipeBij& out) {
	OGDF_ASSERT(u->graphOf() == v->graphOf());
	OGDF_ASSERT(u->degree() == v->degree());
	auto v_adj_it = v->adjEntries.rbegin();
	for (adjEntry u_adj : u->adjEntries) {
		OGDF_ASSERT(v_adj_it != v->adjEntries.rend());
		out.emplaceBack(u_adj, (*v_adj_it));
		v_adj_it++;
	}
	OGDF_ASSERT(v_adj_it == v->adjEntries.rend());
}

void getPipeBijection(node u, node v, AdjEntryArray<adjEntry>& out) {
	OGDF_ASSERT(u->graphOf() == v->graphOf());
	OGDF_ASSERT(u->graphOf() == out.graphOf());
	OGDF_ASSERT(u->degree() == v->degree());
	auto v_adj_it = v->adjEntries.rbegin();
	for (adjEntry u_adj : u->adjEntries) {
		OGDF_ASSERT(v_adj_it != v->adjEntries.rend());
		out[u_adj] = *v_adj_it;
		out[*v_adj_it] = u_adj;
		v_adj_it++;
	}
	OGDF_ASSERT(v_adj_it == v->adjEntries.rend());
}

void getPipeBijection(node u, node v, EdgeArray<edge>& out) {
	OGDF_ASSERT(u->graphOf() == v->graphOf());
	OGDF_ASSERT(u->graphOf() == out.graphOf());
	OGDF_ASSERT(u->degree() == v->degree());
	auto v_adj_it = v->adjEntries.rbegin();
	for (adjEntry u_adj : u->adjEntries) {
		OGDF_ASSERT(v_adj_it != v->adjEntries.rend());
		out[u_adj->theEdge()] = (*v_adj_it)->theEdge();
		v_adj_it++;
	}
	OGDF_ASSERT(v_adj_it == v->adjEntries.rend());
}

void getFrozenPipeBijection(node u, node v, FrozenPipeBij& bij) {
	OGDF_ASSERT(u->degree() == v->degree());
	auto v_adj_it = v->adjEntries.rbegin();
	for (adjEntry u_adj : u->adjEntries) {
		OGDF_ASSERT(v_adj_it != v->adjEntries.rend());
		bij.emplaceBack(u_adj->theEdge()->index(), (*v_adj_it)->theEdge()->index());
		++v_adj_it;
	}
	OGDF_ASSERT(v_adj_it == v->adjEntries.rend());
}

void freezePipeBijection(const PipeBij& in, FrozenPipeBij& out) {
	for (const PipeBijPair& pair : in) {
		out.emplaceBack(pair.first->theEdge()->index(), pair.second->theEdge()->index());
	}
}

std::pair<node, node> split(Graph& G, PipeBij& bij, const EdgeArray<edge>* new_edges,
		const EdgeArray<bool>* reverse_edges, node src, node tgt) {
	if (src == nullptr) {
		src = G.newNode();
	}
	if (tgt == nullptr) {
		tgt = G.newNode();
	}
	for (auto& pair : bij) {
		OGDF_ASSERT(pair.second == nullptr);
		edge new_edge = new_edges == nullptr ? nullptr : (*new_edges)[pair.first];
		pair.second = splitEdge(G, pair.first, src, tgt, new_edge);
		if (reverse_edges != nullptr && (*reverse_edges)[pair.first]) {
			G.reverseEdge(pair.second->theEdge());
		}
	}
	G.reverseAdjEdges(tgt);
	return std::pair<node, node>(src, tgt);
}

void join(Graph& G, node u, node v, PipeBij& bij, List<bool>* reverse_v) {
	join(
			G, u, v, bij, [&G](node n) { G.delNode(n); }, [&G](edge e) { G.delEdge(e); }, reverse_v);
}

void join(Graph& G, node u, node v, PipeBij& bij, const std::function<void(node)>& deleteNode,
		const std::function<void(edge)>& deleteEdge, List<bool>* reverse_v) {
	OGDF_ASSERT(u->degree() == bij.size());
	OGDF_ASSERT(v->degree() == bij.size());
	for (auto& pair : bij) {
		adjEntry f = pair.first->twin();
		bool rev = joinEdge(G, pair.first, pair.second, u, v, deleteEdge);
		if (reverse_v) {
			reverse_v->pushBack(rev);
		}
		pair.first = f;
		pair.second = nullptr;
	}
	OGDF_ASSERT(u->degree() == 0);
	OGDF_ASSERT(v->degree() == 0);
	deleteNode(u);
	deleteNode(v);
}

}
