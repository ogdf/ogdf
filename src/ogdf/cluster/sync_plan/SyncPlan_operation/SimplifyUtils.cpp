/** \file
 * \brief Utilities for SyncPlan::simplify and its UndoOperation.
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
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/GraphSets.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/pctree/NodePCRotation.h>
#include <ogdf/basic/pctree/PCNode.h>
#include <ogdf/basic/pctree/util/IntrusiveList.h>
#include <ogdf/cluster/sync_plan/SyncPlan.h>
#include <ogdf/cluster/sync_plan/SyncPlanComponents.h>
#include <ogdf/cluster/sync_plan/SyncPlan_operation/Simplify.h>

#include <functional>
#include <ostream>

using ogdf::pc_tree::PCNode;

namespace ogdf::sync_plan::internal {

SimplifyMapping::SimplifyMapping(adjEntry u2Adj, adjEntry uAdj, adjEntry vAdj, adjEntry v2Adj)
	: u2_adj(u2Adj), u_adj(uAdj), v2_adj(v2Adj) {
	if (vAdj != nullptr) {
		v_adj.pushBack(vAdj);
	}
}

std::ostream& operator<<(std::ostream& os, const SimplifyMapping& mapping) {
	os << "u2_e" << (mapping.u2_adj == nullptr ? -1 : mapping.u2_adj->theEdge()->index()) << " = u_e"
	   << (mapping.u_adj == nullptr ? -1 : mapping.u_adj->theEdge()->index()) << " = v_e{";
	for (adjEntry adj : mapping.v_adj) {
		os << adj->theEdge()->index() << " ";
	}
	return os << "} = v2_e" << (mapping.v2_adj == nullptr ? -1 : mapping.v2_adj->theEdge()->index());
}

FrozenSimplifyMapping::FrozenSimplifyMapping(int u2Adj, int uAdj, int vAdj, int v2Adj)
	: u2_adj(u2Adj), u_adj(uAdj), v2_adj(v2Adj) {
	if (vAdj >= 0) {
		v_adj.pushBack(vAdj);
	}
}

std::ostream& operator<<(std::ostream& os, const FrozenSimplifyMapping& mapping) {
	return os << "u2_e" << mapping.u2_adj << " = u_e" << mapping.u_adj << " = v_e{" << mapping.v_adj
			  << "} = v2_e" << mapping.v2_adj;
}

std::ostream& UndoSimplify::print(std::ostream& os) const {
	os << "UndoSimplify"
	   << (v2_idx >= 0 ? (v2_idx == u2_idx ? "Toroidal" : "Transitive") : "Terminal")
	   << "(u2=" << u2_idx << ", u=" << u_idx << ", v=" << v_idx << ", v2=" << v2_idx << ", bij=";
	for (const FrozenSimplifyMapping& pair : bij) {
		os << std::endl << "\t" << pair;
	}
	return os << ")";
}

int findCycles(Graph& G, node u, const std::function<adjEntry(adjEntry)>& mapping,
		List<List<adjEntry>>& cycles, std::ostream& log) {
	AdjEntryArray<bool> checked(G, false);
	int sum = 0;
	log << "Cycles:";
	for (adjEntry adj : u->adjEntries) {
		adjEntry cur_adj = adj;
		if (checked[cur_adj]) {
			continue;
		}
		auto cycle = cycles.emplaceBack();
		log << "[";
		for (int i = 0; i <= u->degree(); i++) {
			OGDF_ASSERT(i < u->degree());
			OGDF_ASSERT(!checked[cur_adj]);
			checked[cur_adj] = true;
			(*cycle).pushBack(cur_adj);
			log << cur_adj->theEdge()->index();
			cur_adj = mapping(cur_adj);
			OGDF_ASSERT(cur_adj->theNode() == u);
			if (cur_adj == adj) {
				OGDF_ASSERT((*cycle).size() == (i + 1));
				sum += i + 1;
				log << "]#" << (i + 1) << " ";
				break;
			} else {
				log << ", ";
			}
		}
		if (cycles.size() > 1 && (*cycle).size() != cycles.front().size()) {
			log << "MISMATCH!" << std::endl;
			return -1;
		}
	}
	log << std::endl;
	return sum;
}

adjEntry continueNodeDFS(Graph& G, adjEntry u_adj, node v, EdgeSet& visited,
		List<adjEntry>& dfs_stack) {
	while (!dfs_stack.empty()) {
		OGDF_ASSERT(dfs_stack.back()->twinNode() != u_adj->theNode());
		OGDF_ASSERT(dfs_stack.back()->theNode() != v);
		OGDF_ASSERT(visited.isMember(dfs_stack.back()->theEdge()));
		if (dfs_stack.back()->twinNode() == v) {
			return dfs_stack.back();
		}

		adjEntry next = nullptr;
		for (adjEntry twin_adj : dfs_stack.back()->twinNode()->adjEntries) {
			if (!visited.isMember(twin_adj->theEdge())) {
				next = twin_adj;
				break;
			}
		}

		while (next == nullptr) {
			adjEntry failed = dfs_stack.popBackRet();
			if (dfs_stack.empty()) {
				return nullptr;
			}
			next = failed->cyclicSucc();
			while (next != failed) {
				if (!visited.isMember(next->theEdge())) {
					break;
				} else {
					next = next->cyclicSucc();
				}
			}
			if (next == failed) {
				next = nullptr;
			}
		}

		visited.insert(next->theEdge());
		dfs_stack.pushBack(next);
	}
	return nullptr;
}

adjEntry startNodeDFS(Graph& G, adjEntry u_adj, node v, EdgeSet& visited, List<adjEntry>& dfs_stack) {
	dfs_stack.pushBack(u_adj);
	visited.insert(u_adj->theEdge());
	adjEntry found = continueNodeDFS(G, u_adj, v, visited, dfs_stack);
	OGDF_ASSERT(found != nullptr);
	OGDF_ASSERT(dfs_stack.front() == u_adj);
	OGDF_ASSERT(dfs_stack.back()->twinNode() == v);
	return found;
}

int exhaustiveNodeDFS(Graph& G, adjEntry u_adj, node v, EdgeSet& visited, List<adjEntry>& out) {
	List<adjEntry> queue;
	adjEntry found = startNodeDFS(G, u_adj, v, visited, queue);
	out.pushBack(found->twin());
	if (u_adj->theNode()->degree() != v->degree()) { // degree_mismatch
		while (true) {
			OGDF_ASSERT(out.size() <= v->degree());
			queue.popBack();
			found = continueNodeDFS(G, u_adj, v, visited, queue);
			if (found == nullptr) {
				break;
			}
			out.pushBack(found->twin());
		}
	}
	return out.size();
}

#ifdef OGDF_DEBUG

bool compareWithExhaustiveNodeDFS(Graph& G, adjEntry u_adj, node v, const EdgeSet& visited,
		const List<adjEntry>& found) {
	List<adjEntry> reference_list;
	EdgeSet reference_visited(G);
	exhaustiveNodeDFS(G, u_adj, v, reference_visited, reference_list);
	OGDF_ASSERT(reference_list.size() == found.size());
	for (adjEntry adj : found) {
		OGDF_ASSERT(reference_visited.isMember(adj->theEdge()));
	}
	for (adjEntry adj : reference_list) {
		OGDF_ASSERT(visited.isMember(adj->theEdge()));
	}
	return true;
}

void validatePartnerPCTree(const NodePCRotation* u_pc, const NodePCRotation* v_pc) {
	OGDF_ASSERT(u_pc->getTrivialPartnerPole() == v_pc->getNode());
	if (u_pc->getNode()->degree() == v_pc->getNode()->degree()) { // degrees match
		OGDF_ASSERT(v_pc->getTrivialPartnerPole() == u_pc->getNode());
		OGDF_ASSERT(v_pc->getLeaves().size() == u_pc->getLeaves().size());

		if (u_pc->knowsPartnerEdges()) {
			EdgeArray<PCNode*> u_leaf_for_inc_edge(*u_pc->getGraph(), nullptr);
			u_pc->generateLeafForIncidentEdgeMapping(u_leaf_for_inc_edge);
			for (auto v_leaf : v_pc->getLeaves()) {
				List<edge> u_edges = v_pc->getPartnerEdgesForLeaf(v_leaf);
				OGDF_ASSERT(u_edges.size() == 1);
				PCNode* u_leaf = u_leaf_for_inc_edge[u_edges.front()];
				const List<edge>& v_edges = u_pc->getPartnerEdgesForLeaf(u_leaf);
				OGDF_ASSERT(v_edges.size() == 1);
				OGDF_ASSERT(v_edges.front() == v_pc->getIncidentEdgeForLeaf(v_leaf));
			}
		}
	} else {
		OGDF_ASSERT(u_pc->getNode()->degree() < v_pc->getNode()->degree());
		// room for improvement: we could validate the incident edge for leaf mapping here,
		// but that is tedious and always worked alright up to now
	}
}

bool validateCollectedAdjs(node v, node u, List<SimplifyMapping>& bij_list, EdgeSet& visited,
		SyncPlanComponents& components) {
	int collected_v_adjs = 0;
	for (auto& entry : bij_list) {
		for (adjEntry adj : entry.v_adj) {
			OGDF_ASSERT(adj->theNode() == v);
			OGDF_ASSERT(visited.isMember(adj->theEdge()));
			collected_v_adjs++;
		}
	}
	OGDF_ASSERT(visited.size() >= collected_v_adjs);
	if (!components.isCutVertex(v)) {
		OGDF_ASSERT(collected_v_adjs == v->degree());
		for (adjEntry adj : v->adjEntries) {
			OGDF_ASSERT(visited.isMember(adj->theEdge()));
		}
	} else {
		OGDF_ASSERT(collected_v_adjs < v->degree());
		int v_adjs_in_u_bc = 0;
		node u_bc = components.biconnectedComponent(u);
		for (adjEntry adj : v->adjEntries) {
			node neigh_bc = components.biconnectedComponent(adj->twinNode());
			bool adj_in_bc = neigh_bc == u_bc;
			if (components.isCutVertex(adj->twinNode())) {
				adj_in_bc = components.bcTree().searchEdge(neigh_bc, u_bc) != nullptr;
			}
			if (adj_in_bc) {
				v_adjs_in_u_bc++;
				OGDF_ASSERT(visited.isMember(adj->theEdge()));
			} else {
				OGDF_ASSERT(!visited.isMember(adj->theEdge()));
			}
		}
		OGDF_ASSERT(v_adjs_in_u_bc == collected_v_adjs);
	}
	return true;
}

#endif // OGDF_DEBUG

}
