/** \file
 * \brief Implementation of the SyncPlan::simplify operation and its UndoOperation.
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
#include <ogdf/basic/Logger.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/pctree/NodePCRotation.h>
#include <ogdf/cluster/sync_plan/PMatching.h>
#include <ogdf/cluster/sync_plan/SyncPlan.h>
#include <ogdf/cluster/sync_plan/SyncPlanComponents.h>
#include <ogdf/cluster/sync_plan/SyncPlan_operation/Simplify.h>
#include <ogdf/cluster/sync_plan/basic/GraphUtils.h>
#include <ogdf/cluster/sync_plan/basic/Iterators.h>
#include <ogdf/cluster/sync_plan/utils/Bijection.h>
#include <ogdf/cluster/sync_plan/utils/Logging.h>

#include <algorithm>
#include <functional>
#include <ostream>

namespace ogdf::pc_tree {
class PCNode;
} // namespace ogdf::pc_tree

using namespace ogdf::pc_tree;
using namespace ogdf::sync_plan::internal;

namespace ogdf::sync_plan {
using internal::operator<<;

UndoSimplify::UndoSimplify(const List<SimplifyMapping>& in_bij, node u2, node u, node v, node v2)
	: u2_idx(u2->index())
	, u_idx(u->index())
	, v_idx(v->index())
	, v2_idx(v2 == nullptr ? -1 : v2->index()) {
	for (const SimplifyMapping& pair : in_bij) {
		FrozenSimplifyMapping& entry =
				*bij.emplaceBack(pair.u2_adj->theEdge()->index(), pair.u_adj->theEdge()->index());
		for (adjEntry adj : pair.v_adj) {
			entry.v_adj.pushBack(adj->theEdge()->index());
		}
		OGDF_ASSERT(!entry.v_adj.empty());
		if (v2 != nullptr) {
			entry.v2_adj = pair.v2_adj->theEdge()->index();
		}
	}
}

void UndoSimplify::undo(SyncPlan& pq) {
	// SYNCPLAN_PROFILE_START("undo-simplify")
	node u2 = pq.nodeFromIndex(u2_idx), u = pq.nodeFromIndex(u_idx), v = pq.nodeFromIndex(v_idx),
		 v2 = nullptr;
	pq.log.lout(Logger::Level::Minor) << "Current orders:" << std::endl
									  << "u2: " << printIncidentEdges(u2->adjEntries) << std::endl
									  << "u:  " << printIncidentEdges(u->adjEntries) << std::endl
									  << "v:  " << printIncidentEdges(v->adjEntries) << std::endl;
	if (v2_idx >= 0) {
		v2 = pq.nodeFromIndex(v2_idx);
		pq.log.lout(Logger::Level::Minor)
				<< "v2: " << printIncidentEdges(v2->adjEntries) << std::endl;
	}

	List<SimplifyMapping> bij_list;
	AdjEntryArray<SimplifyMapping*> bij_map(*pq.G, nullptr);
	for (const FrozenSimplifyMapping& pair : bij) {
		SimplifyMapping& entry = *bij_list.emplaceBack(pq.edgeFromIndex(pair.u2_adj)->getAdj(u2),
				pq.edgeFromIndex(pair.u_adj)->getAdj(u));
		for (int adj_idx : pair.v_adj) {
			adjEntry adj = pq.edgeFromIndex(adj_idx)->getAdj(v);
			entry.v_adj.pushBack(adj);
			OGDF_ASSERT(bij_map[adj] == nullptr);
			bij_map[adj] = &entry;
		}
		OGDF_ASSERT(!entry.v_adj.empty());
		if (v2 != nullptr) {
			entry.v2_adj = pq.edgeFromIndex(pair.v2_adj)->getAdj(v2);
		}
		OGDF_ASSERT(bij_map[entry.u2_adj] == nullptr);
		bij_map[entry.u2_adj] = &entry;
	}

	List<adjEntry> u_order;
	List<adjEntry> v_order;
	for (adjEntry u2_adj : u2->adjEntries) {
		SimplifyMapping* entry = bij_map[u2_adj];
		OGDF_ASSERT(entry != nullptr);
		OGDF_ASSERT(entry->u2_adj == u2_adj);
		u_order.pushFront(entry->u_adj);
		if (entry->v_adj.size() == 1) {
			v_order.pushBack(entry->v_adj.front());
		} else {
			List<adjEntry> v_sub_order;
			for (auto it = v->adjEntries.rbegin(); it != v->adjEntries.rend(); it++) {
				adjEntry adj = *it;
				if (bij_map[adj] == entry) {
					v_sub_order.pushFront(adj);
				} else if (bij_map[adj] != nullptr) {
					// ignore any edges to other blocks (they are collected later and can be moved arbitrarily)
					// but break if we find an edge belonging to another bridge of this block
					break;
				}
			}
			if (v_sub_order.size() < entry->v_adj.size()) {
				for (adjEntry adj : v->adjEntries) {
					if (bij_map[adj] == entry) {
						v_sub_order.pushBack(adj);
					} else if (bij_map[adj] != nullptr && !v_sub_order.empty()) {
						// again ignore any edges to other blocks
						// but only break on an edge belonging to another bridge of this block if
						// they interrupt a continuous sequence of v edges
						break;
					}
				}
			}
			OGDF_ASSERT(v_sub_order.size() == entry->v_adj.size());
			v_order.conc(v_sub_order);
		}
	}
	pq.G->sort(u, u_order);
	if (u2 == v) {
		OGDF_ASSERT(compareCyclicOrder(v, v_order) == OrderComp::SAME);
	} else {
		if (v2 == nullptr && v->degree() != v_order.size()) {
			// v was a cut vertex, which is only allowed in the terminal case
			// we now need to collect the missing edges leading to other bicon components
			for (adjEntry adj : v->adjEntries) {
				if (bij_map[adj] == nullptr) {
					v_order.pushBack(adj);
				}
			}
		}
		pq.G->sort(v, v_order);
	}
	pq.log.lout(Logger::Level::Minor) << "New orders:" << std::endl
									  << "u2: " << printIncidentEdges(u2->adjEntries) << std::endl
									  << "u:  " << printIncidentEdges(u->adjEntries) << std::endl
									  << "v:  " << printIncidentEdges(v->adjEntries) << std::endl;
	if (v2) {
		pq.log.lout(Logger::Level::Minor)
				<< "v2: " << printIncidentEdges(v2->adjEntries) << std::endl;
	}

#ifdef OGDF_DEBUG
	if (v2 != nullptr) {
		List<adjEntry> v2_order;
		for (adjEntry u2_adj : u2->adjEntries) {
			v2_order.pushFront(bij_map[u2_adj]->v2_adj);
		}
		OGDF_ASSERT(compareCyclicOrder(v2, v2_order) == OrderComp::SAME);
	}
#endif

	if (v2 != nullptr) {
		pq.matchings.removeMatching(v2, u2);
		pq.matchings.matchNodes(v, v2);
	}
	pq.matchings.matchNodes(u2, u);
	// SYNCPLAN_PROFILE_STOP("undo-simplify")
}

class UndoSimplifyToroidal : public SyncPlan::UndoOperation {
	int u_idx, v_idx, u_first_adj_idx, v_last_adj_idx;
#ifdef OGDF_DEBUG
	FrozenPipeBij bij;
#endif

public:
	UndoSimplifyToroidal(node u, node v)
		: u_idx(u->index())
		, v_idx(v->index())
		, u_first_adj_idx(u->adjEntries.head()->theEdge()->index())
		, v_last_adj_idx(v->adjEntries.tail()->theEdge()->index()) {
#ifdef OGDF_DEBUG
		getFrozenPipeBijection(u, v, bij);
#endif
	}

	void undo(SyncPlan& pq) override {
		node u = pq.nodeFromIndex(u_idx);
		node v = pq.nodeFromIndex(v_idx);
		edge ue = pq.edgeFromIndex(u_first_adj_idx);
		edge ve = pq.edgeFromIndex(v_last_adj_idx);
		moveAdjToFront(*pq.G, ue->getAdj(u));
		moveAdjToBack(*pq.G, ve->getAdj(v));
		pq.matchings.matchNodes(u, v);
#ifdef OGDF_DEBUG
		pq.verifyPipeBijection(u, v, bij);
#endif
	}

	std::ostream& print(std::ostream& os) const override {
		return os << "UndoSimplifyToroidal(u_idx=" << u_idx << ", v_idx=" << v_idx
				  << ", u_first_adj_idx=" << u_first_adj_idx
				  << ", v_last_adj_idx=" << v_last_adj_idx << ")";
	}
};

SyncPlan::Result SyncPlan::simplify(node u, const NodePCRotation* pc) {
	OGDF_ASSERT(matchings.isMatchedPVertex(u));
	OGDF_ASSERT(!components.isCutVertex(u));
	if (!pc->isTrivial()) {
		return SyncPlan::Result::NOT_APPLICABLE;
	}
	OGDF_ASSERT(pc->getNode() == u);
	node v = pc->getTrivialPartnerPole();
	OGDF_ASSERT(v != nullptr);
	node u2 = matchings.getTwin(u);
	node v2 = matchings.getTwinOrNull(v);
	bool degree_mismatch = v->degree() != u->degree();
	log.lout(Logger::Level::High)
			<< "SIMPLIFY " << (v2 == nullptr ? "TERMINAL" : (v2 == u ? "TOROIDAL" : "TRANSITIVE"))
			<< " (" << (components.isCutVertex(v) ? "v is cut, " : "")
			<< (degree_mismatch ? "degree mismatch" : "degrees match") << ")"
			<< " u'=" << fmtPQNode(u2) << " < - > u=" << fmtPQNode(u) << " ==="
			<< " v=" << fmtPQNode(v) << " < - > v'=" << fmtPQNode(v2) << std::endl;

#ifdef SYNCPLAN_OPSTATS
	tp start = tpc::now();
	printOPStatsStart(matchings.getPipe(u),
			v2 == nullptr ? Operation::SIMPLIFY_TERMINAL
						  : (v2 == u ? Operation::SIMPLIFY_TOROIDAL : Operation::SIMPLIFY_TRANSITIVE),
			pc);
	stats_out << "\"pole_matched\":" << matchings.isMatchedPVertex(v)
			  << ",\"pole_cv\":" << components.isCutVertex(v) << ",\"pole_degree\":" << v->degree()
			  << ",";
	if (v2 != nullptr) {
		stats_out << "\"pole_twin_cv\":" << components.isCutVertex(v2)
				  << ",\"pole_twin_blocks\":" << components.biconnectedComponent(v2)->degree()
				  << ",\"pole_twin_bicon_size\":"
				  << components.bcSize(components.biconnectedComponent(v2))
				  << ",\"pole_twin_bc_id\":" << components.biconnectedId(v2)
				  << ",\"pole_twin_cc_id\":" << components.connectedId(v2) << ",";
	}
#endif

	if (v2 != nullptr) {
		int old_prio = matchings.getPipe(v)->pipe_priority;
		if (components.isCutVertex(v)) {
			log.lout() << "Parallel Pole v is a matched cut-vertex, can't Simplify" << std::endl;
			matchings.makePriority(v);
			log.lout() << "Increased priority of pipe (v, v') from " << old_prio << " to "
					   << matchings.getPipe(v)->pipe_priority << std::endl;
#ifdef SYNCPLAN_OPSTATS
			printOPStatsEnd(false, 0);
#endif
			return SyncPlan::Result::NOT_APPLICABLE;
		}
		if (degree_mismatch) {
			log.lout()
					<< "Parallel Pole v is matched but not trivial, can't Simplify before PropagatePQ(v')"
					<< std::endl;
			matchings.makePriority(v);
			log.lout() << "Increased priority of pipe (v, v') from " << old_prio << " to "
					   << matchings.getPipe(v)->pipe_priority << std::endl;
#ifdef SYNCPLAN_OPSTATS
			printOPStatsEnd(false, 0);
#endif
			return SyncPlan::Result::NOT_APPLICABLE;
		}
	}

#ifdef OGDF_DEBUG
	{
		OGDF_ASSERT(components.isCutVertex(v)
				|| components.biconnectedId(v) == components.biconnectedId(u));
		// either u and v are block vertices of the same bicon component, or v is a cut vertex
		// and we only care about the edges in the bicon component of u
		BiconnectedIsolation iso(components, components.biconnectedComponent(u));
		NodePCRotation v_pc(*G, v, true);
		validatePartnerPCTree(pc, &v_pc);
	}
#endif

	// SYNCPLAN_PROFILE_START("simplify")
	List<SimplifyMapping> bij_list;
	EdgeArray<SimplifyMapping*> bij_map;
	if (v2 != nullptr) {
		bij_map.init(*G, nullptr);
	}

	EdgeArray<PCNode*> leaf_for_u_inc_edge;
	if (pc->knowsPartnerEdges()) {
		leaf_for_u_inc_edge.init(*G, nullptr);
		pc->generateLeafForIncidentEdgeMapping(leaf_for_u_inc_edge);
	}
	EdgeSet visited(*G);

	// SYNCPLAN_PROFILE_START("simplify-bondmap")
	for (const PipeBijPair& pair : matchings.getIncidentEdgeBijection(u2)) {
		SimplifyMapping* entry = &(*bij_list.emplaceBack(pair.first, pair.second));

		if (pc->knowsPartnerEdges()) {
			const List<edge>& found = pc->getPartnerEdgesForLeaf(leaf_for_u_inc_edge[pair.second]);
			for (edge e : found) {
				OGDF_ASSERT(!visited.isMember(e));
				visited.insert(e);
				entry->v_adj.pushBack(e->getAdj(v));
			}
			OGDF_ASSERT(compareWithExhaustiveNodeDFS(*G, pair.second, v, visited, entry->v_adj));
		} else {
			exhaustiveNodeDFS(*G, pair.second, v, visited, entry->v_adj);
		}

		if (v2 != nullptr) {
			for (adjEntry adj : entry->v_adj) {
				bij_map[adj] = entry;
			}
			if (v2 != u) {
				bij_map[entry->u2_adj] = entry;
			}
		}
	}
	OGDF_ASSERT(validateCollectedAdjs(v, u, bij_list, visited, components));
	// SYNCPLAN_PROFILE_STOP("simplify-bondmap")

	if (v2 == u) {
		OGDF_ASSERT(u2 == v);

		List<List<adjEntry>> cycles;
		AdjEntryArray<adjEntry> pipe_bij(*G, nullptr);
		matchings.getIncidentEdgeBijection(u, pipe_bij);
		AdjEntryArray<adjEntry> mapping(*G, nullptr);
		for (adjEntry adj : v->adjEntries) {
			mapping[adj] = pipe_bij[bij_map[adj]->u_adj];
		}

		findCycles(*G, v, mapping, cycles, log.lout(Logger::Level::Minor) << "\t");
		int cycleLength = cycles.front().size();
		longestSimplifyToroidalCycle = max(cycleLength, longestSimplifyToroidalCycle);
		log.lout(Logger::Level::Medium) << "\tToroidal cycle length: " << cycleLength << std::endl;
#ifdef SYNCPLAN_OPSTATS
		stats_out << "\"cycle_len\":" << cycleLength << ",";
#endif
		for (auto& c : cycles) {
			if (c.size() != cycleLength) {
				// SYNCPLAN_PROFILE_STOP("simplify")
#ifdef SYNCPLAN_OPSTATS
				printOPStatsEnd(false, dur_ns(tpc::now() - start));
#endif
				return SyncPlan::Result::INVALID_INSTANCE;
			}
		}

		matchings.removeMatching(u, v);
		formatNode(u);
		formatNode(v);
	} else if (v2 != nullptr) {
		// SYNCPLAN_PROFILE_START("simplify-trans")
		for (const PipeBijPair& pair : matchings.getIncidentEdgeBijection(v)) {
			OGDF_ASSERT(bij_map[pair.first] != nullptr);
			bij_map[pair.first]->v2_adj = pair.second;
		}
		List<adjEntry> v2_adjs;
		for (adjEntry u2_adj : u2->adjEntries) {
			OGDF_ASSERT(bij_map[u2_adj] != nullptr);
			v2_adjs.pushFront(bij_map[u2_adj]->v2_adj);
		}
		G->sort(v2, v2_adjs);

		matchings.removeMatching(u, u2);
		matchings.removeMatching(v, v2);
		matchings.matchNodes(u2, v2);
		formatNode(u);
		formatNode(u2);
		formatNode(v);
		formatNode(v2);
		// SYNCPLAN_PROFILE_STOP("simplify-trans")
	} else {
		matchings.removeMatching(u, u2);
		formatNode(u);
		formatNode(u2);
	}

	if (v2 != u) {
		pushUndoOperation(new UndoSimplify(bij_list, u2, u, v, v2));
	} else {
		pushUndoOperation(new UndoSimplifyToroidal(u, v));
	}
#ifdef SYNCPLAN_OPSTATS
	printOPStatsEnd(true, dur_ns(tpc::now() - start));
#endif
	// SYNCPLAN_PROFILE_STOP("simplify")
	return SyncPlan::Result::SUCCESS;
}

}
