/** \file
 * \brief Implementation of the SyncPlan::propagate operation and its UndoOperation.
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
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/GraphSets.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/Logger.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/pctree/NodePCRotation.h>
#include <ogdf/basic/pctree/PCEnum.h>
#include <ogdf/basic/pctree/PCNode.h>
#include <ogdf/basic/pctree/PCRegistry.h>
#include <ogdf/basic/pctree/PCTreeIterators.h>
#include <ogdf/basic/pctree/util/FilteringBFS.h>
#include <ogdf/basic/pctree/util/IntrusiveList.h>
#include <ogdf/cluster/sync_plan/PMatching.h>
#include <ogdf/cluster/sync_plan/QPartitioning.h>
#include <ogdf/cluster/sync_plan/SyncPlan.h>
#include <ogdf/cluster/sync_plan/SyncPlanComponents.h>
#include <ogdf/cluster/sync_plan/basic/GraphUtils.h>
#include <ogdf/cluster/sync_plan/utils/Bijection.h>
#include <ogdf/cluster/sync_plan/utils/Logging.h>

#include <memory>
#include <sstream>
#include <string>
#include <utility>

using namespace ogdf::pc_tree;
using namespace ogdf::sync_plan::internal;

namespace ogdf::sync_plan {
using internal::operator<<;

void smoothLeafToEdge(Graph* G, adjEntry former_adj, node leaf) {
	OGDF_ASSERT(leaf->outdeg() == 1);
	OGDF_ASSERT(leaf->indeg() == 0);
	edge adj_edge = former_adj->theEdge();
	edge tree_edge = leaf->adjEntries.head()->theEdge();
	bool reverse = former_adj->isSource();
	if (reverse) {
		G->reverseEdge(adj_edge);
	}
	OGDF_ASSERT(!former_adj->isSource());
	G->moveTarget(adj_edge, leaf);
	G->unsplit(adj_edge, tree_edge);
	if (reverse) {
		G->reverseEdge(adj_edge);
	}
}

class UndoPropagate : public SyncPlan::UndoOperation {
public:
	List<int> pct_u, pct_v;
	int u_idx, v_idx;
	int u_first_adj_idx, v_last_adj_idx;
	int degree;
#ifdef OGDF_DEBUG
	FrozenPipeBij bij;
#endif

	UndoPropagate(node u, node v)
		: u_idx(u->index())
		, v_idx(v->index())
		, u_first_adj_idx(u->adjEntries.head()->theEdge()->index())
		, v_last_adj_idx(v->adjEntries.tail()->theEdge()->index())
		, degree(u->degree()) {
#ifdef OGDF_DEBUG
		getFrozenPipeBijection(u, v, bij);
#endif
	}

	void undo(SyncPlan& pq) override {
		// SYNCPLAN_PROFILE_START("undo-propagatePQ")
		OGDF_ASSERT(pct_u.size() == pct_v.size());
		pq.log.lout(Logger::Level::High) << "UNDO PROPAGATE PQ degree " << degree << " with "
										 << pct_u.size() << " tree nodes." << std::endl;
		Logger::Indent _(&pq.log);

		NodeArray<int> markers(*pq.G, 0);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
		for (int n : pct_u) {
			markers[n] = 1;
		}
		for (int n : pct_v) {
			OGDF_ASSERT(markers[n] == 0);
			markers[n] = 2;
		}
#pragma GCC diagnostic pop

		node u = pq.nodeFromIndex(u_idx);
		OGDF_ASSERT(pq.deletedNodes.remove(u));
		edge ue = pq.G->newEdge(u, pq.nodeFromIndex(pct_u.front()));
		pq.log.lout() << "Collapsing tree at u into " << pq.fmtPQNode(u, false) << " via edge #"
					  << ue->index() << " " << ue << "." << std::endl;
		int un = collapseTree(pq, u, 1, markers);
		moveAdjToFront(*pq.G, pq.edgeFromIndex(u_first_adj_idx)->getAdj(u));
		pq.log.lout(Logger::Level::Minor) << printIncidentEdges(u->adjEntries) << std::endl;
		OGDF_ASSERT(un == pct_u.size());
		OGDF_ASSERT(u->degree() == degree);

		node v = pq.nodeFromIndex(v_idx);
		OGDF_ASSERT(pq.deletedNodes.remove(v));
		edge ve = pq.G->newEdge(v, pq.nodeFromIndex(pct_v.front()));
		pq.log.lout() << "Collapsing tree at v into " << pq.fmtPQNode(v, false) << " via edge #"
					  << ve->index() << " " << ve << "." << std::endl;
		int vn = collapseTree(pq, v, 2, markers);
		moveAdjToBack(*pq.G, pq.edgeFromIndex(v_last_adj_idx)->getAdj(v));
		pq.log.lout(Logger::Level::Minor) << printIncidentEdges(v->adjEntries) << std::endl;
		OGDF_ASSERT(vn == pct_v.size());
		OGDF_ASSERT(v->degree() == degree);

		pq.matchings.matchNodes(u, v);
		pq.log.lout(Logger::Level::Minor) << pq.matchings.printBijection(u) << std::endl;
		pq.log.lout(Logger::Level::Medium) << "Revived matched nodes " << pq.fmtPQNode(u, false)
										   << " and " << pq.fmtPQNode(v, false) << "." << std::endl;
#ifdef OGDF_DEBUG
		pq.verifyPipeBijection(u, v, bij);
#endif
		// SYNCPLAN_PROFILE_STOP("undo-propagatePQ")
	}

	int collapseTree(SyncPlan& pq, node root, int mark, NodeArray<int>& markers) const {
		Logger::Indent _(&pq.log);
		int count = 0;
		adjEntry adj = root->adjEntries.head();
		while (adj != nullptr) {
			OGDF_ASSERT(adj->theNode() == root);
			node n = adj->twinNode();
			if (markers[n] != mark) {
				adj = adj->succ();
				continue;
			}
			pq.log.lout(Logger::Level::Medium) << "Contracting edge #" << adj->theEdge()->index()
											   << " " << adj << "." << std::endl;
			pq.log.lout(Logger::Level::Minor)
					<< "\t" << printIncidentEdges(n->adjEntries) << std::endl;

			if (pq.matchings.isMatchedPVertex(n)) {
				pq.matchings.removeMatching(n);
			} else if (pq.partitions.isQVertex(n)) {
				pq.partitions.releaseQVertex(n);
			}

			adjEntry pred = adj->pred();
			if (pred == adj) {
				pred = nullptr;
			}

			if (!adj->isSource()) {
				pq.G->reverseEdge(adj->theEdge());
			}
			node ret = pq.G->contract(adj->theEdge());
			OGDF_ASSERT(ret == root);

			if (pred == nullptr) {
				adj = root->adjEntries.head();
			} else {
				adj = pred->succ();
			}
			count++;
		}
#ifdef OGDF_HEAVY_DEBUG
		for (adjEntry a : root->adjEntries) {
			OGDF_ASSERT(markers[a->twinNode()] != mark);
		}
#endif
		if (pq.matchings.isMatchedPVertex(root)) {
			pq.matchings.removeMatching(root);
		} else if (pq.partitions.isQVertex(root)) {
			pq.partitions.releaseQVertex(root);
		}
		return count;
	}

	std::ostream& print(std::ostream& os) const override {
		return os << "UndoPropagate(u=" << u_idx << ", v=" << v_idx << ", degree=" << degree
				  << ", pct_u=" << pct_u << ", pct_v=" << pct_v << ")";
	}
};

SyncPlan::Result SyncPlan::propagatePQ(node u, NodePCRotation* pct, NodePCRotation* pct_v) {
	// SYNCPLAN_PROFILE_START("propagatePQ")
	OGDF_ASSERT(matchings.isMatchedPVertex(u));
	OGDF_ASSERT(!components.isCutVertex(u));
	OGDF_ASSERT(!pct->isTrivial());
	OGDF_ASSERT(pct->getNode() == u);
	node u_bc = components.biconnectedComponent(u);
	node v = matchings.getTwin(u);
	bool v_was_cut = components.isCutVertex(v);
	List<node> v_rays, make_wheels;
	log.lout(Logger::Level::High)
			<< "PROPAGATE PQ into " << (v_was_cut ? "cut" : "biconnected") << std::endl;
#ifdef SYNCPLAN_OPSTATS
	tp start = tpc::now();
	printOPStatsStart(matchings.getPipe(u),
			v_was_cut ? Operation::PROPAGATE_CUT : Operation::PROPAGATE_BICON, pct);
	int64_t v_pc_dur = 0;
#endif
	Logger::Indent _(&log);

	if (v_was_cut) {
		Result result = encapsulate(v);
		OGDF_ASSERT(result == SyncPlan::Result::SUCCESS);
		for (node n : FilteringBFS(*G, {v})) {
			if (n != v) {
				v_rays.pushBack(n);
			}
		}
		// v_rays may not contain duplicates (they would be contracted twice), BFS nicely prevents that

		// the encapsulation broke the incident edge for leaf mapping if u and v are adjacent, so fix the mapping
		for (PCNode* pctn : pct->allNodes()) {
			if (pctn->isLeaf()) {
				edge u_e = pct->getIncidentEdgeForLeaf(pctn);
				if (!u_e->isIncident(u)) {
					pct->setIncidentEdgeForLeaf(pctn,
							matchings.translateIncidentEdge(u_e->getAdj(v)->twin())->theEdge());
				}
			}
		}
	} else if (intersect_trees) {
		bool intersected;
#ifdef SYNCPLAN_OPSTATS
		tp v_pc_start = tpc::now();
#endif
		// SYNCPLAN_PROFILE_START("propagatePQ-makePCv")
		std::unique_ptr<NodePCRotation> computed_pct_v;
		if (pct_v == nullptr) {
			BiconnectedIsolation iso(components, components.biconnectedComponent(v));
			computed_pct_v = std::make_unique<NodePCRotation>(*G, v);
			pct_v = computed_pct_v.get();
		}

		log.lout() << "Intersecting with " << (computed_pct_v ? "computed " : "provided ")
				   << "PC-Tree of block vertex v: " << pct_v << std::endl;
#ifdef SYNCPLAN_OPSTATS
		stats_out << "\"v_p_nodes\":" << pct_v->getPNodeCount()
				  << ",\"v_c_nodes\":" << pct_v->getCNodeCount()
				  << ",\"v_p_node_degs\":" << sumPNodeDegrees(*pct_v) << ",";
#endif
		OGDF_ASSERT(pct->getLeafCount() == pct_v->getLeafCount());

		EdgeArray<edge> v_edge_to_u_edge(*G, nullptr);
		matchings.getIncidentEdgeBijection(v, v_edge_to_u_edge);
		EdgeArray<PCNode*> edge_to_pct_leaf(*G, nullptr);
		pct->generateLeafForIncidentEdgeMapping(edge_to_pct_leaf);
		PCTreeNodeArray<PCNode*> mapping(*pct_v, nullptr);
		for (auto leaf : pct_v->getLeaves()) {
			mapping[leaf] = edge_to_pct_leaf[v_edge_to_u_edge[pct_v->getIncidentEdgeForLeaf(leaf)]];
			OGDF_ASSERT(mapping[leaf] != nullptr);
		}

		intersected = pct->intersect(*pct_v, mapping);
		// SYNCPLAN_PROFILE_STOP("propagatePQ-makePCv")
#ifdef SYNCPLAN_OPSTATS
		v_pc_dur = dur_ns(tpc::now() - v_pc_start);
		stats_out << "\"v_pc_time_ns\":" << v_pc_dur << ",";
#endif

		if (!intersected) {
#ifdef SYNCPLAN_OPSTATS
			printOPStatsEnd(false, dur_ns(tpc::now() - start));
#endif
			log.lout() << "Intersecting failed!" << std::endl;
			// SYNCPLAN_PROFILE_STOP("propagatePQ")
			return SyncPlan::Result::INVALID_INSTANCE; // TODO we need to proceed with the insertion if we want to find a kuratowksi
		}
		log.lout() << "Intersected PC-Tree: " << *pct << std::endl;
#ifdef SYNCPLAN_OPSTATS
		stats_out << "\"i_p_nodes\":" << pct->getPNodeCount()
				  << ",\"i_c_nodes\":" << pct->getCNodeCount()
				  << ",\"i_p_node_degs\":" << sumPNodeDegrees(*pct) << ",";
#endif
	}
	// encapsulate changes v_bc and also potentially edges incident to u
	node v_bc = components.biconnectedComponent(v);
	auto undo_op = new UndoPropagate(u, v);

	Graph pcg;
	node inner_pcg_node = nullptr;
	NodeArray<PCNode*> pcg_to_pct(pcg, nullptr);
	PCTreeNodeArray<ogdf::node> pct_to_pcg(*pct, nullptr);
	pct->getTree(pcg, nullptr, pct_to_pcg);
	pcg.reverseAllEdges();
	for (PCNode* pctn : pct->allNodes()) {
		pcg_to_pct[pct_to_pcg[pctn]] = pctn;
	}

	NodeArray<node> pcg_to_u(pcg, nullptr);
	EdgeArray<edge> dummy(pcg, nullptr);
	G->insert(pcg, pcg_to_u, dummy);

	NodeArray<node> pcg_to_v(pcg, nullptr);
	dummy.fill(nullptr);
	G->insert(pcg, pcg_to_v, dummy);

	log.lout()
			<< "Inserted two trees with " << pcg.numberOfNodes()
			<< " nodes each, now rewiring the 2*" << pct->getLeafCount() << " leaves." << std::endl;
	log.lout(Logger::Level::Medium) << matchings.printBijection(u) << std::endl;

	AdjEntryArray<adjEntry> u_adj_to_v_adj(*G, nullptr);
	matchings.getIncidentEdgeBijection(u, u_adj_to_v_adj);
	matchings.removeMatching(u, v);

	for (node pcgn : pcg.nodes) {
		PCNode* pctn = pcg_to_pct[pcgn];
		OGDF_ASSERT(pctn != nullptr);
		if (pctn->isLeaf()) {
			OGDF_ASSERT(pcgn->degree() == 1);
			edge u_e = pct->getIncidentEdgeForLeaf(pctn);
			OGDF_ASSERT(u_e != nullptr);
			adjEntry u_adj = u_e->getAdj(u);
			adjEntry v_adj = u_adj_to_v_adj[u_adj];
			smoothLeafToEdge(G, u_adj, pcg_to_u[pcgn]);
			smoothLeafToEdge(G, v_adj, pcg_to_v[pcgn]);
		}
	}

	log.lout() << "Fixing the remaining 2*" << (pcg.numberOfNodes() - pct->getLeafCount())
			   << " inner nodes." << std::endl;
	for (node pcgn : pcg.nodes) {
		PCNode* pctn = pcg_to_pct[pcgn];
		if (pctn == nullptr || pctn->getNodeType() != PCNodeType::Leaf) {
			Logger::Indent __(&log);
			OGDF_ASSERT(pcgn->degree() > 2);
			inner_pcg_node = pcgn;
			components.nodeInserted(pcg_to_u[pcgn], u_bc);
			components.nodeInserted(pcg_to_v[pcgn], v_bc);
			undo_op->pct_u.pushBack(pcg_to_u[pcgn]->index());
			undo_op->pct_v.pushBack(pcg_to_v[pcgn]->index());
			bool fixed_rotation = pctn == nullptr || pctn->getNodeType() == PCNodeType::CNode
					|| pcgn->degree() <= 3;
			G->reverseAdjEdges(pcg_to_v[pcgn]);

			if (GA != nullptr) {
				std::stringstream ss;
				ss << "Prop " << (fixed_rotation ? "Q" : "P") << " " << pcg_to_u[pcgn]->index()
				   << " [" << u->index() << "]." << pcgn->index();
				GA->label(pcg_to_u[pcgn]) = ss.str();
				ss.str("");
				ss << "Prop " << (fixed_rotation ? "Q" : "P") << " " << pcg_to_v[pcgn]->index()
				   << " [" << u->index() << "]." << pcgn->index() << "->[" << v->index() << "]";
				GA->label(pcg_to_v[pcgn]) = ss.str();
			}

			auto& s = log.lout(Logger::Level::Medium)
					<< "PC-Tree " << (fixed_rotation ? "Q" : "P") << " node #" << pcgn->index()
					<< " °" << pcgn->degree() << ": ";
			if (fixed_rotation) {
				int p = partitions.makeQVertex(pcg_to_u[pcgn]);
				partitions.makeQVertex(pcg_to_v[pcgn], p);
				if (pcgn->degree() > 3) {
					make_wheels.pushBack(pcg_to_v[pcgn]);
					if (!v_was_cut) { // the intersection might have generated Q-nodes from v, which are not ensured at u
						make_wheels.pushBack(pcg_to_u[pcgn]);
					}
				}
				s << "Q-nodes " << fmtPQNode(pcg_to_u[pcgn]) << " and " << fmtPQNode(pcg_to_v[pcgn])
				  << " are in same partition cell " << p << " and "
				  << (v_was_cut ? "the v node" : "both nodes")
				  << (pcgn->degree() > 3 ? " will" : " are too small to") << " be made a wheel."
				  << std::endl;
			} else {
				OGDF_ASSERT(pctn->getNodeType() == PCNodeType::PNode);
				matchings.matchNodes(pcg_to_u[pcgn], pcg_to_v[pcgn]);
				s << "P-node " << fmtPQNode(pcg_to_u[pcgn]) << " is matched with "
				  << fmtPQNode(pcg_to_v[pcgn]) << "." << std::endl;
			}
			log.lout(Logger::Level::Minor)
					<< printBijection(getPipeBijection(pcg_to_u[pcgn], pcg_to_v[pcgn])) << std::endl;
		}
	}

	OGDF_ASSERT(inner_pcg_node != nullptr);
	if (u_bc->degree() == 0 && components.bcRepr(u_bc) == u) {
		components.makeRepr(u_bc, pcg_to_u[inner_pcg_node]);
	}

	if (v_was_cut) {
		log.lout() << "Relabeling exploded star with " << v_rays.size() << " rays: " << v_rays
				   << std::endl;
		components.relabelExplodedStar(v_bc, nullptr, v_rays);
	} else {
		if (v_bc->degree() == 0 && components.bcRepr(v_bc) == v) {
			components.makeRepr(v_bc, pcg_to_v[inner_pcg_node]);
		}
		components.bc_size[v_bc]--; // room for improvement: move to components callback function
	}
	components.bc_size[u_bc]--;

	OGDF_ASSERT(u->degree() == 0);
	OGDF_ASSERT(v->degree() == 0);
#ifdef OGDF_DEBUG
	deletedNodes.insert(u);
	deletedNodes.insert(v);
#endif
	pushUndoOperationAndCheck(undo_op);

	for (node w : make_wheels) {
		makeWheel(w, true);
	}

	if (v_was_cut) {
		log.lout() << "Undoing encapsulation of v by contracting the " << v_rays.size() << " rays."
				   << std::endl;
		for (node ray : v_rays) {
			// this contraction is always allowed, no matter what allow_contract_bb_pipe says
			bool ac = true;
			std::swap(ac, allow_contract_bb_pipe);
			Result result = contract(ray);
			std::swap(ac, allow_contract_bb_pipe);
			OGDF_ASSERT(result == SyncPlan::Result::SUCCESS);
		}
	}

#ifdef SYNCPLAN_OPSTATS
	printOPStatsEnd(true, dur_ns(tpc::now() - start) - v_pc_dur);
#endif
	// SYNCPLAN_PROFILE_STOP("propagatePQ")
	return SyncPlan::Result::SUCCESS;
}

}
