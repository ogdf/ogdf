/** \file
 * \brief Implementation of the SyncPlan::batchSPQR operation and its UndoOperation.
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
#include <ogdf/basic/List.h>
#include <ogdf/basic/Logger.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/pctree/NodePCRotation.h>
#include <ogdf/basic/pctree/PCEnum.h>
#include <ogdf/basic/pctree/util/FilteringBFS.h>
#include <ogdf/cluster/sync_plan/PMatching.h>
#include <ogdf/cluster/sync_plan/SyncPlan.h>
#include <ogdf/cluster/sync_plan/SyncPlanComponents.h>
#include <ogdf/cluster/sync_plan/basic/OverlappingGraphCopies.h>
#include <ogdf/cluster/sync_plan/utils/Logging.h>
#include <ogdf/cluster/sync_plan/utils/NodeTricRotation.h>

#include <cmath>
#include <cstddef>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

using namespace ogdf::pc_tree;
using namespace ogdf::sync_plan::internal;

namespace ogdf::sync_plan {

class EmbeddingTrees {
public:
	SyncPlan& pq;
	OverlappingGraphCopies OGC_base;
	NodeArray<NodeSSPQRRotation*> embtrees;
	std::vector<NodeSSPQRRotation*> e_todelete;
	NodeArray<SimpleSPQRTree*> spqrtrees;
	std::vector<SimpleSPQRTree*> s_todelete;

	explicit EmbeddingTrees(SyncPlan& _pq)
		: pq(_pq)
		, OGC_base(*_pq.G)
		, embtrees(*_pq.G, nullptr)
		, spqrtrees(_pq.getComponents().bcTree(), nullptr) {
		e_todelete.reserve(_pq.matchings.getPipeCount());
		s_todelete.reserve(_pq.components.bcTree().numberOfNodes());
	}

	~EmbeddingTrees() {
		for (auto ptr : e_todelete) {
			delete ptr;
		}
		for (auto ptr : s_todelete) {
			ptr->GC.breakLinkForMasterDeconstruction();
			delete ptr;
		}
	}

	NodePCRotation* operator[](node n) { return embtrees[n]; }

	NodePCRotation* makeTree(node n) {
		if (embtrees[n] != nullptr) {
			return embtrees[n];
		}
		using internal::operator<<;

		SyncPlanComponents& components = pq.components;
		node bc = components.biconnectedComponent(n);
		if (spqrtrees[bc] == nullptr) {
			spqrtrees[bc] = new SimpleSPQRTree(OGC_base);
			s_todelete.push_back(spqrtrees[bc]);
		}
		SimpleSPQRTree& spqr = *spqrtrees[bc];
		if (spqr.GC.empty()) {
#ifdef SYNCPLAN_OPSTATS
			tp start = tpc::now();
			if (!pq.stats_first_in_array) {
				pq.stats_out << ",";
			} else {
				pq.stats_first_in_array = false;
			}
			pq.stats_out << "{\"op\":\"MAKE_SPQR\""
						 << ",\"bicon_size\":" << components.bcSize(bc)
						 << ",\"bc_id\":" << bc->index()
						 << ",\"cc_id\":" << components.bcConnectedId(bc) << ",";
#endif
			// SYNCPLAN_PROFILE_START("batchSPQR-makeSPQR")
			pq.log.lout()
					<< "Creating SPQR information for " << components.fmtBCNode(bc) << std::endl;
			FilteringBFS bfs = components.nodesInBiconnectedComponent(bc);
			for (node u : bfs) {
				spqr.GC.newNode(u);
				for (adjEntry adj : u->adjEntries) {
					if (bfs.hasVisited(adj->twinNode())) {
						spqr.GC.newEdge(adj->theEdge());
					}
				}
			}
			OGDF_ASSERT(!bfs.valid());
			pq.log.lout() << "Copied " << spqr.GC.numberOfNodes() << " nodes and "
						  << spqr.GC.numberOfEdges() << " edges" << std::endl;
			OGDF_ASSERT(spqr.GC.numberOfNodes() == components.bcSize(bc));
#ifdef SYNCPLAN_OPSTATS
			pq.stats_out << "\"nodes\":" << spqr.GC.numberOfNodes()
						 << ",\"edges\":" << spqr.GC.numberOfEdges() << ",";
#endif
			spqr.init();
			// SYNCPLAN_PROFILE_STOP("batchSPQR-makeSPQR")
#ifdef SYNCPLAN_OPSTATS
			pq.stats_out << "\"op_time_ns\":" << dur_ns(tpc::now() - start) << "}";
#endif
		}
		if (!spqr.planar) {
			return nullptr;
		}

		// SYNCPLAN_PROFILE_START("batchSPQR-makeTree")
#ifdef SYNCPLAN_OPSTATS
		tp start = tpc::now();
#endif
		pq.log.lout() << "Computing NodeSPQRRotation for vertex " << pq.fmtPQNode(n) << " in block "
					  << components.fmtBCOf(n) << std::endl;
		NodeSSPQRRotation* pc = new NodeSSPQRRotation(spqr, n);
		pq.log.lout() << "PC-Tree with " << pc->getPNodeCount() << " P-nodes and "
					  << pc->getCNodeCount() << " C-nodes" << std::endl;
		if (pq.log.is_lout(ogdf::Logger::Level::Minor)) {
			pq.log.lout() << "PC-Tree: " << *pc << " (" << (pc->isTrivial() ? "" : "non-")
						  << "trivial)*" << pc->possibleOrders<size_t>() << "/"
						  << round(tgamma(pc->getNode()->degree())) << std::endl;
			if (pc->isTrivial()) {
				pq.log.lout() << "Partner vertex is " << pq.fmtPQNode(pc->getTrivialPartnerPole())
							  << " in block " << components.fmtBCOf(pc->getTrivialPartnerPole())
							  << std::endl;
			}
		}
		pc->mapPartnerEdges();
		e_todelete.emplace_back(pc);
		embtrees[n] = pc;
#ifdef SYNCPLAN_OPSTATS
		pq.stats_pc_time += dur_ns(tpc::now() - start);
#endif
		// SYNCPLAN_PROFILE_STOP("batchSPQR-makeTree")
#ifdef OGDF_HEAVY_DEBUG
		{
			BiconnectedIsolation iso(components, components.biconnectedComponent(n));
			NodePCRotation pc2(*pq.G, n);
			if (!pc->isEqual(pc2)) {
				pq.log.lout(Logger::Level::Alarm)
						<< "NodeSPQRRotation: " << pc->uniqueID(pc->uidPrinter(), pc->uidComparer())
						<< std::endl;
				pq.log.lout(Logger::Level::Alarm)
						<< "NodePCRotation:   " << pc2.uniqueID(pc2.uidPrinter(), pc2.uidComparer())
						<< std::endl;
				OGDF_ASSERT(false);
			}
		}
#endif
		return pc;
	}
};

struct SimplePipe {
	node block_vertex, other_vertex;
	bool both_block = false;

	SimplePipe(const Pipe& p, const SyncPlanComponents& components)
		: block_vertex(p.node1), other_vertex(p.node2) {
		if (!components.isCutVertex(block_vertex)) {
			both_block = !components.isCutVertex(other_vertex);
		} else {
			using std::swap;
			swap(block_vertex, other_vertex);
		}
		OGDF_ASSERT(!components.isCutVertex(block_vertex));
	}
};

#ifdef SYNCPLAN_OPSTATS
#	define RETURN_INVALID                                  \
		printOPStatsEnd(false, dur_ns(tpc::now() - start)); \
		return SyncPlan::Result::INVALID_INSTANCE;
// SYNCPLAN_PROFILE_STOP("batchSPQR")
#else
#	define RETURN_INVALID return SyncPlan::Result::INVALID_INSTANCE;
// SYNCPLAN_PROFILE_STOP("batchSPQR")
#endif

SyncPlan::Result SyncPlan::batchSPQR() {
	// SYNCPLAN_PROFILE_START("batchSPQR")
	log.lout(Logger::Level::High) << "BATCH SPQR" << std::endl;
#ifdef SYNCPLAN_OPSTATS
	tp start = tpc::now();
	if (!stats_first_in_array) {
		stats_out << ",";
	}
	stats_out << "{\"op\":\"BATCH_SPQR\",\"rem_pipes\":" << matchings.getPipeCount() << ",\"ops\":[";
	stats_first_in_array = true;
#endif
	Logger::Indent _(log);
	SimpleSPQRTree::log = log; // copy settings
	SimpleSPQRTree::log.indent();
	EmbeddingTrees embtrees(*this);
	bool changed = false, simplified = true;

	std::vector<SimplePipe> pipes;
	pipes.reserve(matchings.getPipeCount());
	while (simplified) {
		simplified = false;
		for (const Pipe& p : matchings.getPipes()) {
			if (p.degree() <= 3 || canContract(&p)) {
				continue;
			}
			pipes.emplace_back(p, components);
		}
		for (const SimplePipe& p : pipes) {
			if (!matchings.isMatchedPVertex(p.block_vertex)
					|| !matchings.isMatchedPVertex(p.other_vertex)) {
				continue; // pipe was simplified transitively in the meantime...
			}

			NodePCRotation* pc = embtrees.makeTree(p.block_vertex);
			if (pc == nullptr) {
				RETURN_INVALID
			}
			if (!pc->isTrivial() && p.both_block && intersect_trees) {
				if (embtrees.makeTree(p.other_vertex) == nullptr) {
					RETURN_INVALID
				}
			}

			// SYNCPLAN_PROFILE_START("batchSPQR-simplify")
			Result r = simplify(p.block_vertex, pc);
			// SYNCPLAN_PROFILE_STOP("batchSPQR-simplify")
			if (r == SyncPlan::Result::INVALID_INSTANCE) {
				RETURN_INVALID
			} else if (r == SyncPlan::Result::SUCCESS) {
				changed = simplified = true;
			} else {
				OGDF_ASSERT(r == SyncPlan::Result::NOT_APPLICABLE);

				if (p.both_block) {
					// SYNCPLAN_PROFILE_START("batchSPQR-simplify")
					NodePCRotation* opc = embtrees.makeTree(p.other_vertex);
					r = SyncPlan::Result::INVALID_INSTANCE;
					if (opc != nullptr) {
						r = simplify(p.other_vertex, opc);
					}
					// SYNCPLAN_PROFILE_STOP("batchSPQR-simplify")
					if (r == SyncPlan::Result::INVALID_INSTANCE) {
						RETURN_INVALID
					} else if (r == SyncPlan::Result::SUCCESS) {
						changed = simplified = true;
					}
				}
			}
		}
		pipes.clear();
	}

	for (const Pipe& p : matchings.getPipes()) {
		if (p.degree() <= 3 || canContract(&p)) {
			continue;
		}
		pipes.emplace_back(p, components);
	}
	for (const SimplePipe& p : pipes) {
		NodePCRotation* pc = embtrees[p.block_vertex];
		NodePCRotation* pc_v = embtrees[p.other_vertex];
		OGDF_ASSERT(pc != nullptr);
		if (pc->isTrivial() && (!pc_v || pc_v->isTrivial())) {
			continue;
		}

		Result r;
		// SYNCPLAN_PROFILE_START("batchSPQR-propagatePQ")
		if (pc_v) {
			if (pc->isTrivial()
					|| (!intersect_trees
							&& pc_v->possibleOrders<size_t>() < pc->possibleOrders<size_t>())) {
				r = propagatePQ(p.other_vertex, pc_v, pc);
			} else {
				r = propagatePQ(p.block_vertex, pc, pc_v);
			}
		} else {
			r = propagatePQ(p.block_vertex, pc);
		}
		// SYNCPLAN_PROFILE_STOP("batchSPQR-propagatePQ")
		if (r == SyncPlan::Result::INVALID_INSTANCE) {
			RETURN_INVALID
		} else if (r == SyncPlan::Result::SUCCESS) {
			changed = true;
		}
	}

	// OGDF_ASSERT(matchings.getPipeCount() < pipe_cnt); // number of pipes might actually have increased,
	// // but their total degree decreased (a value which we unfortunately do not maintain)
	// SYNCPLAN_PROFILE_STOP("batchSPQR")
#ifdef SYNCPLAN_OPSTATS
	stats_out << "],";
	stats_first_in_array = false;
	printOPStatsEnd(changed, dur_ns(tpc::now() - start));
#endif
	return changed ? SyncPlan::Result::SUCCESS : SyncPlan::Result::NOT_APPLICABLE;
}

}
