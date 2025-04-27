/** \file
 * \brief Implementation of the SyncPlan::embed operation
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
#include <ogdf/cluster/sync_plan/PMatching.h>
#include <ogdf/cluster/sync_plan/SyncPlan.h>
#include <ogdf/cluster/sync_plan/SyncPlanConsistency.h>
#include <ogdf/cluster/sync_plan/utils/Bijection.h>
#include <ogdf/cluster/sync_plan/utils/Logging.h>

#include <sstream>
#include <string>

using namespace ogdf::sync_plan::internal;

namespace ogdf::sync_plan {

class UpdateGraphReg : public GraphObserver {
	NodeArray<node>* node_reg;
	EdgeArray<edge>* edge_reg;

public:
	UpdateGraphReg(const Graph* g, NodeArray<node>* nodeReg, EdgeArray<edge>* edgeReg)
		: GraphObserver(), node_reg(nodeReg), edge_reg(edgeReg) {
		reregister(g);
		nodeReg->init(*g, nullptr);
		edgeReg->init(*g, nullptr);
		for (node n : g->nodes) {
			(*node_reg)[n] = n;
		}
		for (edge e : g->edges) {
			(*edge_reg)[e] = e;
		}
	}

	~UpdateGraphReg() override {
		node_reg->init();
		edge_reg->init();
	}

	void nodeDeleted(node v) override { (*node_reg)[v] = nullptr; }

	void nodeAdded(node v) override { (*node_reg)[v] = v; }

	void edgeDeleted(edge e) override { (*edge_reg)[e] = nullptr; }

	void edgeAdded(edge e) override { (*edge_reg)[e] = e; }

	void cleared() override { }
};

node SyncPlan::nodeFromIndex(int idx) const {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
	return node_reg[idx];
#pragma GCC diagnostic pop
}

edge SyncPlan::edgeFromIndex(int idx) const {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
	return edge_reg[idx];
#pragma GCC diagnostic pop
}

void SyncPlan::thawPipeBijection(node u, node v, const FrozenPipeBij& in, PipeBij& out) const {
	for (const FrozenPipeBijPair& pair : in) {
		out.emplaceBack(edgeFromIndex(pair.first)->getAdj(u), edgeFromIndex(pair.second)->getAdj(v));
	}
}

bool SyncPlan::verifyPipeBijection(node u, node v, const FrozenPipeBij& bij) const {
	PipeBij thawed_bij;
	thawPipeBijection(u, v, bij, thawed_bij);

	PipeBij new_bij;
	matchings.getIncidentEdgeBijection(u, new_bij);

	const PipeBijCmp& cmp = PipeBijCmp();
	thawed_bij.quicksort(cmp);
	new_bij.quicksort(cmp);
	bool bijection_broke = thawed_bij != new_bij;
	if (bijection_broke) {
		log.lout(Logger::Level::Alarm) << "old_bij: " << printBijection(thawed_bij) << std::endl;
		log.lout(Logger::Level::Alarm) << "new_bij: " << printBijection(new_bij) << std::endl;
	}
	OGDF_ASSERT(!bijection_broke);
	return !bijection_broke;
}

void SyncPlan::embed() {
	OGDF_ASSERT(G->representsCombEmbedding());

	// SYNCPLAN_PROFILE_START("embed")
	int undo_cnt = undo_stack.size();
	log.lout(Logger::Level::High) << undo_cnt << " Operations to undo" << std::endl;
	UpdateGraphReg updater(G, &node_reg, &edge_reg);
	for (edge e : deletedEdges) {
		edge_reg[e] = e;
	}
	OGDF_ASSERT(matchings.isReduced());
	matchings.setPipeQueue(nullptr);
	while (!undo_stack.empty()) {
		// SYNCPLAN_PROFILE_START("embed-step")
		UndoOperation* op = undo_stack.popBackRet();
		log.lout(Logger::Level::High) << (undo_cnt - undo_stack.size())
#ifdef OGDF_DEBUG
									  << " (" << op->consistency_nr << ")"
#endif
									  << ": " << *op << std::endl;

#ifdef SYNCPLAN_OPSTATS
		std::chrono::time_point<std::chrono::high_resolution_clock> start = tpc::now();
#endif
		op->undo(*this);
#ifdef SYNCPLAN_OPSTATS
		stats_out << (stats_first_in_array ? "" : ",") << "{\"op\":\"embed-" << op->name() << "\""
				  << ",\"op_time_ns\":" << dur_ns(tpc::now() - start) << "}";
		stats_first_in_array = false;
#endif

#ifdef OGDF_DEBUG
		if (consistency.doWriteOut) {
			std::stringstream ss;
			ss << "undoOp" << op->consistency_nr;
			consistency.writeOut(ss.str(), false, false);
		}
#endif
		delete op;
		OGDF_ASSERT(G->representsCombEmbedding());
		// SYNCPLAN_PROFILE_STOP("embed-step")
	}
	// SYNCPLAN_PROFILE_STOP("embed")
}

}
