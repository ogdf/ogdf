/** \file
 * \brief Implementation of the SyncPlan::makeWheel and SyncPlan::contractWheel operations.
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
#include <ogdf/basic/List.h>
#include <ogdf/basic/Logger.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/basic.h>
#include <ogdf/cluster/sync_plan/SyncPlan.h>
#include <ogdf/cluster/sync_plan/SyncPlanComponents.h>
#include <ogdf/cluster/sync_plan/utils/Logging.h>

#include <sstream>
#include <string>

using namespace ogdf::sync_plan::internal;

namespace ogdf::sync_plan {
using internal::operator<<;

class UndoMakeWheel : public SyncPlan::UndoOperation {
public:
	int centre_idx;

	UndoMakeWheel(node centre) : centre_idx(centre->index()) { }

	void undo(SyncPlan& pq) override { pq.contractWheel(pq.nodeFromIndex(centre_idx)); }

	std::ostream& print(std::ostream& os) const override {
		return os << "UndoMakeWheel(" << centre_idx << ")";
	}
};

void SyncPlan::makeWheel(node centre, bool update_cuts) {
	// SYNCPLAN_PROFILE_START("makeWheel")
	bool is_cut = update_cuts && components.isCutVertex(centre);
	log.lout(Logger::Level::High) << "MAKE WHEEL centre " << fmtPQNode(centre) << " (update cut "
								  << update_cuts << is_cut << ")" << std::endl;
	Logger::Indent _(&log);
	log.lout(Logger::Level::Medium) << printIncidentEdges(centre->adjEntries) << std::endl;
	NodeArray<SList<adjEntry>> block_neigh;
	if (is_cut) {
		block_neigh.init(components.bcTree());
	}

	node prev_ray = nullptr, first_ray = nullptr;
	node centre_bc = components.biconnectedComponent(centre);
	List<adjEntry> adjEntries;
	centre->allAdjEntries(adjEntries);
	for (adjEntry adj : adjEntries) {
		node neigh = adj->twinNode();
		edge the_edge = adj->theEdge();
		bool reverse = adj->isSource();
		if (reverse) {
			G->reverseEdge(the_edge);
		}
		edge new_edge = G->split(the_edge);
		OGDF_ASSERT(neigh == the_edge->source());
		OGDF_ASSERT(the_edge->target() == new_edge->source());
		OGDF_ASSERT(new_edge->target() == centre);
		node spike = new_edge->source();
		is_wheel[spike] = true;
		if (reverse) {
			G->reverseEdge(the_edge);
			G->reverseEdge(new_edge);
		}

		if (GA != nullptr) {
			std::ostringstream ss;
			ss << "WheelSpike " << spike->index() << " for " << centre->index() << "-"
			   << neigh->index();
			GA->label(spike) = ss.str();
		}
		if (is_cut) {
			node neigh_bc = components.biconnectedComponent(neigh);
			if (components.isCutComponent(neigh_bc)) {
				block_neigh[components.findCommonBiconComp(centre_bc, neigh_bc)].pushBack(adj);
			} else {
				block_neigh[neigh_bc].pushBack(adj);
			}
		}
		components.nodeInserted(spike, centre_bc);

		if (prev_ray != nullptr) {
			G->newEdge(prev_ray, spike);
		}
		prev_ray = spike;
		if (first_ray == nullptr) {
			first_ray = spike;
		}
	}
	if (first_ray != nullptr) {
		G->newEdge(prev_ray, first_ray);
	}
	is_wheel[centre] = true;
	log.lout(Logger::Level::Medium) << printIncidentEdges(centre->adjEntries) << std::endl;

	if (is_cut) {
		components.cutReplacedByWheel(centre,
				block_neigh); // room for improvement: use union-find to merge blocks faster
	}

	pushUndoOperation(new UndoMakeWheel(centre));
	// SYNCPLAN_PROFILE_STOP("makeWheel")
}

void SyncPlan::contractWheel(node centre) {
	// SYNCPLAN_PROFILE_START("contractWheel")
	OGDF_ASSERT(is_wheel[centre]);
	log.lout(Logger::Level::High)
			<< "CONTRACT WHEEL centre " << fmtPQNode(centre, false) << std::endl;
	Logger::Indent _(&log);
	log.lout(Logger::Level::Medium) << printIncidentEdges(centre->adjEntries) << std::endl;
	List<adjEntry> adjEntries;
	centre->allAdjEntries(adjEntries);
	for (adjEntry adj : adjEntries) {
		OGDF_ASSERT(is_wheel[adj->twinNode()]);
		if (!adj->isSource()) {
			G->reverseEdge(adj->theEdge());
		}
		node res = G->contract(adj->theEdge());
		OGDF_ASSERT(res == centre);
	}
	is_wheel[centre] = false;
	log.lout(Logger::Level::Medium) << printIncidentEdges(centre->adjEntries) << std::endl;
	// SYNCPLAN_PROFILE_STOP("contractWheel")
}

}
