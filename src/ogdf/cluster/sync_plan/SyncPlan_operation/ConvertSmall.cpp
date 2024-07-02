/** \file
 * \brief Implementation of the SyncPlan::converSmall operation and its UndoOperation.
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
#include <ogdf/basic/Logger.h>
#include <ogdf/basic/basic.h>
#include <ogdf/cluster/sync_plan/PMatching.h>
#include <ogdf/cluster/sync_plan/QPartitioning.h>
#include <ogdf/cluster/sync_plan/SyncPlan.h>
#include <ogdf/cluster/sync_plan/SyncPlanComponents.h>
#include <ogdf/cluster/sync_plan/basic/GraphUtils.h>
#include <ogdf/cluster/sync_plan/utils/Bijection.h>
#include <ogdf/cluster/sync_plan/utils/Logging.h>

#include <ostream>

using namespace ogdf::sync_plan::internal;

namespace ogdf::sync_plan {
using internal::operator<<;

class UndoConvertSmall : public SyncPlan::UndoOperation {
public:
	int small_idx, twin_idx;
	int small_first_adj_idx, twin_last_adj_idx;
#ifdef OGDF_DEBUG
	FrozenPipeBij bij;
#endif

	UndoConvertSmall(node small, node twin) : small_idx(small->index()), twin_idx(twin->index()) {
		if (small->degree() > 0) {
			small_first_adj_idx = small->adjEntries.head()->theEdge()->index();
			twin_last_adj_idx = twin->adjEntries.tail()->theEdge()->index();
		} else {
			small_first_adj_idx = twin_last_adj_idx = -1;
		}
#ifdef OGDF_DEBUG
		getFrozenPipeBijection(small, twin, bij);
#endif
	}

	void undo(SyncPlan& pq) override {
		// SYNCPLAN_PROFILE_START("undo-convertSmall")
		node small = pq.nodeFromIndex(small_idx);
		node twin = pq.nodeFromIndex(twin_idx);
		if (small_first_adj_idx >= 0) {
			moveAdjToFront(*pq.G, pq.edgeFromIndex(small_first_adj_idx)->getAdj(small));
			moveAdjToBack(*pq.G, pq.edgeFromIndex(twin_last_adj_idx)->getAdj(twin));
		}
		if (pq.partitions.isQVertex(small)) {
			OGDF_ASSERT(small->degree() > 2);
			pq.partitions.releaseQVertex(small);
			pq.partitions.releaseQVertex(twin);
		} else {
			OGDF_ASSERT(small->degree() <= 2);
		}
		pq.matchings.matchNodes(small, twin);

#ifdef OGDF_DEBUG
		pq.verifyPipeBijection(small, twin, bij);
#endif
		// SYNCPLAN_PROFILE_STOP("undo-convertSmall")
	}

	std::ostream& print(std::ostream& os) const override {
		return os << "UndoConvertSmall(small=" << small_idx << ", twin=" << twin_idx
#ifdef OGDF_DEBUG
				  << ", bij=" << printFrozenBijection(bij)
#endif
				  << ")";
	}
};

SyncPlan::Result SyncPlan::convertSmall(node small) {
	if (small->degree() > 4 || partitions.isQVertex(small)) {
		return SyncPlan::Result::NOT_APPLICABLE;
	} else if (matchings.isMatchedPVertex(small)) {
		// SYNCPLAN_PROFILE_START("convertSmall")
		log.lout(Logger::Level::High) << "CONVERT SMALL degree " << small->degree() << std::endl;
		log.lout(Logger::Level::Minor) << matchings.printBijection(small) << std::endl;
		node twin = matchings.removeMatching(small);
		if (small->degree() > 2) {
			// this forces the edge bijection to stay valid
			int partition = partitions.makeQVertex(small);
			partitions.makeQVertex(twin, partition);
		}
		pushUndoOperation(new UndoConvertSmall(small, twin));
		if (small->degree() > 2) {
			if (components.isCutVertex(small)) {
				makeWheel(small);
			}
			if (components.isCutVertex(twin)) {
				makeWheel(twin);
			}
		}
		formatNode(twin);
		// SYNCPLAN_PROFILE_STOP("convertSmall")
	}
	formatNode(small);
	return SyncPlan::Result::SUCCESS;
}

}
