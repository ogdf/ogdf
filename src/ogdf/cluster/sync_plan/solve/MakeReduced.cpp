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
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/Logger.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/pctree/NodePCRotation.h>
#include <ogdf/basic/pctree/PCEnum.h>
#include <ogdf/cluster/sync_plan/PMatching.h>
#include <ogdf/cluster/sync_plan/PQPlanarity.h>
#include <ogdf/cluster/sync_plan/PQPlanarityComponents.h>
#include <ogdf/cluster/sync_plan/PQPlanarityConsistency.h>
#include <ogdf/cluster/sync_plan/PipeOrder.h>
#include <ogdf/cluster/sync_plan/utils/Logging.h>

#include <cmath>
#include <cstddef>
#include <memory>
#include <ostream>

PQPlanarity::Result PQPlanarity::checkPCTree(node u) {
	try {
#ifdef SYNCPLAN_OPSTATS
		tp pc_start = tpc::now();
#endif
		// SYNCPLAN_PROFILE_START("checkPCTree")
		BiconnectedIsolation iso(components, components.biconnectedComponent(u));
		NodePCRotation pc(*G, u, true);
		iso.restore(); // TODO Make function of components to get EmbeddingTree interface
		log.lout() << "PC-Tree with " << pc.getPNodeCount() << " P-nodes and " << pc.getCNodeCount()
				   << " C-nodes" << std::endl;
		if (log.is_lout(ogdf::Logger::Level::Minor)) {
			log.lout() << "PC-Tree: " << pc << " (" << (pc.isTrivial() ? "" : "non-") << "trivial)*"
					   << pc.possibleOrders<size_t>() << "/" << round(tgamma(u->degree()))
					   << std::endl;
		}
		// SYNCPLAN_PROFILE_STOP("checkPCTree")
#ifdef SYNCPLAN_OPSTATS
		stats_pc_time = dur_ns(tpc::now() - pc_start);
#endif
		Result result;
		if (!pc.isTrivial()) {
			result = propagatePQ(u, &pc);
		} else {
			result = simplify(u, &pc);
		}
		return result;
	} catch (pc_tree::GraphNotPlanarException&) {
		log.lout(Logger::Level::Alarm)
				<< "Instance became non-planar during reduction, so creating a PC-Tree failed!"
				<< std::endl;
		return PQPlanarity::Result::INVALID_INSTANCE;
	}
}

bool PQPlanarity::makeReduced(int check_planarity_every) {
	// SYNCPLAN_PROFILE_START("makeReduced")
	if (!indices_saved) {
		undo_stack.pushBack(new ResetIndices(*this));
	}
	if (!matchings.getPipeQueue()) {
		log.lout(Logger::Level::Minor) << "Using default PipeQueueByDegree" << std::endl;
		matchings.setPipeQueue(std::make_unique<PipeQueueByDegree>());
	}
	pushUndoOperationAndCheck(new VerifyPipeBijections(*this));
	int steps = 0;
	while (!matchings.isReduced()) {
		if (check_planarity_every > 0 && steps % check_planarity_every == 0 && !isPlanar(*G)) {
			log.lout(Logger::Level::Alarm)
					<< "Instance became non-planar during reduction!" << std::endl;
			return false;
		}
		// SYNCPLAN_PROFILE_START("makeReduced-step")
		steps++;
		const Pipe& pipe = matchings.getTopPipe();

		int prio_pipe_cnt = matchings.getPriorityPipeCount(), top_prio = pipe.pipe_priority;
		if (log.is_lout(Logger::Level::High)) {
			log.lout(Logger::Level::High)
					<< "Step " << steps << ": " << matchings.getPipeCount() << " pipes left"
					<< (log.is_lout(Logger::Level::Medium) ? ":" : "") << std::endl;
			log.lout(Logger::Level::Medium) << matchings.getPipes() << std::endl;
			log.lout(Logger::Level::High)
					<< "Processing pipe matching " << fmtPQNode(pipe.node1) << std::endl
					<< "                    with " << fmtPQNode(pipe.node2) << std::endl
					<< "           of components " << components.fmtBCOf(pipe.node1) << std::endl
					<< "                     and " << components.fmtBCOf(pipe.node2) << "."
					<< std::endl;
			if (pipe.pipe_priority >= 0) {
				log.lout(Logger::Level::High)
						<< "Pipe has priority " << top_prio << ", priority pipe count is "
						<< prio_pipe_cnt << std::endl;
			} else {
				OGDF_ASSERT(prio_pipe_cnt == 0);
			}
		}
		Logger::Indent _(&log);

		if (pipe.degree() <= 3) {
			Result result = convertSmall(pipe.node1);
			OGDF_ASSERT(result == PQPlanarity::Result::SUCCESS);
			// SYNCPLAN_PROFILE_STOP("makeReduced-step")
			continue;
		}

		if (canContract(&pipe)) {
#ifdef SYNCPLAN_OPSTATS
			tp contract_start = tpc::now();
			printOPStatsStart(matchings.getPipe(pipe.node1),
					components.isCutVertex(pipe.node1) ? Operation::ENCAPSULATE_CONTRACT
													   : Operation::CONTRACT_BICON);
#endif
			Result contract_result = contract(pipe.node1);
			OGDF_ASSERT(contract_result == PQPlanarity::Result::SUCCESS);
#ifdef SYNCPLAN_OPSTATS
			printOPStatsEnd(true, dur_ns(tpc::now() - contract_start));
#endif
			// SYNCPLAN_PROFILE_STOP("makeReduced-step")
			continue;
		}

		if (batch_spqr) {
			Result batch_result = batchSPQR();
			if (batch_result == PQPlanarity::Result::INVALID_INSTANCE) {
				// SYNCPLAN_PROFILE_STOP("makeReduced-step")
				return false;
			} else if (batch_result == PQPlanarity::Result::SUCCESS) {
				// SYNCPLAN_PROFILE_STOP("makeReduced-step")
				continue;
			}
		}

		node block_vertex;
		if (!components.isCutVertex(pipe.node1)) {
			block_vertex = pipe.node1;
		} else {
			OGDF_ASSERT(!components.isCutVertex(pipe.node2));
			block_vertex = pipe.node2;
		}

		Result tree_result = checkPCTree(block_vertex);
		// SYNCPLAN_PROFILE_STOP("makeReduced-step")
		if (tree_result == PQPlanarity::Result::NOT_APPLICABLE) {
			OGDF_ASSERT(matchings.getTopPipe().pipe_priority > top_prio);
			continue;
		} else if (tree_result == PQPlanarity::Result::SUCCESS) {
			continue;
		} else if (tree_result == PQPlanarity::Result::INVALID_INSTANCE) {
			return false;
		}
	}
	OGDF_ASSERT(consistency.consistencyCheck());
	// SYNCPLAN_PROFILE_STOP("makeReduced")
	return true;
}
