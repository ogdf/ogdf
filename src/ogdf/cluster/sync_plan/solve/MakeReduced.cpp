#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/cluster/sync_plan/PQPlanarity.h>
#include <ogdf/cluster/sync_plan/PipeOrder.h>
#include <ogdf/cluster/sync_plan/utils/Logging.h>

#include <NodePCRotation.h>

PQPlanarity::Result PQPlanarity::checkPCTree(node u) {
	try {
#ifdef PQ_OPSTATS
		tp pc_start = tpc::now();
#endif
		PQ_PROFILE_START("checkPCTree")
		BiconnectedIsolation iso(components, components.biconnectedComponent(u));
		NodePCRotation pc(*G, u, true);
		iso.restore(); // TODO Make function of components to get EmbeddingTree interface
		log.lout() << "PC-Tree with " << pc.getPNodeCount() << " P-nodes and " << pc.getCNodeCount()
				   << " C-nodes" << endl;
		if (log.is_lout(ogdf::Logger::Level::Minor)) {
			log.lout() << "PC-Tree: " << pc << " (" << (pc.isTrivial() ? "" : "non-") << "trivial)*"
					   << pc.possibleOrders() << "/" << round(tgamma(u->degree())) << endl;
		}
		PQ_PROFILE_STOP("checkPCTree")
#ifdef PQ_OPSTATS
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
		return INVALID_INSTANCE;
	}
}

bool PQPlanarity::makeReduced(int check_planarity_every) {
	PQ_PROFILE_START("makeReduced")
	if (!indices_saved) {
		undo_stack.pushBack(new ResetIndices(*this));
	}
	if (!matchings.getPipeQueue()) {
		log.lout(Logger::Level::Minor) << "Using default PipeQueueByDegree" << std::endl;
		matchings.setPipeQueue(make_unique<PipeQueueByDegree>());
	}
	pushUndoOperationAndCheck(new VerifyPipeBijections(*this));
	int steps = 0;
	while (!matchings.isReduced()) {
		if (check_planarity_every > 0 && steps % check_planarity_every == 0 && !isPlanar(*G)) {
			log.lout(Logger::Level::Alarm)
					<< "Instance became non-planar during reduction!" << std::endl;
			return false;
		}
		PQ_PROFILE_START("makeReduced-step")
		steps++;
		const Pipe& pipe = matchings.getTopPipe();

		int prio_pipe_cnt = matchings.getPriorityPipeCount(), top_prio = pipe.pipe_priority;
		if (log.is_lout(Logger::Level::High)) {
			log.lout(Logger::Level::High)
					<< "Step " << steps << ": " << matchings.getPipeCount() << " pipes left"
					<< (log.is_lout(Logger::Level::Medium) ? ":" : "") << endl;
			log.lout(Logger::Level::Medium) << matchings.getPipes() << endl;
			log.lout(Logger::Level::High)
					<< "Processing pipe matching " << fmtPQNode(pipe.node1) << endl
					<< "                    with " << fmtPQNode(pipe.node2) << endl
					<< "           of components " << components.fmtBCOf(pipe.node1) << endl
					<< "                     and " << components.fmtBCOf(pipe.node2) << "." << endl;
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
			OGDF_ASSERT(result == SUCCESS);
			PQ_PROFILE_STOP("makeReduced-step")
			continue;
		}

		if (canContract(&pipe)) {
#ifdef PQ_OPSTATS
			tp contract_start = tpc::now();
			printOPStatsStart(matchings.getPipe(pipe.node1),
					components.isCutVertex(pipe.node1) ? Operation::ENCAPSULATE_CONTRACT
													   : Operation::CONTRACT_BICON);
#endif
			Result contract_result = contract(pipe.node1);
			OGDF_ASSERT(contract_result == SUCCESS);
#ifdef PQ_OPSTATS
			printOPStatsEnd(true, dur_ns(tpc::now() - contract_start));
#endif
			PQ_PROFILE_STOP("makeReduced-step")
			continue;
		}

		if (batch_spqr) {
			Result batch_result = batchSPQR();
			if (batch_result == INVALID_INSTANCE) {
				PQ_PROFILE_STOP("makeReduced-step")
				return false;
			} else if (batch_result == SUCCESS) {
				PQ_PROFILE_STOP("makeReduced-step")
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
		PQ_PROFILE_STOP("makeReduced-step")
		if (tree_result == NOT_APPLICABLE) {
			OGDF_ASSERT(matchings.getTopPipe().pipe_priority > top_prio);
			continue;
		} else if (tree_result == SUCCESS) {
			continue;
		} else if (tree_result == INVALID_INSTANCE) {
			return false;
		}
	}
	OGDF_ASSERT(consistency.consistencyCheck());
	PQ_PROFILE_STOP("makeReduced")
	return true;
}
