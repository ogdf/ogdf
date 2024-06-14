#include <ogdf/basic/extended_graph_alg.h>

#include "PQPlanarity.h"
#include "solve/BlockEmbedding.h"
#include "utils/GraphUtils.h"
#include "utils/TwoSAT.h"

#ifdef PQ_OPSTATS
#	define STEP(op, meta)                                                                         \
		stats_out << (stats_first_in_array ? "" : ",") << "{\"op\":\"solvedReduced-" << op << "\"" \
				  << ",\"op_time_ns\":" << dur_ns(tpc::now() - start) meta << "}";                 \
		stats_first_in_array = false;                                                              \
		start = tpc::now()
#else
#	define STEP(op, meta)
#endif

bool PQPlanarity::solveReduced(bool fail_fast) {
	PQ_PROFILE_START("solveReduced")
#ifdef PQ_OPSTATS
	std::chrono::time_point<std::chrono::high_resolution_clock> start = tpc::now();
#endif
	OGDF_ASSERT(matchings.isReduced());
	// ensure that all Q-node are surrounded by wheels // TODO makeWheel could also be replaced by a Q-vertex-aware embedding tree generator
	int wheels = 0;
	for (int part = 0; part < partitions.partitionCount(); part++) {
		for (node u : partitions.nodesInPartition(part)) {
			if (!is_wheel[u]) {
				// also replace for u->degree() == 3, otherwise multiple following assertions break
				makeWheel(u);
				wheels++;
			}
		}
	}
	STEP("makeWheels", << ",\"wheels\":" << wheels);
	if (fail_fast && !isPlanar(*G)) {
		PQ_PROFILE_STOP("solvedReduced")
		return false;
	}
	log.lout(Logger::Level::High) << "SOLVE REDUCED" << endl;

	// generate SPQR trees for all blocks
	log.lout(Logger::Level::High) << "Generating SPQR Trees for up to "
								  << components.bcTree().numberOfNodes() << " blocks" << endl;
	GnMultiArray Gn_to_subgraph(*G);
	NodeArray<BlockEmbedding> blocks(components.bcTree(), BlockEmbedding(Gn_to_subgraph));
	EdgeArray<edge> Ge_to_subgraph(*G, nullptr);
	EdgeArray<BlockEmbedding*> Ge_to_block(*G, nullptr);
	int block_cnt = 0;
	for (node block : components.bcTree().nodes) {
		if (!components.isCutComponent(block)) {
			blocks[block].init(*G, components, block, Ge_to_subgraph, Ge_to_block);
			block_cnt++;
		}
	}
	STEP("deriveSPQR", << ",\"blocks\":" << block_cnt);

	// generate 2-SAT instance
	log.lout(Logger::Level::High) << "Generating 2-SAT instances for "
								  << partitions.partitionCount() << " partitions" << endl;
	TwoSAT sat;
	PartitionArray<twosat_var> variables(partitions);
	int q_vertices = 0;
	for (int part = 0; part < partitions.partitionCount(); part++) {
		twosat_var part_var = variables[part] = sat.newVariable();
		for (node u : partitions.nodesInPartition(part)) {
			OGDF_ASSERT(partitions.isQVertex(u));
			OGDF_ASSERT(!components.isCutVertex(u));
			q_vertices++;
			if (!blocks[components.biconnectedComponent(u)].addQVertex(u, Ge_to_subgraph, sat,
						part_var)) {
				PQ_PROFILE_STOP("solvedReduced")
				return false;
			}
		}
	}
	STEP("deriveSAT", << ",\"q_vertices\":" << q_vertices);

	log.lout(Logger::Level::High) << "Solving 2-SAT..." << endl;
	if (!sat.solve()) {
		PQ_PROFILE_STOP("solvedReduced")
		return false;
	}
	STEP("solveSAT", );

	// embed rigids according to 2-SAT result
	log.lout(Logger::Level::High) << "Embedding rigids according to 2-SAT solution" << endl;
	int rigids = 0;
	for (node bc : components.bcTree().nodes) {
		BlockEmbedding& block = blocks[bc];
		if (block.spqr) {
			for (node rigid : block.spqr->tree().nodes) {
				if (!fail_fast && !isPlanar(block.spqr->skeleton(rigid).getGraph())) {
					PQ_PROFILE_STOP("solvedReduced")
					return false;
				}
				if (block.rigid_vars[rigid] != TwoSAT_Var_Undefined
						&& !sat.getAssignment(block.rigid_vars[rigid])) {
					block.spqr->reverse(rigid);
				}
				rigids++;
			}
			block.spqr->embed(block.subgraph);
		} else {
			bool planar = planarEmbed(block.subgraph);
			OGDF_ASSERT(planar);
		}
		OGDF_ASSERT(block.subgraph.representsCombEmbedding());
	}
	STEP("embedSPQR", << ",\"rigids\":" << rigids);

	// copy embedding back into G
	log.lout(Logger::Level::High) << "Copying embedding back to graph" << endl;
	for (node n : G->nodes) {
		List<adjEntry> new_adj;
		if (components.isCutVertex(n)) {
			OGDF_ASSERT(!partitions.isQVertex(n));
			for (adjEntry bc_adj : components.biconnectedComponent(n)->adjEntries) {
				BlockEmbedding& block = blocks[bc_adj->twinNode()];
				for (adjEntry adj : block.Gn_to_subgraph(n, &block)->adjEntries) {
					new_adj.pushBack(block.subgraph_to_Ge[adj->theEdge()]->getAdj(n));
				}
			}
		} else {
			BlockEmbedding& block = blocks[components.biconnectedComponent(n)];
			for (adjEntry adj : block.Gn_to_subgraph(n, &block)->adjEntries) {
				new_adj.pushBack(block.subgraph_to_Ge[adj->theEdge()]->getAdj(n));
			}
			if (partitions.isQVertex(n)) {
				if (sat.getAssignment(variables[partitions.getPartitionOf(n)])) {
					OGDF_ASSERT(compareCyclicOrder(n, new_adj) == SAME);
				} else {
					OGDF_ASSERT(compareCyclicOrder(n, new_adj) == REVERSED);
				}
			}
		}
		G->sort(n, new_adj);
	}
	OGDF_ASSERT(consistency.consistencyCheck());
	OGDF_ASSERT(G->representsCombEmbedding());
	STEP("applyEmbedding", );
	PQ_PROFILE_STOP("solveReduced")

	return true;
}
