/** \file
 * \brief Implementation of the SyncPlan::solve operation
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
#include <ogdf/basic/basic.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/cluster/sync_plan/PMatching.h>
#include <ogdf/cluster/sync_plan/QPartitioning.h>
#include <ogdf/cluster/sync_plan/SyncPlan.h>
#include <ogdf/cluster/sync_plan/SyncPlanComponents.h>
#include <ogdf/cluster/sync_plan/SyncPlanConsistency.h>
#include <ogdf/cluster/sync_plan/SyncPlan_solve/BlockEmbedding.h>
#include <ogdf/cluster/sync_plan/basic/GraphUtils.h>
#include <ogdf/cluster/sync_plan/basic/TwoSAT.h>
#include <ogdf/decomposition/Skeleton.h>
#include <ogdf/decomposition/StaticPlanarSPQRTree.h>

#include <memory>
#include <ostream>

using namespace ogdf::sync_plan::internal;

namespace ogdf::sync_plan {

#ifdef SYNCPLAN_OPSTATS
#	define SYNCPLAN_OPSTATS_STEP(op, meta)                                                        \
		stats_out << (stats_first_in_array ? "" : ",") << "{\"op\":\"solvedReduced-" << op << "\"" \
				  << ",\"op_time_ns\":" << dur_ns(tpc::now() - start) meta << "}";                 \
		stats_first_in_array = false;                                                              \
		start = tpc::now()
#else
#	define SYNCPLAN_OPSTATS_STEP(op, meta)
#endif

bool SyncPlan::solveReduced(bool fail_fast) {
	// SYNCPLAN_PROFILE_START("solveReduced")
#ifdef SYNCPLAN_OPSTATS
	std::chrono::time_point<std::chrono::high_resolution_clock> start = tpc::now();
#endif
	OGDF_ASSERT(matchings.isReduced());
	// ensure that all Q-node are surrounded by wheels
	// room for improvement: makeWheel could also be replaced by a Q-vertex-aware embedding tree generator
#ifdef SYNCPLAN_OPSTATS
	int wheels = 0;
#endif
	for (int part = 0; part < partitions.partitionCount(); part++) {
		for (node u : partitions.nodesInPartition(part)) {
			if (!is_wheel[u]) {
				// also replace for u->degree() == 3, otherwise multiple following assertions break
				makeWheel(u);
#ifdef SYNCPLAN_OPSTATS
				wheels++;
#endif
			}
		}
	}
	SYNCPLAN_OPSTATS_STEP("makeWheels", << ",\"wheels\":" << wheels);
	if (fail_fast && !isPlanar(*G)) {
		// SYNCPLAN_PROFILE_STOP("solvedReduced")
		return false;
	}
	log.lout(Logger::Level::High) << "SOLVE REDUCED" << std::endl;

	// generate SPQR trees for all blocks
	log.lout(Logger::Level::High) << "Generating SPQR Trees for up to "
								  << components.bcTree().numberOfNodes() << " blocks" << std::endl;
	GnMultiArray Gn_to_subgraph(*G);
	NodeArrayP<BlockEmbedding> blocks(components.bcTree());
	for (node n : components.bcTree().nodes) {
		blocks[n] = std::make_unique<BlockEmbedding>(BlockEmbedding(Gn_to_subgraph));
	}
	EdgeArray<edge> Ge_to_subgraph(*G, nullptr);
	EdgeArray<BlockEmbedding*> Ge_to_block(*G, nullptr);
#ifdef SYNCPLAN_OPSTATS
	int block_cnt = 0;
#endif
	for (node block : components.bcTree().nodes) {
		if (!components.isCutComponent(block)) {
			blocks[block]->init(*G, components, block, Ge_to_subgraph, Ge_to_block);
#ifdef SYNCPLAN_OPSTATS
			block_cnt++;
#endif
		}
	}
	SYNCPLAN_OPSTATS_STEP("deriveSPQR", << ",\"blocks\":" << block_cnt);

	// generate 2-SAT instance
	log.lout(Logger::Level::High) << "Generating 2-SAT instances for "
								  << partitions.partitionCount() << " partitions" << std::endl;
	TwoSAT sat;
	PartitionArray<twosat_var> variables(partitions, TwoSAT_Var_Undefined);
#ifdef SYNCPLAN_OPSTATS
	int q_vertices = 0;
#endif
	for (int part = 0; part < partitions.partitionCount(); part++) {
		twosat_var part_var = variables[part] = sat.newVariable();
		for (node u : partitions.nodesInPartition(part)) {
			OGDF_ASSERT(partitions.isQVertex(u));
			OGDF_ASSERT(!components.isCutVertex(u));
#ifdef SYNCPLAN_OPSTATS
			q_vertices++;
#endif
			if (!blocks[components.biconnectedComponent(u)]->addQVertex(u, Ge_to_subgraph, sat,
						part_var)) {
				// SYNCPLAN_PROFILE_STOP("solvedReduced")
				return false;
			}
		}
	}
	SYNCPLAN_OPSTATS_STEP("deriveSAT", << ",\"q_vertices\":" << q_vertices);

	log.lout(Logger::Level::High) << "Solving 2-SAT..." << std::endl;
	if (!sat.solve()) {
		// SYNCPLAN_PROFILE_STOP("solvedReduced")
		return false;
	}
	SYNCPLAN_OPSTATS_STEP("solveSAT", );

	// embed rigids according to 2-SAT result
	log.lout(Logger::Level::High) << "Embedding rigids according to 2-SAT solution" << std::endl;
#ifdef SYNCPLAN_OPSTATS
	int rigids = 0;
#endif
	for (node bc : components.bcTree().nodes) {
		BlockEmbedding& block = *blocks[bc];
		if (block.spqr) {
			for (node rigid : block.spqr->tree().nodes) {
				if (!fail_fast && !isPlanar(block.spqr->skeleton(rigid).getGraph())) {
					// SYNCPLAN_PROFILE_STOP("solvedReduced")
					return false;
				}
				if (block.rigid_vars[rigid] != TwoSAT_Var_Undefined
						&& !sat.getAssignment(block.rigid_vars[rigid])) {
					block.spqr->reverse(rigid);
				}
#ifdef SYNCPLAN_OPSTATS
				rigids++;
#endif
			}
			block.spqr->embed(block.subgraph);
		} else {
			bool planar = planarEmbed(block.subgraph);
			OGDF_ASSERT(planar);
		}
		OGDF_ASSERT(block.subgraph.representsCombEmbedding());
	}
	SYNCPLAN_OPSTATS_STEP("embedSPQR", << ",\"rigids\":" << rigids);

	// copy embedding back into G
	log.lout(Logger::Level::High) << "Copying embedding back to graph" << std::endl;
	for (node n : G->nodes) {
		if (n->degree() < 2) {
			continue;
		}
		OGDF_ASSERT(!deletedNodes.isMember(n));
		List<adjEntry> new_adj;
		if (components.isCutVertex(n)) {
			OGDF_ASSERT(!partitions.isQVertex(n));
			for (adjEntry bc_adj : components.biconnectedComponent(n)->adjEntries) {
				BlockEmbedding& block = *blocks[bc_adj->twinNode()];
				for (adjEntry adj : block.Gn_to_subgraph(n, &block)->adjEntries) {
					new_adj.pushBack(block.subgraph_to_Ge[adj->theEdge()]->getAdj(n));
				}
			}
		} else {
			BlockEmbedding& block = *blocks[components.biconnectedComponent(n)];
			for (adjEntry adj : block.Gn_to_subgraph(n, &block)->adjEntries) {
				new_adj.pushBack(block.subgraph_to_Ge[adj->theEdge()]->getAdj(n));
			}
			if (partitions.isQVertex(n)) {
				if (sat.getAssignment(variables[partitions.getPartitionOf(n)])) {
					OGDF_ASSERT(compareCyclicOrder(n, new_adj) == OrderComp::SAME);
				} else {
					OGDF_ASSERT(compareCyclicOrder(n, new_adj) == OrderComp::REVERSED);
				}
			}
		}
		G->sort(n, new_adj);
	}
	OGDF_ASSERT(consistency.consistencyCheck());
	OGDF_ASSERT(G->representsCombEmbedding());
	SYNCPLAN_OPSTATS_STEP("applyEmbedding", );
	// SYNCPLAN_PROFILE_STOP("solveReduced")

	return true;
}

}
