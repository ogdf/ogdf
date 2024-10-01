/** \file
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
#include <ogdf/basic/Array.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/GraphSets.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/Logger.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/comparer.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/cluster/sync_plan/PMatching.h>
#include <ogdf/cluster/sync_plan/QPartitioning.h>
#include <ogdf/cluster/sync_plan/SyncPlan.h>
#include <ogdf/cluster/sync_plan/SyncPlanComponents.h>
#include <ogdf/cluster/sync_plan/SyncPlanConsistency.h>
#include <ogdf/cluster/sync_plan/SyncPlanDrawer.h>
#include <ogdf/cluster/sync_plan/basic/GraphUtils.h>
#include <ogdf/cluster/sync_plan/utils/Bijection.h>
#include <ogdf/cluster/sync_plan/utils/Logging.h>
#include <ogdf/decomposition/BCTree.h>
#include <ogdf/fileformats/GraphIO.h>

#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>

using namespace ogdf::sync_plan::internal;

namespace ogdf::sync_plan {
using internal::operator<<;

bool SyncPlanConsistency::doWriteOut = false;

void normalize(List<adjEntry>& adjs) {
	adjEntry min_adj = adjs.front();
	for (adjEntry adj : adjs) {
		if (adj->theEdge()->index() < min_adj->theEdge()->index()) {
			min_adj = adj;
		}
	}
	while (adjs.front() != min_adj) {
		adjs.pushBack(adjs.popFrontRet());
	}
}

void SyncPlanConsistency::writeOut(std::string name, bool format, bool components) {
	if (name.empty()) {
		std::stringstream ss;
		ss << "consistencyCheck" << checkCounter;
		name = ss.str();
	}
	draw.layout(format, components);
	if (pq.GA != nullptr) {
		std::stringstream ss;
		ss << name << ".svg";
		std::ofstream os(ss.str());
		GraphIO::drawSVG(*pq.GA, os, draw.getSvg());
	}
	if (components) {
		std::stringstream ss;
		ss << name << "Tree.svg";
		std::ofstream os(ss.str());
		GraphIO::drawSVG(draw.getBC_GA(), os, draw.getSvg());
	}
	draw.cleanUp();
	{
		std::stringstream ss;
		ss << name << ".gml";
		std::ofstream os(ss.str());
		if (pq.GA != nullptr) {
			GraphIO::writeGML(*pq.GA, os);
		} else {
			GraphIO::writeGML(*pq.G, os);
		}
	}
	if (components) {
		std::stringstream ss;
		ss << name << "Tree.gml";
		std::ofstream os(ss.str());
		GraphIO::writeGML(draw.getBC_GA(), os);
	}
#if 0
	{
		std::stringstream ss;
		ss << name << ".json";
		std::ofstream os(ss.str());
		nlohmann::json json;
		SyncPlanOptions::generateConfigJSON(pq, json);
		os << json.dump(4) << std::endl;
	}
#endif
	{
		std::stringstream ss;
		ss << name << ".edges.txt";
		std::ofstream os(ss.str());
		List<edge> edges;
		pq.G->allEdges(edges);
		TargetComparer<EdgeElement, EdgeElement> comp;
		edges.quicksort(comp);
		for (edge e : edges) {
			os << "#" << e->index() << " (" << e->source() << ", " << e->target() << ")" << std::endl;
		}
	}
	{
		std::stringstream ss;
		ss << name << ".nodes.txt";
		std::ofstream os(ss.str());
		List<node> nodes;
		pq.G->allNodes(nodes);
		TargetComparer<NodeElement, NodeElement> comp;
		nodes.quicksort(comp);
		for (node n : nodes) {
			os << n << std::endl;
		}
	}
	{
		std::stringstream ss;
		ss << name << ".pipes.txt";
		std::ofstream os(ss.str());
		List<node> pipe_nodes;
		for (const Pipe& p : pq.matchings) {
			pipe_nodes.pushBack((p.node1->index() > p.node2->index()) ? p.node2 : p.node1);
		}
		TargetComparer<NodeElement, NodeElement> comp;
		pipe_nodes.quicksort(comp);
		for (node pipe_node : pipe_nodes) {
			os << "(" << pipe_node << ", " << pq.matchings.getTwin(pipe_node) << ") {";
			PipeBij bij;
			pq.matchings.getIncidentEdgeBijection(pipe_node, bij);
			bij.quicksort(PipeBijCmp());
			for (auto edges : bij) {
				os << "(" << edges.first->theEdge()->index() << ", "
				   << edges.second->theEdge()->index() << "), ";
			}
			os << "}" << std::endl;
		}
	}
	{
		std::stringstream ss;
		ss << name << ".partitions.txt";
		std::ofstream os(ss.str());
		for (int i = 0; i < pq.partitions.partitionCount(); i++) {
			if (pq.partitions.nodesInPartition(i).empty()) {
				continue;
			}
			List<node> nodes(pq.partitions.nodesInPartition(i));
			TargetComparer<NodeElement, NodeElement> comp;
			nodes.quicksort(comp);
			int flip = 0;
			for (node n : nodes) {
				List<adjEntry> adjs;
				n->allAdjEntries(adjs);
				normalize(adjs);
				if (flip == 0) {
					flip = (*adjs.begin().succ())->theEdge()->index()
									> adjs.back()->theEdge()->index()
							? -1
							: 1;
				}
				if (flip == -1) {
					adjs.reverse();
					adjs.pushFront(adjs.popBackRet());
					OGDF_ASSERT(compareCyclicOrder(n, adjs) == OrderComp::REVERSED);
				} else {
					OGDF_ASSERT(compareCyclicOrder(n, adjs) == OrderComp::SAME);
				}
				os << "p" << i << " n" << n->index() << " " << printIncidentEdges(adjs) << std::endl;
			}
		}
	}
	pq.log.lout(Logger::Level::High) << ">>>>> Wrote " << name << "(Tree).gml/svg <<<<<" << std::endl;
}

bool SyncPlanConsistency::consistencyCheck(bool force_check_components) {
	if (doWriteOut) {
		writeOut();
	}

#ifdef OGDF_DEBUG
	pq.G->consistencyCheck();
	pq.components.bcTree().consistencyCheck();

	NodeArray<node> node_reg(*pq.G, nullptr);
	NodeArray<node> bc_reg(pq.components.BC, nullptr);
	for (node bc : pq.components.BC.nodes) {
		bc_reg[bc] = bc;
	}

	int q_count = 0, p_count = 0;
	for (node n : pq.G->nodes) {
		if (pq.deletedNodes.isMember(n)) {
			continue;
		}
		node_reg[n] = n;
		OGDF_ASSERT(!(pq.matchings.isMatchedPVertex(n) && pq.partitions.isQVertex(n)));
		// OGDF_ASSERT(!(pq.partitions.isQVertex(n) && pq.components.isCutVertex(n) && n->degree() > 3));
		if (pq.partitions.isQVertex(n)) {
			q_count++;
		}
		if (pq.matchings.isMatchedPVertex(n)) {
			p_count++;
		}
		pq.matchings.getTwinOrNull(n);
		node bicon = pq.components.biconnectedComponent(n);
		OGDF_ASSERT(bc_reg[bicon] == bicon);
		int conn_id = pq.components.bcConnectedId(bicon);
		OGDF_ASSERT(conn_id >= 0);
		for (adjEntry adj : n->adjEntries) {
			node adj_bicon = pq.components.biconnectedComponent(adj->twinNode());
			if (bicon != adj_bicon) {
				OGDF_ASSERT(conn_id == pq.components.bcConnectedId(adj_bicon));
				OGDF_ASSERT(pq.components.isCutComponent(bicon)
						|| pq.components.isCutComponent(adj_bicon));
				const std::pair<edge, edge>& bcEdge =
						pq.components.graphEdgeToBCEdge(bicon, adj_bicon);
				if (pq.components.isCutComponent(bicon) && pq.components.isCutComponent(adj_bicon)) {
					OGDF_ASSERT(bcEdge.second != nullptr);
					node common = bcEdge.first->commonNode(bcEdge.second);
					OGDF_ASSERT(common != nullptr);
					OGDF_ASSERT(!pq.components.isCutComponent(common));
				}
			}
		}
	}
	OGDF_ASSERT(q_count == pq.partitions.qVertexCount());
	OGDF_ASSERT(p_count == pq.matchings.getPipeCount() * 2);

	for (node bc : pq.components.BC.nodes) {
		node bc_repr_n = pq.components.bcRepr(bc);
		BCTree::BNodeType type = pq.components.bcNodeType(bc);
		OGDF_ASSERT(bc_repr_n == nullptr || node_reg[bc_repr_n] == bc_repr_n);
		OGDF_ASSERT(type == BCTree::BNodeType::BComp
				|| (type == BCTree::BNodeType::CComp && bc->degree() >= 2));
	}

	q_count = 0;
	for (int part = 0; part < pq.partitions.partitionCount(); part++) {
		for (node u : pq.partitions.nodesInPartition(part)) {
			if (pq.deletedNodes.isMember(u)) {
				continue;
			}
			OGDF_ASSERT(node_reg[u] == u);
			OGDF_ASSERT(pq.partitions.getPartitionOf(u) == part);
			// OGDF_ASSERT(!pq.components.isCutVertex(u) || u->degree() <= 3);
			q_count++;
		}
	}
	OGDF_ASSERT(q_count == pq.partitions.qVertexCount());

	// SimplePipeQueue *pipeCmp = nullptr;
	if (pq.matchings.queue) {
		OGDF_ASSERT(pq.matchings.pipes_list.size() == pq.matchings.queue->size());
		// pipeCmp = (SimplePipeQueue *) pq.matchings.queue.get();
	}
	// const Pipe* biggest_pipe = nullptr;
	for (const auto& pipe : pq.matchings) {
		OGDF_ASSERT(node_reg[pipe.node1] == pipe.node1);
		OGDF_ASSERT(node_reg[pipe.node2] == pipe.node2);
		OGDF_ASSERT(pipe.node1 != pipe.node2);
		OGDF_ASSERT(pq.matchings.getTwin(pipe.node1) == pipe.node2);
		// OGDF_ASSERT(pq.matchings.pipes_heap.value(pipe.heap_entry) == &pipe);
		// if (biggest_pipe == nullptr || (pipeCmp && (*pipeCmp)(&pipe, biggest_pipe)))
		//     biggest_pipe = &pipe;
	}
	//    OGDF_ASSERT(*pq.matchings.getTopPipe().heap_entry == 1);
	//     if (!pq.matchings.isReduced() && pipeCmp) {
	// const Pipe *top = &pq.matchings.getTopPipe();
	// OGDF_ASSERT(pipeCmp->checkOrder(biggest_pipe, top));
	// pq.matchings.rebuildHeap();
	// const Pipe *top2 = &pq.matchings.getTopPipe();
	// OGDF_ASSERT(top2 == top || !(*pipeCmp)(top2, top));
	// }

	if (force_check_components || checkCounter % 10 == 0) { // save time by only sometimes checking
		checkComponentRegeneration();
	}
#endif
	checkCounter++;
	return true;
}

void SyncPlanConsistency::checkComponentRegeneration() {
#ifdef OGDF_DEBUG
	BCTree ref_bc(*pq.G, true);
	NodeArray<int> ref_conn(ref_bc.bcTree(), 0);
	int ref_conn_count = connectedComponents(ref_bc.bcTree(), ref_conn);
	{
		NodeArray<bool> bc_seen(ref_bc.bcTree(), false);
		NodeArray<node> bc_store(ref_bc.bcTree(), nullptr);
		NodeArray<node> bc_ref(ref_bc.bcTree(), nullptr);

		for (node g_n : pq.G->nodes) {
			if (g_n->degree() == 0) {
				if (!pq.deletedNodes.isMember(g_n)) {
					node bc = pq.components.biconnectedComponent(g_n);
					OGDF_ASSERT(bc->degree() == 0);
					OGDF_ASSERT(pq.components.bcRepr(bc) == g_n);
					OGDF_ASSERT(pq.components.bc_size[bc] == 1);
					OGDF_ASSERT(!pq.components.isCutComponent(bc));
				}
				continue;
			}
			node v1 = ref_bc.bcproper(g_n);
			node v2 = pq.components.biconnectedComponent(g_n);
			if (bc_seen[v1]) {
				if (bc_store[v1] != v2) {
					pq.log.lout(Logger::Level::Alarm)
							<< "From node " << pq.fmtPQNode(bc_ref[v1]) << " I learned the mapping "
							<< "{{" << ref_bc.typeOfBNode(v1) << " #" << v1->index() << " °"
							<< v1->degree() << " @" << ref_conn(v1) << "}}"
							<< " => " << pq.components.fmtBCNode(bc_store[v1]) << ". "
							<< "For node " << pq.fmtPQNode(g_n)
							<< " I got the same key, but the value now is "
							<< pq.components.fmtBCNode(v2) << "." << std::endl;
					OGDF_ASSERT(false);
				}
			} else {
				OGDF_ASSERT(ref_bc.numberOfNodes(v1) == pq.components.bcSize(v2));
				bc_seen[v1] = true;
				bc_store[v1] = v2;
				bc_ref[v1] = g_n;
			}
		}
	}
	OGDF_ASSERT(ref_bc.bcTree().numberOfNodes()
			== pq.components.bcTree().numberOfNodes() + pq.deletedNodes.size());
	{
		Array<bool> conn_seen(0, ref_conn_count - 1, false);
		Array<int> conn_store(0, ref_conn_count - 1, -1);
		Array<node> conn_ref(0, ref_conn_count - 1, nullptr);

		for (node g_n : pq.G->nodes) {
			if (g_n->degree() == 0) {
				continue;
			}
			int v1 = ref_conn[ref_bc.bcproper(g_n)];
			int v2 = pq.components.connectedId(g_n);
			if (conn_seen[v1]) {
				if (conn_store[v1] != v2) {
					pq.log.lout(Logger::Level::Alarm)
							<< "From node " << pq.fmtPQNode(conn_ref[v1])
							<< " I learned the mapping " << v1 << " => " << conn_store[v1] << ". "
							<< "For node " << pq.fmtPQNode(g_n)
							<< " I got the same key, but the value now is " << v2 << ".";
					OGDF_ASSERT(false);
				}
			} else {
				conn_seen[v1] = true;
				conn_store[v1] = v2;
				conn_ref[v1] = g_n;
			}
		}
	}
	OGDF_ASSERT(ref_conn_count == pq.components.connectedCount() + pq.deletedNodes.size());
#endif
}

}
