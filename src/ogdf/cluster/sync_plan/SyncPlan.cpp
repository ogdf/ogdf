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
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/GraphSets.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/graphics.h>
#include <ogdf/basic/pctree/PCEnum.h>
#include <ogdf/basic/pctree/PCNode.h>
#include <ogdf/basic/pctree/PCTree.h>
#include <ogdf/basic/pctree/PCTreeIterators.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/cluster/sync_plan/PMatching.h>
#include <ogdf/cluster/sync_plan/QPartitioning.h>
#include <ogdf/cluster/sync_plan/SyncPlan.h>
#include <ogdf/cluster/sync_plan/SyncPlanComponents.h>
#include <ogdf/cluster/sync_plan/SyncPlanDrawer.h>
#include <ogdf/cluster/sync_plan/basic/Drawing.h>
#include <ogdf/cluster/sync_plan/utils/Bijection.h>
#include <ogdf/decomposition/BCTree.h>

#include <array>
#include <functional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <typeinfo>

using namespace ogdf::sync_plan::internal;

namespace ogdf::sync_plan {

SyncPlan::SyncPlan(Graph* g, GraphAttributes* ga)
	: G(g)
	, matchings(G)
	, partitions(G)
	, components(G)
	, deletedEdges(*G)
#ifdef OGDF_DEBUG
	, deletedNodes(*G)
#endif
	, GA(ga)
	, is_wheel(*G, false)
#ifdef OGDF_DEBUG
	, consistency(*this)
#endif
{
	initComponents();
}

void SyncPlan::initComponents() {
	OGDF_ASSERT(isLoopFree(*G));
	BCTree bc(*G, true);
	components.reset();
	components.insert(bc);
	if (bc.auxiliaryGraph().numberOfNodes() != G->numberOfNodes() + bc.bcTree().numberOfEdges()) {
		components.labelIsolatedNodes();
	}
}

void SyncPlan::formatNode(node n) const {
	if (GA == nullptr) {
		return;
	}
#ifdef OGDF_DEBUG
	if (deletedNodes.isMember(n)) {
		GA->width(n) = 0;
		GA->height(n) = 0;
		return;
	}
#endif
	ogdf::sync_plan::formatNode(n, GA, components.biconnectedId(n));
	if (matchings.isMatchedPVertex(n)) {
		GA->width(n) = 10;
		GA->height(n) = 10;
		GA->fillColor(n) = Color::Name::White;
	} else if (partitions.isQVertex(n)) {
		GA->shape(n) = Shape::Rect;
		GA->width(n) = 10;
		GA->height(n) = 10;
		GA->fillColor(n) = colors[partitions.getPartitionOf(n) % colors.size()];
	}
}

PipeType SyncPlan::getPipeType(const Pipe* p) const {
	OGDF_ASSERT(p != nullptr);
	if (components.isCutVertex(p->node1)) {
		if (components.isCutVertex(p->node2)) {
			return PipeType::CutCut;
		} else {
			return PipeType::BlockCut;
		}
	} else {
		if (components.isCutVertex(p->node2)) {
			return PipeType::BlockCut;
		} else {
			return PipeType::BlockBlock;
		}
	}
}

std::ostream& operator<<(std::ostream& os, Operation op) {
	switch (op) {
	case Operation::ENCAPSULATE_CONTRACT:
		return os << "ENCAPSULATE_CONTRACT";
	case Operation::CONTRACT_BICON:
		return os << "CONTRACT_BICON";
	case Operation::PROPAGATE_CUT:
		return os << "PROPAGATE_CUT";
	case Operation::PROPAGATE_BICON:
		return os << "PROPAGATE_BICON";
	case Operation::SIMPLIFY_TERMINAL:
		return os << "SIMPLIFY_TERMINAL";
	case Operation::SIMPLIFY_TRANSITIVE:
		return os << "SIMPLIFY_TRANSITIVE";
	case Operation::SIMPLIFY_TOROIDAL:
		return os << "SIMPLIFY_TOROIDAL";
	case Operation::BATCH_SPQR:
		return os << "BATCH_SPQR";
	default:
		return os << "Operation::???";
	}
}

int sumPNodeDegrees(const pc_tree::PCTree& pct) {
	int deg = 0;
	for (pc_tree::PCNode* node : pct.innerNodes()) {
		if (node->getNodeType() == pc_tree::PCNodeType::PNode) {
			deg += node->getDegree();
		}
	}
	return deg;
}

#ifdef SYNCPLAN_OPSTATS

void SyncPlan::printOPStatsStart(const Pipe* p, Operation op, const NodePCRotation* pct) {
	if (!stats_first_in_array) {
		stats_out << ",";
	} else {
		stats_first_in_array = false;
	}
	stats_out << "{\"op\":\"" << op << "\""
			  << ",\"rem_pipes\":" << matchings.getPipeCount() << ",\"deg\":" << p->degree()
			  << ",\"u_cv\":" << components.isCutVertex(p->node1)
			  << ",\"u_blocks\":" << components.biconnectedComponent(p->node1)->degree()
			  << ",\"u_bicon_size\":" << components.bcSize(components.biconnectedComponent(p->node1))
			  << ",\"u_bc_id\":" << components.biconnectedId(p->node1)
			  << ",\"u_cc_id\":" << components.connectedId(p->node1)
			  << ",\"v_cv\":" << components.isCutVertex(p->node2)
			  << ",\"v_blocks\":" << components.biconnectedComponent(p->node2)->degree()
			  << ",\"v_bicon_size\":" << components.bcSize(components.biconnectedComponent(p->node2))
			  << ",\"v_bc_id\":" << components.biconnectedId(p->node2)
			  << ",\"v_cc_id\":" << components.connectedId(p->node2) << ",";
	if (pct != nullptr) {
		if (stats_pc_time != 0) {
			stats_out << "\"pc_time_ns\":" << stats_pc_time << ",";
			stats_pc_time = 0;
		}
		stats_out << "\"p_nodes\":" << pct->getPNodeCount()
				  << ",\"c_nodes\":" << pct->getCNodeCount()
				  << ",\"p_node_degs\":" << sumPNodeDegrees(*pct) << ",";
	}
}

void SyncPlan::printOPStatsEnd(bool success, int64_t time_ns) {
	if (!success) {
		stats_out << "\"suc\":false,";
	}
	stats_out << "\"op_time_ns\":" << time_ns << "}";
}

#endif

bool SyncPlan::canContract(const Pipe* p) const {
	if (p == nullptr) {
		return false;
	}
	if (components.nodeType(p->node1) != components.nodeType(p->node2)) {
		return false;
	}
	if (components.nodeType(p->node1) == BCTree::GNodeType::Normal) {
		if (allow_contract_bb_pipe) {
			return components.connectedId(p->node1) != components.connectedId(p->node2);
		} else {
			return false;
		}
	} else {
		return true;
	}
}

std::ostream& operator<<(std::ostream& os, const SyncPlan& pq) {
	os << "SyncPlan Instance with " << pq.G->numberOfNodes() << " nodes, " << pq.G->numberOfEdges()
	   << " edges, " << pq.matchings.getPipeCount() << " pipes";
	PipeQueue* queue = pq.matchings.getPipeQueue();
	if (queue) {
		os << " (max degree " << queue->getTop()->degree() << ")";
	}
	return os << ", " << pq.partitions.qVertexCount() << " Q-Vertices in "
			  << pq.partitions.partitionCount() << " partitions and "
			  << pq.components.connectedCount() << " connected components";
}

std::ostream& operator<<(std::ostream& os, const SyncPlan::UndoOperation& undo_op) {
	return undo_op.print(os);
}

std::function<std::ostream&(std::ostream&)> SyncPlan::fmtPQNode(node n, bool include_comp) const {
	OGDF_ASSERT(n == nullptr || n->graphOf() == G);
	return [n, include_comp, this](std::ostream& ss) -> std::ostream& {
		if (n != nullptr) {
			ss << "["
			   << (matchings.isMatchedPVertex(n) ? "mP" : (partitions.isQVertex(n) ? "Q" : "uP"));
			if (include_comp) {
				ss << (components.isCutVertex(n) ? "C" : "B");
			}
			ss << " #" << n->index() << " °" << n->degree();
			if (include_comp) {
				ss << " @" << components.biconnectedId(n) << "/" << components.connectedId(n);
			}
			if (GA != nullptr) {
				ss << " \"" << GA->label(n) << "\"";
			}
			ss << "]";
		} else {
			ss << "[NULL]";
		}
		return ss;
	};
}

SyncPlan::VerifyPipeBijections::VerifyPipeBijections(SyncPlan& pq) {
	for (const Pipe& pipe : pq.matchings) {
		std::tuple<int, int, FrozenPipeBij>& tuple = *pipes.emplaceBack();
		std::get<0>(tuple) = pipe.node1->index();
		std::get<1>(tuple) = pipe.node2->index();
		PipeBij bij;
		pq.matchings.getIncidentEdgeBijection(pipe.node1, bij);
		freezePipeBijection(bij, std::get<2>(tuple));
	}
}

void SyncPlan::VerifyPipeBijections::undo(SyncPlan& pq) {
	if (pq.matchings.getPipeCount() != pipes.size()) {
		throw std::runtime_error("number of pipes changed!");
	}
	for (std::tuple<int, int, FrozenPipeBij>& tuple : pipes) {
		node u = pq.nodeFromIndex(std::get<0>(tuple));
		node v = pq.nodeFromIndex(std::get<1>(tuple));
		if (pq.matchings.getTwin(u) != v) {
			throw std::runtime_error("pipe endpoint changed!");
		}
		if (!pq.verifyPipeBijection(u, v, std::get<2>(tuple))) {
			throw std::runtime_error("pipe bijection broke!");
		}
	}
}

std::ostream& SyncPlan::VerifyPipeBijections::print(std::ostream& os) const {
	return os << "VerifyPipeBijections";
}

SyncPlan::ResetIndices::ResetIndices(SyncPlan& pq)
	: max_node(pq.G->maxNodeIndex())
	, max_edge(pq.G->maxEdgeIndex())
	, count_node(pq.G->numberOfNodes())
	, count_edge(pq.G->numberOfEdges()) {
	OGDF_ASSERT(!pq.indices_saved);
	pq.indices_saved = true;
}

void SyncPlan::ResetIndices::undo(SyncPlan& pq) {
	if (count_node != pq.G->numberOfNodes()) {
		throw std::runtime_error("number of nodes changed!");
	}
	if (count_edge != pq.G->numberOfEdges()) {
		throw std::runtime_error("number of edges changed!");
	}
	pq.G->resetNodeIdCount(max_node);
	pq.G->resetEdgeIdCount(max_edge);
	pq.indices_saved = false;
}

std::ostream& SyncPlan::ResetIndices::print(std::ostream& os) const {
	return os << "ResetIndices(max_node " << max_node << ", max_edge " << max_edge
			  << ", count_node " << count_node << ", count_edge " << count_edge << ")";
}

std::string SyncPlan::UndoOperation::name() const { return typeid(*this).name(); }

}
