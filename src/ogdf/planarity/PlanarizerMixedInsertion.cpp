/** \file
 * \brief Implementation of class PlanarizerMixedInsertion.
 *
 * \author Max Ilsen
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

#include <ogdf/planarity/PlanarizerMixedInsertion.h>
#include <ogdf/planarity/StarInserter.h>
#include <ogdf/planarity/PlanarSubgraphFast.h>
#include <ogdf/planarity/VariableEmbeddingInserter.h>
#include <ogdf/planarity/embedder/CrossingStructure.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <set>

namespace ogdf {

PlanarizerMixedInsertion::PlanarizerMixedInsertion()
	: CrossingMinimizationModule()
	, m_nodeSelectionMethod{NodeSelectionMethod::HigherDegree}
{
	PlanarSubgraphFast<int> *subgraph = new PlanarSubgraphFast<int>();
	m_subgraph.reset(subgraph);
}

PlanarizerMixedInsertion::PlanarizerMixedInsertion(
	const PlanarizerMixedInsertion &planarizer)
	: CrossingMinimizationModule()
{
	m_subgraph.reset(planarizer.m_subgraph->clone());
	m_nodeSelectionMethod = planarizer.m_nodeSelectionMethod;
}

CrossingMinimizationModule *PlanarizerMixedInsertion::clone() const {
	return new PlanarizerMixedInsertion(*this);
}

PlanarizerMixedInsertion &PlanarizerMixedInsertion::operator=(
	const PlanarizerMixedInsertion &planarizer)
{
	m_subgraph.reset(planarizer.m_subgraph->clone());
	m_nodeSelectionMethod = planarizer.m_nodeSelectionMethod;
	return *this;
}

Module::ReturnType PlanarizerMixedInsertion::doCall(
	PlanRep &pr,
	int cc,
	const EdgeArray<int> *pCostOrig,
	const EdgeArray<bool> *pForbiddenOrig,
	const EdgeArray<uint32_t> *pEdgeSubGraphs,
	int &crossingNumber)
{
	OGDF_ASSERT(isSimpleUndirected(pr));
	crossingNumber = 0;

	int64_t startTime;
	System::usedRealTime(startTime);
	int64_t stopTime {m_timeLimit >= 0 ?
		startTime + int64_t(1000 * m_timeLimit) : -1};

	pr.initCC(cc);

	// Compute edges to delete for planar subgraph.
	m_subgraph->timeLimit(m_timeLimit);
	List<edge> delEdges;
	ReturnType retValue {pCostOrig ?
		m_subgraph->call(pr.original(), *pCostOrig, delEdges) :
		m_subgraph->call(pr.original(), delEdges)};

	if (!isSolution(retValue)) {
		return retValue;
	}

	// Delete edges such that only planar subgraph remains.
	PlanRepLight prl {pr};
	prl.initCC(cc);
	for (edge eDel : delEdges) {
		prl.removeEdgePath(eDel);
	}

#ifdef OGDF_DEBUG
	bool embeddedSuccessfully =
#endif
		planarEmbed(prl);
	OGDF_ASSERT(embeddedSuccessfully);

	// If the graph is planar, return with a planar embedding.
	if (delEdges.empty()) {
		embedder::CrossingStructure cs;
		cs.init(prl, -1);
		cs.restore(pr, cc);
		return ReturnType::Optimal;
	}

	// Each delEdge with two cut vertices as endpoints must be inserted via an
	// edge inserter (the star inserter would remove an endpoint and
	// disconnect the graph). Collect these edges in an array first.
	ArrayBuffer<node> cutVertices;
	findCutVertices(prl, cutVertices);
	NodeArray<bool> isCutVertex {prl, false};
	for (node v : cutVertices) {
		isCutVertex[v] = true;
	}

	SListPure<edge> cutVertexEdges;
	ListIterator<edge> itSucc;
	for (ListIterator<edge> it = delEdges.begin(); it.valid(); it = itSucc) {
		itSucc = it.succ();
		edge eDel {*it};
		if (isCutVertex[prl.copy(eDel->source())] &&
		    isCutVertex[prl.copy(eDel->target())]) {
			cutVertexEdges.pushBack(eDel);
			delEdges.del(it);
		}
	}
	Array<edge> edgeInsertionEdges{cutVertexEdges.size()};
	int edgeNum {0};
	for (edge e : cutVertexEdges) {
		edgeInsertionEdges[edgeNum++] = e;
	}

	// Insert edges whose endpoints are both cut vertices.
	// We must use RemoveReinsertType::None since it otherwise tries to use
	// postprocessing for edges which are not yet in the graph (i.e. delEdges).
	VariableEmbeddingInserter edgeInserter;
	edgeInserter.removeReinsert(RemoveReinsertType::None);
	ReturnType ret {edgeInserter.callEx(prl, edgeInsertionEdges,
		pCostOrig, pForbiddenOrig, pEdgeSubGraphs)};
	if (!isSolution(ret)) {
		return ret;
	}

	// If needed: For each node, collect number of incident edges not in the
	// planar subgraph.
	NodeArray<int> nonPlanarDegree;
	if (m_nodeSelectionMethod == NodeSelectionMethod::HigherNonPlanarDegree ||
		m_nodeSelectionMethod == NodeSelectionMethod::LowerNonPlanarDegree) {
		nonPlanarDegree.init(prl.original(), 0);
		for (edge eDel : delEdges) {
			nonPlanarDegree[eDel->source()]++;
			nonPlanarDegree[eDel->target()]++;
		}
	}

	// Get stars from deleted edges (cheap vertex cover).
	SListPure<node> delNodes;
	std::set<node> delNodesSet;
	for (edge eDel : delEdges) {
		node src {eDel->source()};
		node tgt {eDel->target()};

		// If one of the nodes, the other one has to be chosen.
		if (isCutVertex[prl.copy(src)]) {
			OGDF_ASSERT(!isCutVertex[prl.copy(tgt)]);
			delNodesSet.insert(tgt);
			continue;
		}

		if (isCutVertex[prl.copy(tgt)]) {
			OGDF_ASSERT(!isCutVertex[prl.copy(src)]);
			delNodesSet.insert(src);
			continue;
		}

		// Neither node is a cut vertex.
		switch (m_nodeSelectionMethod) {
			case NodeSelectionMethod::Random:
				delNodesSet.insert(randomNumber(0, 1) == 0 ? src : tgt);
				break;

			case NodeSelectionMethod::HigherDegree:
				delNodesSet.insert(src->degree() > tgt->degree() ? src : tgt);
				break;

			case NodeSelectionMethod::LowerDegree:
				delNodesSet.insert(src->degree() < tgt->degree() ? src : tgt);
				break;

			case NodeSelectionMethod::HigherNonPlanarDegree:
				delNodesSet.insert(nonPlanarDegree[src] > nonPlanarDegree[tgt] ?
					src : tgt);
				break;

			case NodeSelectionMethod::LowerNonPlanarDegree:
				delNodesSet.insert(nonPlanarDegree[src] < nonPlanarDegree[tgt] ?
					src : tgt);
				break;

			default:
				delNodesSet.insert(src);
				delNodesSet.insert(tgt);
		}
	}

	ArrayBuffer<node> nodesToReinsert;
	for (node v : delNodesSet) {
		nodesToReinsert.push(v);
	}

	// Reinsert stars.
	StarInserter inserter;
	CombinatorialEmbedding emb {prl};
	DynamicDualGraph dualGraph {emb};
	for (node nodeToReinsert : nodesToReinsert) {
		// Remove remaining edge paths of edges incident to nodeToReinsert.
		for (adjEntry adj : nodeToReinsert->adjEntries) {
			edge e {adj->theEdge()};
			if (!prl.chain(e).empty()) {
				prl.GraphCopy::removeEdgePathEmbedded(emb, dualGraph, e);
			}
		}

		// Reinsert nodeToReinsert.
		inserter.call(prl, dualGraph, nodeToReinsert, pCostOrig);

		// If the time limit is reached before all deleted nodes could be
		// inserted, the current solution is infeasible.
		if (stopTime >= 0 && System::realTime() >= stopTime) {
			return ReturnType::TimeoutInfeasible;
		}
	}

	// Restore planarization in pr.
	embedder::CrossingStructure cs;
	cs.init(prl, -1);
	cs.restore(pr, cc);

	// Remove pseudo crossings and recompute crossing number.
#ifdef OGDF_DEBUG
	bool planar =
#endif
		planarEmbed(pr);
	OGDF_ASSERT(planar);
	pr.removePseudoCrossings();
	crossingNumber = computeCrossingNumber(pr, pCostOrig, pEdgeSubGraphs);

	OGDF_ASSERT(isPlanar(pr));
	OGDF_ASSERT(!pr.hasNonSimpleCrossings());
	return ReturnType::Feasible;
}

}
