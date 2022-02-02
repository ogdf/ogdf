/** \file
 * \brief Implementation of class PlanarizerStarReinsertion.
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

#include <ogdf/planarity/PlanarizerStarReinsertion.h>
#include <ogdf/planarity/SubgraphPlanarizer.h>
#include <ogdf/planarity/FixedEmbeddingInserter.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/simple_graph_alg.h>

namespace ogdf {

PlanarizerStarReinsertion::PlanarizerStarReinsertion()
	: m_inserter{}
	, m_setTimeout{true}
	, m_maxIterations{-1}
	, m_stopTime{-1}
{
	SubgraphPlanarizer* heuristic = new SubgraphPlanarizer();
	heuristic->setInserter(new FixedEmbeddingInserter);
	heuristic->permutations(1);
	m_planarization.reset(heuristic);
}

PlanarizerStarReinsertion::PlanarizerStarReinsertion(const PlanarizerStarReinsertion &planarizer)
	: CrossingMinimizationModule(planarizer)
	, m_inserter{}
	, m_setTimeout{planarizer.m_setTimeout}
	, m_maxIterations{planarizer.m_maxIterations}
	, m_stopTime{-1}
{
	m_planarization.reset(planarizer.m_planarization->clone());
}

CrossingMinimizationModule *PlanarizerStarReinsertion::clone() const {
	return new PlanarizerStarReinsertion(*this);
}

PlanarizerStarReinsertion &PlanarizerStarReinsertion::operator=(const PlanarizerStarReinsertion &planarizer)
{
	m_timeLimit = planarizer.m_timeLimit;
	m_planarization.reset(planarizer.m_planarization->clone());

	m_setTimeout = planarizer.m_setTimeout;
	m_maxIterations = planarizer.m_maxIterations;
	m_stopTime = -1;
	return *this;
}

bool PlanarizerStarReinsertion::reinsertStar(
	GraphCopy &currentPlanarization,
	DynamicDualGraph &dualGraph,
	node nodeToReinsert,
	CrossingStructure &bestCS,
	const EdgeArray<int> *pCostOrig,
	const EdgeArray<bool> *pForbiddenOrig,
	const EdgeArray<uint32_t> *pEdgeSubGraphs)
{
	// Delete previous insertion paths of star from planarization.
	for (adjEntry adj : nodeToReinsert->adjEntries) {
		currentPlanarization.removeEdgePathEmbedded(
			dualGraph.getPrimalEmbedding(), dualGraph, adj->theEdge());
	}

	// Reinsert the node.
	m_inserter.call(currentPlanarization, dualGraph, nodeToReinsert, pCostOrig);
	OGDF_ASSERT(!currentPlanarization.hasNonSimpleCrossings());

	// Compute the new number of crossings.
	int cr {computeCrossingNumber(currentPlanarization, pCostOrig, pEdgeSubGraphs)};
	OGDF_ASSERT(cr <= bestCS.weightedCrossingNumber());

	bool crChange {cr != bestCS.weightedCrossingNumber()};
	bestCS.init(currentPlanarization, cr);
	return crChange;
}

Module::ReturnType PlanarizerStarReinsertion::mainLoop(
	const PlanRep &pr,
	CrossingStructure &bestCS,
	const EdgeArray<int> *pCostOrig,
	const EdgeArray<bool> *pForbiddenOrig,
	const EdgeArray<uint32_t> *pEdgeSubGraphs)
{
	GraphCopy currentPlanarization {pr};
	CombinatorialEmbedding emb {currentPlanarization};
	DynamicDualGraph dual {emb};
	bool solutionUpdated {true};

	int seed {rand()};
	std::minstd_rand rng(seed);
	SListPure<node> nodesToReinsert;
	pr.allNodes(nodesToReinsert);

	// Try to reinsert every node until there is no improvement possible
	// anymore, i.e. until the drawing is locally crossing-optimal, or the
	// maximum number of iterations is reached.
	int i {0};
	while (solutionUpdated && i != m_maxIterations) {
		i++;
		solutionUpdated = false;
		nodesToReinsert.permute(rng);

		for (node nodeToReinsert : nodesToReinsert) {
			if (!pr.isDummy(nodeToReinsert)) {
				nodeToReinsert = pr.original(nodeToReinsert);
				solutionUpdated = reinsertStar(currentPlanarization, dual,
					nodeToReinsert, bestCS, pCostOrig, pForbiddenOrig,
					pEdgeSubGraphs);

				if (m_stopTime >= 0 && System::realTime() >= m_stopTime) {
					return ReturnType::TimeoutFeasible;
				}

				if (solutionUpdated) {
					break;
				}
			}
		}
	}

	return ReturnType::Feasible;
}

Module::ReturnType PlanarizerStarReinsertion::doCall(
	PlanRep &pr,
	int cc,
	const EdgeArray<int> *pCostOrig,
	const EdgeArray<bool> *pForbiddenOrig,
	const EdgeArray<uint32_t> *pEdgeSubGraphs,
	int &crossingNumber)
{
	OGDF_ASSERT(isBiconnected(pr));
	OGDF_ASSERT(isSimpleUndirected(pr));

	int64_t startTime;
	System::usedRealTime(startTime);
	m_stopTime = m_timeLimit >= 0 ? startTime + int64_t(1000*m_timeLimit) : -1;

	// Get initial (usually bad but cheap) planarization of the graph.
	if (m_setTimeout) {
		m_planarization->timeLimit(m_timeLimit);
	}
	pr.initCC(cc);
	m_planarization->call(pr, cc, crossingNumber,
		pCostOrig, pForbiddenOrig, pEdgeSubGraphs);
	pr.removeNonSimpleCrossings();
	OGDF_ASSERT(!pr.hasNonSimpleCrossings());
	if (crossingNumber == 0) {
		return ReturnType::Optimal;
	}

	// Main loop.
	CrossingStructure bestCS;
	bestCS.init(pr, crossingNumber);
	ReturnType retVal {mainLoop(pr, bestCS, pCostOrig, pForbiddenOrig,
			pEdgeSubGraphs)};
	if (!isSolution(retVal)) {
		return retVal;
	}

	// Restore best solution in pr.
	// Remove pseudo crossings and recompute crossing number.
	pr.init(pr.original());
	bestCS.restore(pr, cc);
#ifdef OGDF_DEBUG
	bool planar =
#endif
		planarEmbed(pr);
	OGDF_ASSERT(planar);
	pr.removePseudoCrossings();
	crossingNumber = computeCrossingNumber(pr, pCostOrig, pEdgeSubGraphs);

	OGDF_ASSERT(isPlanar(pr));
	return retVal;
}

}
