/** \file
 * \brief Declaration of class PlanarizerStarReinsertion.
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

#pragma once

#include <ogdf/planarity/CrossingMinimizationModule.h>
#include <ogdf/planarity/StarInserter.h>
#include <ogdf/planarity/embedder/CrossingStructure.h>
#include <memory>


namespace ogdf {

using embedder::CrossingStructure;

//! The star (re-)insertion approach for crossing minimization.
/**
 * @ingroup ga-crossmin
 *
 * This crossing minimization module represents a customizable implementation of
 * the star insertion approach. This approach consists of two phases.
 * In the first phase, a planarization of the graph is computed, i.e. a
 * planarized representation with crossings replaced by dummy nodes of degree 4.
 * In the second phase, a set of original nodes (each with its star insertion
 * paths) is re-inserted with as few crossings as possible, until no improvement
 * is possible anymore (the configuration is "locally optimal" or a maximum
 * number of iterations is reached.
 *
 * The first step, i.e. the computation of the planarization, is implemented via
 * a module option while the maximum number of iterations can be set directly.
 *
 * More details on the star insertion approach can be found in
 *
 * K. Clancy, M. Haythorpe, A. Newcombe:
 * <i>An effective crossing minimisation heuristic based on star insertion</i>.
 * J. Graph Algorithms Appl. 23(2): 135-166 (2019)
 *
 * <H3>Optional parameters</H3>
 *
 * <table>
 *   <tr>
 *     <th><i>Option</i><th><i>Type</i><th><i>Default</i><th><i>Description</i>
 *   </tr><tr>
 *     <td><i>setTimeout</i><td>bool<td>true
 *     <td>If set to true, the time limit is also passed to submodules;
 *     otherwise, a timeout might be checked late when a submodule requires a
 *     lot of runtime.
 *   </tr><tr>
 *     <td><i>maxIterations</i><td>int<td>-1
 *     <td>Maximum number of iterations. If negative, the algorithm stops when
 *     local optimality is reached. See #maxIterations(int).
 *   </tr>
 * </table>
 *
 * <H3>%Module options</H3>
 * The various phases of the algorithm can be exchanged by setting
 * module options allowing flexible customization. The algorithm provides
 * the following module options:
 *
 * <table>
 *   <tr>
 *     <th><i>Option</i><th><i>Type</i><th><i>Default</i><th><i>Description</i>
 *   </tr><tr>
 *     <td><i>planarization</i><td>CrossingMinimizationModule
 *     <td>SubgraphPlanarizer using FixedEmbeddingInserter
 *     <td>The module for the computation of the planar subgraph.
 *   </tr>
 * </table>
*/
class OGDF_EXPORT PlanarizerStarReinsertion : public CrossingMinimizationModule
{
private:
	//! The initial planarization algorithm.
	std::unique_ptr<CrossingMinimizationModule> m_planarization;

	StarInserter m_inserter; //! Helper to insert stars.

	bool m_setTimeout; //!< The option for setting timeouts in submodules.

	int m_maxIterations; //!< The maximum number of iterations.

	int64_t m_stopTime; //!< When the algorithm should stop.

	/**
	 * Reinserts the star \p nodeToReinsert in \p currentPlanarization.
	 * If the reinsertion leads to a lower crossing number, update \p bestCS.
	 *
	 * @param currentPlanarization planarized representation of the graph.
	 * @param dualGraph dual graph of \p currentPlanarization.
	 * @param nodeToReinsert the node to be reinserted.
	 * @param bestCS is (assigned) a representation of \p currentPlanarization
	 * with the least crossings found so far as well as its crossing number.
	 * @param pCostOrig points to the cost of each original edge.
	 * @param pForbiddenOrig points to an indicator for each original edge
	 * whether it is forbidden to be crossed.
	 * @param pEdgeSubGraphs points to an indicator for each original edge
	 * to which subgraph it belongs.
	 * @returns whether a better solution was found, i.e. whether
	 * \p bestCS was updated.
	 */
	bool reinsertStar(
		GraphCopy &currentPlanarization,
		DynamicDualGraph &dualGraph,
		node nodeToReinsert,
		CrossingStructure &bestCS,
		const EdgeArray<int> *pCostOrig,
		const EdgeArray<bool> *pForbiddenOrig,
		const EdgeArray<uint32_t> *pEdgeSubGraphs);

	//! Reinserts a specific set of stars until a termination criterion is met.
	/**
	 * Either the timout or the maximum number of iterations is reached, or the
	 * solution is locally optimal.
	 *
	 * @param pr planarized representation of the graph.
	 * @param bestCS is assigned a representation of \p currentPlanarization
	 * with the least crossings found during the loop and its crossing number.
	 * @param pCostOrig points to the cost of each original edge.
	 * @param pForbiddenOrig points to an indicator for each original edge
	 * whether it is forbidden to be crossed.
	 * @param pEdgeSubGraphs points to an indicator for each original edge
	 * to which subgraph it belongs.
	 * @return the return type of the algorithm.
	 */
	ReturnType mainLoop(const PlanRep &pr,
		CrossingStructure &bestCS,
		const EdgeArray<int> *pCostOrig,
		const EdgeArray<bool> *pForbiddenOrig,
		const EdgeArray<uint32_t> *pEdgeSubGraphs);

protected:
	//! Implements the algorithm call.
	/**
	 * @pre \p pr must be biconnected, simple, and already embedded, with pseudo
	 * crossings being removed.
	 */
	virtual ReturnType doCall(PlanRep &pr,
		int cc,
		const EdgeArray<int>      *pCostOrig,
		const EdgeArray<bool>     *pForbiddenOrig,
		const EdgeArray<uint32_t> *pEdgeSubGraphs,
		int &crossingNumber) override;

public:
	//! Creates a PlanarizerStarReinsertion with default settings.
	PlanarizerStarReinsertion();

	//! Creates a PlanarizerStarReinsertion with the same settings as \p planarizer.
	PlanarizerStarReinsertion(const PlanarizerStarReinsertion &planarizer);

	//! Returns a new PlanarizerStarReinsertion with the same option settings.
	virtual CrossingMinimizationModule *clone() const override;

	//! Assignment operator, copies option settings only.
	PlanarizerStarReinsertion &operator=(const PlanarizerStarReinsertion &planarizer);

	//! Sets the module option for the computation of the inital planarization.
	void setPlanarization(CrossingMinimizationModule *pPlanarizationModule) {
		m_planarization.reset(pPlanarizationModule);
	}

	//! @name Optional parameters
	//! @{

	//! Returns the current setting of options <i>setTimeout</i>.
	bool setTimeout() { return m_setTimeout; }

	//! Sets the option <i>setTimeout</i> to \p b.
	void setTimeout(bool b) { m_setTimeout = b; }

	//! Returns the number of maxIterations.
	int maxIterations() { return m_maxIterations; }

	//! Sets the maximum number of iterations.
	/**
	 * The algorithm terminates when the drawing is locally crossing-optimal: It
	 * tries to re-insert all of the nodes of the graph one after another; it
	 * stops when there is no node anymore whose re-insertion would improve the
	 * number of crossings.
	 *
	 * The algorithm can be stopped ealier by setting a maximum number of
	 * iterations (during each of which all reinsertable nodes are reinserted).
	 * A number lower than 0 indicates that the algorithm should not be stopped
	 * before a locally optimal solution is found.
	 */
	void maxIterations(int maxIterations) { m_maxIterations = maxIterations; }

	//! @}

};

}
