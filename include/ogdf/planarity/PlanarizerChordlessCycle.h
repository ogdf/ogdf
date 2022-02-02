/** \file
 * \brief Declaration of class PlanarizerChordlessCycle.
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


namespace ogdf {

/**
 * @ingroup ga-crossmin
 *
 * Starts with a chordless cycle of the graph and then inserts each original
 * node that is adjacent to already inserted ones via the StarInserter. See:
 *
 * K. Clancy, M. Haythorpe, A. Newcombe:
 * <i>An effective crossing minimisation heuristic based on star insertion</i>.
 * J. Graph Algorithms Appl. 23(2): 135-166 (2019)
 *
 * @warning Ignores the time limit set in the CrossingMinimizationModule.
 * However, the heuristic should be very fast in practice anyway.
 */
class OGDF_EXPORT PlanarizerChordlessCycle : public CrossingMinimizationModule
{
protected:
	//! Implements the algorithm call.
	/**
	 * @pre \p pr must be simple.
	 */
	virtual ReturnType doCall(PlanRep &pr,
		int cc,
		const EdgeArray<int>      *pCostOrig,
		const EdgeArray<bool>     *pForbiddenOrig,
		const EdgeArray<uint32_t> *pEdgeSubGraphs,
		int &crossingNumber) override;

public:
	//! Creates a PlanarizerChordlessCycle with default settings.
	PlanarizerChordlessCycle();

	//! Creates a PlanarizerChordlessCycle with the same settings as \p planarizer.
	PlanarizerChordlessCycle(const PlanarizerChordlessCycle &planarizer);

	//! Returns a new PlanarizerChordlessCycle with the same option settings.
	virtual CrossingMinimizationModule *clone() const override;

	//! Assignment operator, copies option settings only.
	PlanarizerChordlessCycle &operator=(const PlanarizerChordlessCycle &planarizer);

private:
	/**
	 * @param G graph in which a chordless cycle should be found.
	 * @param cycle is assigned a list of nodes in \p G which induce a chordless
	 * cycle.
	 */
	bool findChordlessCycle(const Graph &G, List<node> &cycle);

	/**
	 * Creates crossings in \p pr that resemble the crossings in \p copyCopy.
	 * @warning the order of adjEntries around each node is not copied.
	 *
	 * @param pr is assigned the final planarization as it is contained in
	 * copyCopy in terms of crossings (but \p pr must be planarly embedded again
	 * afterwards).
	 * @param graphCopy is a copy of \p pr's original graph and the original of
	 * \p copyCopy.
	 * @param copyCopy is a planarization that contains the crossings that
	 * should be copied to \p pr.
	 */
	void transferToPlanRep(
		PlanRep &pr,
		const GraphCopy &graphCopy,
		const GraphCopy &copyCopy);

	/**
	 * Creates a copy of \p vOrig in \p graphCopy and optimally inserts a copy
	 * of this copy in the planarization \p copyCopy.
	 *
	 * @param graphCopy is the graph copy of the original graph.
	 * @param copyCopy is the graph copy/planarization of \p graphCopy.
	 * @param dual is the dual graph of \p copyCopy.
	 * @param vOrig is the node in the original graph for which a copy is
	 * created in \p graphCopy and \p copyCopy.
	 * @param pCostOrig are the edge costs in the original graph.
	 * @param pCostCopy are the edge costs in \p graphCopy (edge costs for the
	 * newly created edges are added).
	 */
	void addToGraphCopy(
		GraphCopy &graphCopy,
		GraphCopy &copyCopy,
		DynamicDualGraph &dual,
		node vOrig,
		const EdgeArray<int> *pCostOrig,
		EdgeArray<int> *pCostCopy);

	//! The StarInserter used to insert new nodes into the planarization.
	StarInserter m_inserter;
};

}
