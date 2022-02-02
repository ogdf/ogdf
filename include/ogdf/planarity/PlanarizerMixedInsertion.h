/** \file
 * \brief Declaration of class PlanarizerMixedInsertion.
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
#include <ogdf/planarity/PlanarSubgraphModule.h>


namespace ogdf {

/**
 * @ingroup ga-crossmin
 *
 * Computes a planar subgraph of the graph and then re-inserts each original
 * node that is incident to at least one edge not in the subgraph via the
 * StarInserter. Whether all nodes incident to these edges are re-inserted or
 * just some of them, can be configured. If an edge is not in the planar
 * subgraph but both of its endpoints are cut vertices in the planar subgraph,
 * the edge is reinserted via an EdgeInsertionModule. See:
 *
 * M. Chimani, M. Ilsen, T. Wiedera:
 * <i>Star-Struck by Fixed Embeddings: Modern Crossing Number Heuristics.</i>
 * In: Purchase H.C., Rutter I. (eds)  GD 2021.
 * LNCS, vol. 12868, pp. 41--56. Springer, Cham (2021).
 * https://doi.org/10.1007/978-3-030-92931-2_3
 *
 * <H3>Options</H3>
 *
 * <table>
 *   <tr>
 *     <th><i>Option</i><th><i>Type</i><th><i>Default</i><th><i>Description</i>
 *   </tr><tr>
 *     <td><i>subgraph</i><td>PlanarSubgraphModule<td>PlanarSubgraphFast<int>
 *     <td>The module for the computation of the planar subgraph.
 *   </tr><tr>
 *     <td><i>nodeSelectionMethod</i><td>NodeSelectionMethod<td>HigherDegree
 *     <td>Determines which nodes are reinserted.
 *   </tr>
 * </table>
 */
class OGDF_EXPORT PlanarizerMixedInsertion : public CrossingMinimizationModule
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
	// (The algorithm may still work when pr is not simple and the checking
	// assertion is removed as long as only edge insertion is used.)

public:
	//! Determines the node(s) of each deleted edge e which will be reinserted
	//! if neither of them is a cut vertex.
	enum class NodeSelectionMethod {
		Random, //!< Either the source or the target of e, decided randomly.
		HigherDegree, //!< The endpoint of e with the higher degree.
		LowerDegree, //!< The endpoint of e with the lower degree.
		HigherNonPlanarDegree, //!< The endpoint of e with the higher number of
		                       //!< incident edges not in the planar subgraph.
		LowerNonPlanarDegree, //!< The endpoint of e with the lower number of
		                      //!< incident edges not in the planar subgraph.
		BothEndpoints //!< Both the source and the target of e.
	};

	//! Creates a PlanarizerMixedInsertion with default settings.
	PlanarizerMixedInsertion();

	//! Creates a PlanarizerMixedInsertion with the same settings as \p planarizer.
	PlanarizerMixedInsertion(const PlanarizerMixedInsertion &planarizer);

	//! Returns a new PlanarizerMixedInsertion with the same option settings.
	virtual CrossingMinimizationModule *clone() const override;

	//! Assignment operator, copies option settings only.
	PlanarizerMixedInsertion &operator=(const PlanarizerMixedInsertion &planarizer);

	//! Sets the module option for the computation of the planar subgraph.
	void setSubgraph(PlanarSubgraphModule<int> *pSubgraph) {
		m_subgraph.reset(pSubgraph);
	}

	//! Returns the used method of selecting nodes to reinsert.
	NodeSelectionMethod nodeSelectionMethod() { return m_nodeSelectionMethod; }

	//! Sets the used method of selecting nodes to reinsert.
	void nodeSelectionMethod(NodeSelectionMethod method) {
		m_nodeSelectionMethod = method;
	}

private:
	//!< The planar subgraph algorithm.
	std::unique_ptr<PlanarSubgraphModule<int>> m_subgraph;

	//!< How to select the nodes which are reinserted.
	NodeSelectionMethod m_nodeSelectionMethod;
};

}
