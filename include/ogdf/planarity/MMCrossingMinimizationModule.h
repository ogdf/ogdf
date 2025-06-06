/** \file
 * \brief Declaration of MMCrossingMinimization Module, an interface
 * for minor-monotone crossing minimization algorithms.
 *
 * \author Carsten Gutwenger
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

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/Module.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/memory.h>

namespace ogdf {
class PlanRepExpansion;
template<class E>
class List;

/**
 * \brief Interface for minor-monotone crossing minimization algorithms.
 *
 */
class OGDF_EXPORT MMCrossingMinimizationModule : public Module {
public:
	//! Initializes a minor-monotone crossing minimization module.
	MMCrossingMinimizationModule() : m_nodeSplits(0), m_splittedNodes(0) { }

	// destruction
	virtual ~MMCrossingMinimizationModule() { }

	/**
	 * \brief Computes a planarized representation of an expansion of the input graph.
	 *
	 * @param PG represents the input graph as well as the computed planarized
	 *        expansion after the call. \p PG has to be initialzed as a
	 *        PlanRepExpansion of the input graph and is modified to obatain the planarized
	 *        representation (nodes are eventually expanded by splitting the node and
	 *        crossings are replaced by dummy vertices with degree four).
	 * @param cc is the number of the connected component in \p PG that is considered.
	 * @param crossingNumber is assigned the number of crossings.
	 * @param forbid points to an edge array indicating which edges are not allowed
	 *        to be crossed, i.e., (*forbid)[e] = true. If forbid = 0, no edges are
	 *        forbidden.
	 * \return the status of the result.
	 */
	ReturnType call(PlanRepExpansion& PG, int cc, int& crossingNumber,
			const EdgeArray<bool>* forbid = nullptr) {
		return doCall(PG, cc, forbid, crossingNumber, m_nodeSplits, m_splittedNodes);
	};

	/**
	 * \brief Performs minor-monotone crossing minimization on \p G.
	 *
	 * @param G is the input graph.
	 * @param cr is assigned the number of crossings.
	 * @param forbid points to an edge array indicating which edges are not allowed
	 *        to be crossed, i.e., (*forbid)[e] = true. If forbid = 0, no edges are
	 *        forbidden.
	 * \return the status of the result.
	 */
	ReturnType call(const Graph& G, int& cr, const EdgeArray<bool>* forbid = nullptr);

	/**
	 * \brief Performs minor-monotone crossing minimization on \p G for given splittable nodes.
	 *
	 * @param G is the input graph.
	 * @param splittableNodes is the list of nodes that are allowed to be split.
	 * @param cr is assigned the number of crossings.
	 * @param forbid points to an edge array indicating which edges are not allowed
	 *        to be crossed, i.e., (*forbid)[e] = true. If forbid = 0, no edges are
	 *        forbidden.
	 * \return the status of the result.
	 */
	ReturnType call(const Graph& G, const List<node>& splittableNodes, int& cr,
			const EdgeArray<bool>* forbid = nullptr);

	/**
	 * \brief Returns the number of required node splits after the call.
	 */
	int numberOfNodeSplits() const { return m_nodeSplits; }

	int numberOfSplittedNodes() const { return m_splittedNodes; }

protected:
	/**
	 * \brief Actual algorithm call that needs to be implemented by derived classed.
	 *
	 * @param PG represents the input graph as well as the computed planarized expansion
	 *        after the call. \p PG is initialized as a PlanRepExpansion of the input
	 *        graph and needs to be modified to obatain the planarized representation
	 *        (crossings are replaced by dummy vertices with degree four).
	 * @param cc is the number of the connected component in \p PG that is considered.
	 * @param forbid points to an edge array indicating which edges are not allowed
	 *        to be crossed, i.e., (*forbid)[e] = true.
	 * @param crossingNumber needs to be assigned the number of crossings.
	 * @param numNS needs to be assigned the required number of node splits.
	 * @param numSN needs to be assigned the number of splitted nodes.
	 * \return the status of the result.
	 */
	virtual ReturnType doCall(PlanRepExpansion& PG, int cc, const EdgeArray<bool>* forbid,
			int& crossingNumber, int& numNS, int& numSN) = 0;

private:
	int m_nodeSplits; //!< The number of required node splits.
	int m_splittedNodes; //!< The number of nodes that are split.

	OGDF_MALLOC_NEW_DELETE
};

}
