/** \file
 * \brief Declares maximum density subgraph algorithms.
 *
 * \author Finn Stutzenstein
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
#include <ogdf/basic/NodeSet.h>

namespace ogdf {

/**
 * \brief Calculates the maximum density subgraph of \p G
 *
 * Returns the subgraph of \p G given by the nodes of the subgraph. The subgraph
 * with the highest density, which is \f$\frac{|E|}{|V|}\f$ of the subgraph, is
 * returned. The graph is treated as undirected and unweighted.
 *
 * The algorithm is based on: A. V. Goldberg. Finding a maximum density
 * subgraph. EECS Department, University of California, Berkeley. Technical Report
 * No. UCB/CSD-84-171, 1984.
 * https://www2.eecs.berkeley.edu/Pubs/TechRpts/1984/CSD-84-171.pdf
 *
 * Internally, the MinSTCutMaxFlow algorithm (using MaxFlowGoldbergTarjan) is
 * called \f$\mathcal{O}(\log n)\f$ times. Assuming a runtime of
 * \f$\mathcal{O}(mn^2)\f$ for the min cut algorithm, the overall runtime is
 * \f$\mathcal{O}(mn^2\log n)\f$.
 *
 * \par Side Effect
 * Note that the input graph \p G is modified. If this is not wanted, use a
 * GraphCopy or GraphCopySimple.
 *
 * \par resultNodeMap
 * \p resultNodeMap can be used if \p subgraphNodes does not belong to \p G.
 * This helps when passing in a GraphCopy or GraphCopySimple:
 * \code
 * NodeSet<true> subgraphNodes(G);
 * GraphCopySimple GC(G);
 * maximumDensitySubgraph(GC, subgraphNodes, [&](node n){ return GC.original(n);});
 * \endcode
 *
 * @ingroup ga-induced
 *
 * @param G is the input graph.
 * @param subgraphNodes is the set of nodes inducing the subgraph.
 * @param resultNodeMap maps each subgraph node (nodes of G) to some other node
 * @param timelimit set to a value greater than -1 to set a timelimit in milliseconds.
 *                  Note that 0 is a valid timelimit. When encountering a timelimit there is no valid result.
 *
 * @returns true, if the algorithm was successful and did not run into a timeout.
 */
OGDF_EXPORT bool maximumDensitySubgraph(
		Graph& G, NodeSet<true>& subgraphNodes,
		std::function<node(node)> resultNodeMap = [](node v) { return v; }, int64_t timelimit = -1);

}
