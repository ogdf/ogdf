/** \file
 * \brief Declaration of several shortest path algorithms.
 *
 * \author Mark Ortmann, University of Konstanz
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

#include <ogdf/basic/SList.h>
#include <ogdf/basic/Graph_d.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/NodeArray.h>
#include <ogdf/basic/Array.h>


namespace ogdf {

//! Computes all-pairs shortest paths in \p G using breadth-first serach (BFS).
/**
 * @ingroup ga-sp
 *
 * The cost of each edge are \p edgeCost and the result is stored in \p distance.
 */
OGDF_EXPORT
void bfs_SPAP(const Graph& G, NodeArray<NodeArray<double> >& distance,
		double edgeCosts);


//! Computes single-source shortest paths from \p s in \p G using breadth-first search (BFS).
/**
 * @ingroup ga-sp
 *
 * The cost of each edge are \p edgeCost and the result is stored in \p distanceArray.
 */
OGDF_EXPORT
void bfs_SPSS(node s, const Graph& G, NodeArray<double> & distanceArray, double edgeCosts);


//! Computes all-pairs shortest paths in \p GA using Dijkstra's algorithm.
/**
 * @ingroup ga-sp
 *
 * The cost of an edge \a e are given by GA.doubleWeight(\a e) and the result is stored in \p shortestPathMatrix.
 *
 * @return returns the average edge cost
 */
OGDF_EXPORT
double dijkstra_SPAP(const GraphAttributes& GA, NodeArray<NodeArray<double> >& shortestPathMatrix);


//! Computes all-pairs shortest paths in graph \p G using Dijkstra's algorithm.
/**
 * @ingroup ga-sp
 *
 * The cost of an edge are given by \p edgeCosts and the result is stored in \p shortestPathMatrix.
 */
OGDF_EXPORT
void dijkstra_SPAP(
	const Graph& G,
	NodeArray<NodeArray<double> >& shortestPathMatrix,
	const EdgeArray<double>& edgeCosts);


//! Computes single-source shortest paths from node \p s in \p G using Disjkstra's algorithm.
/**
 * @ingroup ga-sp
 *
 * The cost of an edge are given by \p edgeCosts and the result is stored in \p shortestPathMatrix.
 * Note this algorithm equals Dijkstra<T>::call, though it does not
 * compute the predecessors on the path and is not inlined.
 */
OGDF_EXPORT
void dijkstra_SPSS(
	node s,
	const Graph& G,
	NodeArray<double>& shortestPathMatrix,
	const EdgeArray<double>& edgeCosts);


//! Computes all-pairs shortest paths in graph \p G using Floyd-Warshall's algorithm.
/**
 * @ingroup ga-sp
 *
 * Note that the \p shortestPathMatrix has to be initialized and all entries must be positive.
 * The costs of non-adjacent nodes should be set to std::numeric_limits<double>::infinity().
 */
OGDF_EXPORT
void floydWarshall_SPAP(NodeArray<NodeArray<double> >& shortestPathMatrix, const Graph& G);


} /* namespace ogdf */
