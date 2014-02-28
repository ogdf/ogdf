/*
 * $Revision: 3941 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2014-02-28 11:33:37 +0100 (Fr, 28. Feb 2014) $
 ***************************************************************/

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
 * See README.txt in the root directory of the OGDF installation for details.
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
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/

#ifdef _MSC_VER
#pragma once
#endif

#ifndef OGDF_SHORTEST_PATH_ALGORITHMS_H_
#define OGDF_SHORTEST_PATH_ALGORITHMS_H_

#include <ogdf/basic/SList.h>
#include <ogdf/basic/Graph_d.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/NodeArray.h>
#include <ogdf/basic/Array.h>
#include <ogdf/basic/BinaryHeap2.h>


namespace ogdf {

//! BFS to compute shortest path all pairs. The costs for
//! traversing an edge corresponds to /a edgeCosts.
OGDF_EXPORT
void bfs_SPAP(const Graph& G, NodeArray<NodeArray<double> >& distance,
		double edgeCosts);

//! BFS to compute shortest path single source. The costs for
//! traversing an edge corresponds to /a edgeCosts.
OGDF_EXPORT
void bfs_SPSS(const node& v, const Graph& G, NodeArray<double> & distanceArray,
		double edgeCosts);

//! Dijkstra algorithm to compute shortest path all pairs. The costs for traversing edge e
//! corresponds to \a GA.doubleWeight(e)
/**
 * @return returns the average edge costs
 */
OGDF_EXPORT
double dijkstra_SPAP(const GraphAttributes& GA,
		NodeArray<NodeArray<double> >& shortestPathMatrix);

//! Dijkstra algorithm to compute shortest path all pairs. The costs for traversing edge e
//! corresponds to \a edgeCosts[e]
OGDF_EXPORT
void dijkstra_SPAP(const Graph& G,
		NodeArray<NodeArray<double> >& shortestPathMatrix,
		const EdgeArray<double>& edgeCosts);

//! Dijkstra algorithm to compute shortest path single source. The costs for traversing edge e
//! corresponds to \a edgeCosts[e]. Note this algorithm equals Dijkstra<T>::call, though it does not
//! compute the predecessors on the path and is not inlined.
OGDF_EXPORT
void dijkstra_SPSS(node v, const Graph& G,
		NodeArray<double>& shortestPathMatrix,
		const EdgeArray<double>& edgeCosts);

//! Floyd-Wharshall algorithm to compute shortest path all pairs given a weighted graph.
//! Note the shortestPathMatrix has to be initialized and all entries positive. The costs
//! non-adjacent nodes should be set to std::numeric_limits<double>::infinity().
OGDF_EXPORT
void floydWarshall_SPAP(NodeArray<NodeArray<double> >& shortestPathMatrix,
		const Graph& G);

} /* namespace ogdf */
#endif /* SHORTESTPATHALGORITHMS_H_ */
