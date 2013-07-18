/*
 * $Revision: 3643 $
 *
 * last checkin:
 *   $Author: beyer $
 *   $Date: 2013-07-07 14:00:35 +0200 (So, 07. Jul 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of Dijkstra's single source shortest path algorithm
 *
 * \author Matthias Woste
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

#ifndef OGDF_DIJKSTRA_H_
#define OGDF_DIJKSTRA_H_

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/BinaryHeap2.h>


namespace ogdf {

/*!
 * \brief Dijkstra's single source shortest path algorithm.
 *
 * This class implements Dijkstra's algorithm for computing single source shortest path
 * in (undirected or directed) graphs with proper, positive edge weights.
 * It returns a predecessor array as well as the shortest distances from the source node
 * to all others.
 */
template<typename T>
class Dijkstra {
public:

	/*!
	 * \brief Calculates, based on the graph G with corresponding edge costs and source nodes,
	 * the shortest paths and distances to all other nodes by Dijkstra's algorithm.
	 */
	void call(const Graph &G, //!< The original input graph
		  const EdgeArray<T> &weight, //!< The edge weights
		  const List<node> &sources, //!< A list of source nodes
		  NodeArray<edge> &predecessor, //!< The resulting predecessor relation
		  NodeArray<T> &distance, //!< The resulting distances to all other nodes
		  bool directed = false) //!< True iff G should be interpreted as directed graph
	{
		BinaryHeap2<T, node> queue(G.numberOfNodes());
		NodeArray<int> qpos(G);

		// initialization
		node v;
		forall_nodes(v, G) {
			distance[v] = numeric_limits<T>::max();
			predecessor[v] = NULL;
			queue.insert(v, distance[v], &qpos[v]);
		}
		forall_listiterators(node, s, sources) {
			queue.decreaseKey(qpos[*s], (distance[*s] = 0));
		}
#ifdef OGDF_DEBUG
		edge de;
		forall_edges(de, G){
			if (weight[de] <= 0) OGDF_THROW(PreconditionViolatedException);
		}
#endif

		while (!queue.empty()) {
			v = queue.extractMin();
			if (!predecessor[v] && distance[v]) { // v is unreachable, ignore
				continue;
			}
			adjEntry adj;
			forall_adj(adj, v) {
				edge e = adj->theEdge();
				node w = adj->twinNode();
				if (directed && e->target() == v) { // edge is in wrong direction
					continue;
				}
				if (distance[w] > distance[v] + weight[e]) {
					if (numeric_limits<double>::max() - weight[e] < distance[v]) cerr << "Overflow\n";
					if (-numeric_limits<double>::max() - weight[e] > distance[v]) cerr << "Overflow\n";
					queue.decreaseKey(qpos[w], (distance[w] = distance[v] + weight[e]));
					predecessor[w] = e;
				}
			}
		}
	}

	/*!
	 * \brief Calculates, based on the graph G with corresponding edge costs and a source node s,
	 * the shortest paths and distances to all other nodes by Dijkstra's algorithm.
	 */
	void call(const Graph &G, //!< The original input graph
		  const EdgeArray<T> &weight, //!< The edge weights
		  node s, //!< The source node
		  NodeArray<edge> &predecessor, //!< The resulting predecessor relation
		  NodeArray<T> &distance, //!< The resulting distances to all other nodes
		  bool directed = false) //!< True iff G should be interpreted as directed graph
	{
		List<node> sources;
		sources.pushBack(s);
		call(G, weight, sources, predecessor, distance, directed);
	}
};

} // end namespace ogdf

#endif /* OGDF_DIJKSTRA_H_ */
