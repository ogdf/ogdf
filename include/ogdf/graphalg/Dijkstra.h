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
#include <ogdf/basic/EpsilonTest.h>
#include <ogdf/basic/PriorityQueue.h>


namespace ogdf {

/*!
 * \brief %Dijkstra's single source shortest path algorithm.
 *
 * @ingroup ga-sp
 *
 * This class implements Dijkstra's algorithm for computing single source shortest path in
 * (undirected or directed) graphs with proper, positive edge weights.
 * It returns a predecessor array as well as the shortest distances from the source node(s)
 * to all others.
 * It optionally supports early termination if only the shortest path to a specific node
 * is required, or the maximum path length is to be limited.
 */
template<typename T, template<typename P, class C> class H = PairingHeap>
class Dijkstra {
protected:
	EpsilonTest m_eps; //!< For floating point comparisons (if floating point is used)

public:
	//! Calculates, based on the graph G with corresponding edge costs and source nodes,
	//! the shortest paths and distances to all other nodes by Dijkstra's algorithm.
	/**
	 * @param G The original input graph
	 * @param weight The edge weights
	 * @param sources A list of source nodes
	 * @param predecessor The resulting predecessor relation
	 * @param distance The resulting distances to all other nodes
	 * @param directed True iff G should be interpreted as a directed graph
	 * @param arcsReversed True if the arcs should be followed in reverse. It has only
	 * an effect when setting \p directed to true
	 */
	void callUnbound(const Graph &G,
		  const EdgeArray<T> &weight,
		  const List<node> &sources,
		  NodeArray<edge> &predecessor,
		  NodeArray<T> &distance,
		  bool directed = false,
		  bool arcsReversed = false)
	{
		PrioritizedMapQueue<node, T, std::less<T>, H> queue(G);
		distance.init(G, std::numeric_limits<T>::max());
		predecessor.init(G, nullptr);

#ifdef OGDF_DEBUG
		// No source should be given multiple times.
		NodeArray<bool> isSource {G, false};
		for (node s : sources) {
			OGDF_ASSERT(!isSource[s]);
			isSource[s] = true;
		}
#endif

		// initialization
		for (node v : G.nodes) {
			queue.push(v, distance[v]);
		}
		for (node s : sources) {
			queue.decrease(s, (distance[s] = 0));
		}

#ifdef OGDF_DEBUG
		for (edge de : G.edges) {
			OGDF_ASSERT(weight[de] >= 0);
		}
#endif

		while (!queue.empty()) {
			node v = queue.topElement();
			queue.pop();
			if (!predecessor[v]
			 && m_eps.greater(distance[v], static_cast<T>(0))) { // v is unreachable, ignore
				continue;
			}
			for(adjEntry adj : v->adjEntries) {
				edge e = adj->theEdge();
				node w = adj->twinNode();
				if (directed && ((!arcsReversed && e->target() == v) || (arcsReversed && e->target() != v))) {
					continue;
				}

				if (m_eps.greater(distance[w], distance[v] + weight[e])) {
					OGDF_ASSERT(std::numeric_limits<T>::max() - weight[e] >= distance[v]);
					queue.decrease(w, (distance[w] = distance[v] + weight[e]));
					predecessor[w] = e;
				}
			}
		}
	}

	//! @copybrief ::callUnbound(const Graph&, const EdgeArray<T>&, const List<node>&, NodeArray<edge>&, NodeArray<T>&, bool, bool)
	/**
	 * Allows to specify a target node and maximum distance, after reaching which the algorithm
	 * will terminate early.
	 *
	 * This implementation is different from the implementation of callUnbound() as runtime tests have shown the
	 * additional checks to increase run time for the basic use case.
	 *
	 * @copydetails ::callUnbound(const Graph&, const EdgeArray<T>&, const List<node>&, NodeArray<edge>&, NodeArray<T>&, bool, bool)
	 *
	 * @param target A target node. Terminate once the shortest path to this node is found
	 * @param maxLength Upper bound on path length
	 */
	void callBound(const Graph &G,
		  const EdgeArray<T> &weight,
		  const List<node> &sources,
		  NodeArray<edge> &predecessor,
		  NodeArray<T> &distance,
		  bool directed,
		  bool arcsReversed,
		  node target,
		  T maxLength = std::numeric_limits<T>::max())
	{
		PrioritizedMapQueue<node, T, std::less<T>, H> queue(G);
		distance.init(G, std::numeric_limits<T>::max());
		predecessor.init(G, nullptr);

		// initialization
		for (node s : sources) {
			queue.push(s, (distance[s] = 0));
		}

#ifdef OGDF_DEBUG
		for (edge de : G.edges) {
			OGDF_ASSERT(weight[de] >= 0);
		}
#endif

		while (!queue.empty()) {
			node v = queue.topElement();
			if (v == target) { // terminate early if this is our sole target
				break;
			}

			queue.pop();
			if (!predecessor[v]
			 && m_eps.greater(distance[v], static_cast<T>(0))) { // v is unreachable, ignore
				continue;
			}
			for(adjEntry adj : v->adjEntries) {
				edge e = adj->theEdge();
				node w = adj->twinNode();
				if (directed && ((!arcsReversed && e->target() == v) || (arcsReversed && e->target() != v))) {
					continue;
				}

				const T newDistance = distance[v] + weight[e];
				if (m_eps.greater(newDistance, maxLength)) {
					// using this edge would result in a path length greater than our upper bound
					continue;
				}
				if (m_eps.greater(distance[w], newDistance)) {
					OGDF_ASSERT(std::numeric_limits<T>::max() - weight[e] >= distance[v]);
					distance[w] = newDistance;
					if (queue.contains(w)) {
						queue.decrease(w, distance[w]);
					} else {
						queue.push(w, distance[w]);
					}
					predecessor[w] = e;
				}
			}
		}
	}

	//! Calculates, based on the graph G with corresponding edge costs and a source node s,
	//! the shortest paths and distances to all other nodes by Dijkstra's algorithm.
	/**
	 * @param G The original input graph
	 * @param weight The edge weights
	 * @param s The source node
	 * @param predecessor The resulting predecessor relation
	 * @param distance The resulting distances to all other nodes
	 * @param directed True iff G should be interpreted as a directed graph
	 * @param arcsReversed True if the arcs should be followed in reverse. It has only
	 * an effect when setting \p directed to true
	 */
	void callUnbound(const Graph &G,
		  const EdgeArray<T> &weight,
		  node s,
		  NodeArray<edge> &predecessor,
		  NodeArray<T> &distance,
		  bool directed = false,
		  bool arcsReversed = false)
	{
		List<node> sources;
		sources.pushBack(s);
		callUnbound(G, weight, sources, predecessor, distance, directed, arcsReversed);
	}

	//! @copybrief ::callUnbound(const Graph&, const EdgeArray<T>&, node, NodeArray<edge>&, NodeArray<T>&, bool, bool)
	/**
	 * Allows to specify a target node and maximum distance, after reaching which the algorithm
	 * will terminate early.
	 *
	 * This implementation is different from the implementation of callUnbound() as runtime tests have shown the
	 * additional checks to increase run time for the basic use case.
	 *
	 * @copydetails ::callUnbound(const Graph&, const EdgeArray<T>&, node, NodeArray<edge>&, NodeArray<T>&, bool, bool)
	 * @param target A target node. Terminate once the shortest path to this node is found
	 * @param maxLength Upper bound on path length
	 */
	void callBound(const Graph &G,
		  const EdgeArray<T> &weight,
		  node s,
		  NodeArray<edge> &predecessor,
		  NodeArray<T> &distance,
		  bool directed,
		  bool arcsReversed,
		  node target,
		  T maxLength = std::numeric_limits<T>::max())
	{
		List<node> sources;
		sources.pushBack(s);
		callBound(G, weight, sources, predecessor, distance, directed, arcsReversed, target, maxLength);
	}

	//! @copydoc ::callUnbound(const Graph&, const EdgeArray<T>&, const List<node>&, NodeArray<edge>&, NodeArray<T>&, bool, bool)
	/**
	 * @param target A target node. Terminate once the shortest path to this node is found
	 * @param maxLength Upper bound on path length
	 *
	 * @note If no target or maximum distance is given, use the unbound algorithm that runs faster
	 * on most instances. On some types of instances (especially sparse ones) the bound algorithm
	 * tends to run faster. To force its usage, use the \c callBound method directly.
	 *
	 * @see callBound(const Graph&, const EdgeArray<T>&, const List<node>&, NodeArray<edge>&, NodeArray<T>&, bool, node, T)
	 */
	void call(const Graph &G,
		  const EdgeArray<T> &weight,
		  const List<node> &sources,
		  NodeArray<edge> &predecessor,
		  NodeArray<T> &distance,
		  bool directed = false,
		  bool arcsReversed = false,
		  node target = nullptr,
		  T maxLength = std::numeric_limits<T>::max())
	{
		if (target == nullptr && maxLength == std::numeric_limits<T>::max()) {
			callUnbound(G, weight, sources, predecessor, distance, directed, arcsReversed);
		}
		else {
			callBound(G, weight, sources, predecessor, distance, directed, arcsReversed, target, maxLength);
		}
	}

	//! @copydoc ::callUnbound(const Graph&, const EdgeArray<T>&, node, NodeArray<edge>&, NodeArray<T>&, bool, bool)
	/**
	 * @param target A target node. Terminate once the shortest path to this node is found
	 * @param maxLength Upper bound on path length
	 *
	 * @note If no target or maximum distance is given, use the unbound algorithm that runs faster
	 * on most instances. On some types of instances (especially sparse ones) the bound algorithm
	 * tends to run faster. To force its usage, use the \c callBound method directly.
	 *
	 * @see callBound(const Graph&, const EdgeArray<T>&, node, NodeArray<edge>&, NodeArray<T>&, bool, node, T)
	 */
	void call(const Graph &G,
		  const EdgeArray<T> &weight,
		  const node s,
		  NodeArray<edge> &predecessor,
		  NodeArray<T> &distance,
		  bool directed = false,
		  bool arcsReversed = false,
		  node target = nullptr,
		  T maxLength = std::numeric_limits<T>::max())
	{
		if (target == nullptr && maxLength == std::numeric_limits<T>::max()) {
			callUnbound(G, weight, s, predecessor, distance, directed, arcsReversed);
		}
		else {
			callBound(G, weight, s, predecessor, distance, directed, arcsReversed, target, maxLength);
		}
	}

};

}
