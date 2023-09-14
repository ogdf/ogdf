/** \file
 * \brief Implements (non-templated) simple matching functions
 *
 * \author Stephan Beyer, Thomas Klein
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

#include <ogdf/basic/GraphCopy.h>
#include <ogdf/graphalg/Matching.h>

#include <limits.h>

namespace ogdf {
namespace Matching {

void findMaximalMatching(const Graph& graph, ArrayBuffer<edge>& matching) {
	EdgeArray<bool> covered {graph, false};

	for (edge e : graph.edges) {
		if (!covered[e]) {
			matching.push(e);
			for (node v : e->nodes()) {
				for (adjEntry adj : v->adjEntries) {
					covered[adj->theEdge()] = true;
				}
			}
		}
	}
}

/**
 * Auxiliary BFS method for HKK max cardinality matching.
 *
 * @param pairV the partners of nodes in V
 * @param pairU the partners of nodes in U
 * @param dist the distances to each node
 * @param dummy the dummy node
 * @param U the list to start BFS at
 * @return whether an improvement could be found
 */
bool bfsMaximumCardinalityMatching(NodeArray<node>& pairV, NodeArray<node>& pairU,
		NodeArray<int>& dist, node dummy, List<node>& U) {
	SListPure<node> Q;

	for (node u : U) {
		if (pairU[u] == dummy) {
			dist[u] = 0;
			Q.pushBack(u);
		} else {
			dist[u] = std::numeric_limits<int>::max();
		}
	}
	dist[dummy] = std::numeric_limits<int>::max();

	while (!Q.empty()) {
		node u = Q.popFrontRet();
		if (dist[u] < dist[dummy]) {
			for (adjEntry adj : u->adjEntries) {
				node v = adj->twinNode();
				if (dist[pairV[v]] == std::numeric_limits<int>::max()) {
					dist[pairV[v]] = dist[u] + 1;
					Q.pushBack(pairV[v]);
				}
			}
		}
	}
	return dist[dummy] != std::numeric_limits<int>::max();
}

/**
 * Auxiliary DFS method for HKK max cardinality matching.
 *
 * @param pairV the partners of nodes in V
 * @param pairU the partners of nodes in U
 * @param dist the distances to each node
 * @param dummy the dummy node
 * @param u the node to start DFS at
 * @return whether an improvement could be found
 */
bool dfsMaximumCardinalityMatching(NodeArray<node>& pairV, NodeArray<node>& pairU,
		NodeArray<int>& dist, node dummy, node u) {
	if (u != dummy) {
		for (adjEntry adj : u->adjEntries) {
			node v = adj->twinNode();
			if (dist[pairV[v]] == dist[u] + 1) {
				if (dfsMaximumCardinalityMatching(pairV, pairU, dist, dummy, pairV[v])) {
					pairV[v] = u;
					pairU[u] = v;
					return true;
				}
			}
		}
		dist[u] = std::numeric_limits<int>::max();
		return false;
	}
	return true;
}

int findMaximumCardinalityMatching(const Graph& G, const List<node>& U, const List<node>& V,
		EdgeArray<bool>& matching) {
	// create a GraphCopy to work on and add a dummy node
	GraphCopy g(G);
	node dummy = g.newNode();

	// list for the new nodes
	List<node> gU;

	// connect all nodes to the dummy
	for (node u : U) {
		OGDF_ASSERT(u->graphOf() == &G);
		node uCopy = g.copy(u);
		gU.pushBack(uCopy);
		g.newEdge(uCopy, dummy);
	}

#ifdef OGDF_DEBUG
	for (node v : V) {
		OGDF_ASSERT(v->graphOf() == &G);
	}
#endif

	// initialize the arrays
	matching.init(G, false);
	NodeArray<node> pairU(g, dummy);
	NodeArray<node> pairV(g, dummy);
	NodeArray<int> dist(g, 0);

	int size = 0;
	while (bfsMaximumCardinalityMatching(pairV, pairU, dist, dummy, gU)) {
		for (node u : gU) {
			if (pairU[u] == dummy) {
				if (dfsMaximumCardinalityMatching(pairV, pairU, dist, dummy, u)) {
					size++;
				}
			}
		}
	}

	// fill the matching-array
	for (node u : gU) {
		node buddy = pairU[u];
		if (buddy != dummy) {
			edge e = g.searchEdge(u, buddy);
			matching[g.original(e)] = true;
		}
	}

	return size;
}

}
}
