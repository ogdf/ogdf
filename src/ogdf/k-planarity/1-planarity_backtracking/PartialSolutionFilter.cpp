/** \file
 * \brief Implementation of SeparatingCycleFilter methods.
 *
 * \author Matthias Pfretzschner
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

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/GraphSets.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/basic.h>
#include <ogdf/graphalg/MinSTCutMaxFlow.h>
#include <ogdf/k-planarity/1-planarity_backtracking/EdgePairPartition.h>
#include <ogdf/k-planarity/1-planarity_backtracking/OnePlanarization.h>
#include <ogdf/k-planarity/1-planarity_backtracking/PartialSolutionFilter.h>

#include <climits>
#include <deque>
#include <iterator>
#include <unordered_set>
#include <vector>

using namespace ogdf;
using namespace oneplan_backtracking;
using namespace std;

bool SeparatingCycleFilter::canCut(OnePlanarization& pl) {
	for (node cv : pl.crossingVertices()) {
		OGDF_ASSERT(cv->degree() == 4);

		// Find endpoints of the crossed edges corresponding to cv
		unordered_set<edge> crossedEdges;
		for (adjEntry adj : cv->adjEntries) {
			crossedEdges.insert(pl.original(adj->theEdge()));
		}
		OGDF_ASSERT(crossedEdges.size() == 2);
		NodePair p1 = {pl.copy((*crossedEdges.begin())->source()),
				pl.copy((*crossedEdges.begin())->target())};
		NodePair p2 = {pl.copy((*std::next(crossedEdges.begin()))->source()),
				pl.copy((*std::next(crossedEdges.begin()))->target())};

		if (test(pl, cv, p1, p2) || test(pl, cv, p2, p1)) {
			return true;
		}
	}
	return false;
}

bool SeparatingCycleFilter::test(OnePlanarization& pl, node crossingVertex, NodePair cyclePair,
		NodePair otherPair) {
	Graph::HiddenEdgeSet hes(pl);
	List<adjEntry> adjCopy;
	crossingVertex->allAdjEntries(adjCopy);
	for (adjEntry adj : adjCopy) {
		hes.hide(adj->theEdge());
	}

	EdgeSet cycleEdges;
	NodeSet cycleNodes;
	int minCycleLength =
			minFreeEdgeCycle(pl, crossingVertex, cyclePair, otherPair, cycleEdges, cycleNodes);
	// As long as there is a cycle c containing cyclePair, test c and subsequently hide it.
	while (minCycleLength < INT_MAX) {
		for (edge e : cycleEdges) {
			hes.hide(e);
		}

		GraphCopy cp(*reinterpret_cast<Graph*>(&pl));
		node v1 = cp.copy(otherPair.source);
		node v2 = cp.copy(otherPair.target);

		// Contract all edges that cannot cross any edge of the cycle
		contractUncrossableEdges(pl, cp, cycleEdges, cycleNodes, v1, v2);

		if (v1 == v2) {
			return true;
		}
		for (node n : cycleNodes) {
			cp.delNode(cp.copy(n));
		}

		// Test whether the maximum number of edge-disjoint paths between otherPair is greater than minCycleLength
		MinSTCutMaxFlow<int> minCutAlg;
		List<edge> res;
		EdgeArray<int> edgeCosts(cp, 1);
		minCutAlg.call(cp, v1, v2, res);
		if (res.size() > minCycleLength) {
			return true;
		}

		minCycleLength =
				minFreeEdgeCycle(pl, crossingVertex, cyclePair, otherPair, cycleEdges, cycleNodes);
	}
	return false;
}

void SeparatingCycleFilter::contractUncrossableEdges(const OnePlanarization& pl, GraphCopy& cp,
		const EdgeSet& cycleEdges, const NodeSet& cycleNodes, node& v1, node& v2) {
	vector<edge> contractableEdges;
	for (edge e : pl.edges) {
		bool crossableWithCycle = false;
		if (pl.freeEdges().isMember(e)) {
			if (cycleEdges.isMember(e)) {
				continue;
			}
			for (edge cycEdge : cycleEdges) {
				edge origCycEdge = pl.original(cycEdge);
				edge origE = pl.original(e);
				if (origCycEdge != nullptr
						&& pl.getEdgePairPartition()->isFree({origE, origCycEdge})) {
					crossableWithCycle = true;
					break;
				}
			}
		}
		if (!crossableWithCycle) {
			contractableEdges.push_back(e);
		}
	}

	for (edge e : contractableEdges) {
		edge eCopy = cp.copy(e);
		OGDF_ASSERT(eCopy->source()->degree() > 0 && eCopy->target()->degree() > 0);
		if (cycleNodes.isMember(e->source()) || cycleNodes.isMember(e->target())) {
			continue;
		}

		if (eCopy->target() == v1 || eCopy->target() == v2) {
			cp.reverseEdge(eCopy);
		}

		if (eCopy->source() != eCopy->target()) {
			if (eCopy->target() == v1) {
				v1 = eCopy->source();
			}
			if (eCopy->target() == v2) {
				v2 = eCopy->source();
			}
			cp.contract(eCopy, true);
		}
	}
}

int SeparatingCycleFilter::minFreeEdgeCycle(const OnePlanarization& pl, node crossingVertex,
		NodePair& cyclePair, NodePair& otherPair, EdgeSet& cycleEdges, NodeSet& cycleNodes) {
	cycleEdges.init(pl);
	cycleNodes.init(pl);
	node start = cyclePair.source;
	node target = cyclePair.target;
	deque<node> q {start};
	NodeArray<adjEntry> parent(pl, nullptr);
	NodeArray<int> dist(pl, INT_MAX - 1);
	dist[start] = 0;

	while (!q.empty()) { // 0-1 BFS
		node next = q.front();
		q.pop_front();
		if (next == target) {
			break;
		}

		for (adjEntry adj : next->adjEntries) {
			node neigh = adj->twinNode();
			OGDF_ASSERT(neigh != crossingVertex);
			if (neigh == otherPair.source || neigh == otherPair.target) {
				continue;
			}
			int edgeWeight = pl.freeEdges().isMember(adj->theEdge());
			if (dist[next] + edgeWeight < dist[neigh]) {
				dist[neigh] = dist[next] + edgeWeight;
				parent[neigh] = adj->twin();
				if (edgeWeight == 1) {
					q.push_back(neigh);
				} else {
					q.push_front(neigh);
				}
			}
		}
	}

	if (parent[target] == nullptr) {
		return INT_MAX;
	}

	OGDF_ASSERT(dist[target] < INT_MAX - 1);
	node cur = parent[target]->twinNode();
	cycleEdges.insert(parent[target]->theEdge());
	cycleNodes.insert(start);
	cycleNodes.insert(target);
	while (cur != start) {
		cycleEdges.insert(parent[cur]->theEdge());
		cycleNodes.insert(cur);
		cur = parent[cur]->twinNode();
	}

	return dist[target];
}
