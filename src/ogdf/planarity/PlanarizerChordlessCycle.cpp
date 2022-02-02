/** \file
 * \brief Implementation of class PlanarizerChordlessCycle.
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

#include <ogdf/planarity/PlanarizerChordlessCycle.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/Queue.h>
#ifdef OGDF_DEBUG
# include <ogdf/basic/simple_graph_alg.h>
#endif
#include <set>

namespace ogdf {

PlanarizerChordlessCycle::PlanarizerChordlessCycle()
	: CrossingMinimizationModule(), m_inserter{}
	{ }

PlanarizerChordlessCycle::PlanarizerChordlessCycle(
	const PlanarizerChordlessCycle &planarizer)
	: CrossingMinimizationModule()
	{ }

CrossingMinimizationModule *PlanarizerChordlessCycle::clone() const {
	return new PlanarizerChordlessCycle(*this);
}

PlanarizerChordlessCycle &PlanarizerChordlessCycle::operator=(
	const PlanarizerChordlessCycle &planarizer)
{
	return *this;
}

bool PlanarizerChordlessCycle::findChordlessCycle(
	const Graph &G,
	List<node> &cycle)
{
	Queue<node> queue;
	NodeArray<bool> seen {G, false};
	NodeArray<edge> pred {G, nullptr};

	node src {G.firstNode()};
	seen[src] = true;
	queue.append(src);
	edge connectingEdge {nullptr};

	// Breadth-first search: As soon as an already visited node is found, the
	// predecessors of this and the current node form a chordless cycle.
	[&]{
	while (!queue.empty()) {
		node v {queue.pop()};
		for (adjEntry adj : v->adjEntries) {
			node w {adj->twinNode()};
			if (pred[v] == nullptr || w != pred[v]->opposite(v)) {
				if (seen[w]) {
					connectingEdge = adj->theEdge();
					return; // break out of lambda
				}
				if (w != v) {
					seen[w] = true;
					pred[w] = adj->theEdge();
					queue.append(w);
				}
			}
		}
	}
	}();

	// No cycle found:
	if (queue.empty()) {
		return false;
	}

	// Cycle found: Get lowest common ancestor of
	// connectingEdge->source/target() in the bfs-tree.
	node lca {nullptr};
	NodeArray<bool> visited {G, false};
	node cur {connectingEdge->source()};
	while (cur != nullptr) {
		visited[cur] = true;
		cur = pred[cur] == nullptr ? nullptr : pred[cur]->opposite(cur);
	}
	cur = connectingEdge->target();
	while (cur != nullptr) {
		if (visited[cur]) {
			lca = cur;
			break;
		}
		cur = pred[cur] == nullptr ? nullptr : pred[cur]->opposite(cur);
	}
	OGDF_ASSERT(lca != nullptr);

	// From connectingEdge, traverse predecessors of source once, then the
	// predecessors of target until the lca is reached. Collect the nodes since
	// these form the chordless cycle.
	for (bool useSource : {true, false}) {
		node currentNode {useSource ? connectingEdge->source() : connectingEdge->target()};
		while (currentNode != lca) {
			if (useSource) {
				cycle.pushBack(currentNode);
			} else {
				cycle.pushFront(currentNode);
			}
			edge predEdge {pred[currentNode]};
			currentNode = predEdge->opposite(currentNode);
		}
	}
	cycle.pushFront(lca);
	return true;
}


void PlanarizerChordlessCycle::addToGraphCopy(
	GraphCopy &graphCopy,
	GraphCopy &copyCopy,
	DynamicDualGraph &dual,
	node vOrig,
	const EdgeArray<int> *pCostOrig,
	EdgeArray<int> *pCostCopy)
{
	node vCopy {graphCopy.newNode(vOrig)};

	// Add edges in graphCopy between vCopy and other nodes in graphCopy.
	for (adjEntry adj : vOrig->adjEntries) {
		if (graphCopy.copy(adj->twinNode()) != nullptr) {
			edge eOrig {adj->theEdge()};
			edge eCopy {graphCopy.newEdge(eOrig)};
			if (pCostOrig) {
				(*pCostCopy)[eCopy] = (*pCostOrig)[eOrig];
			}
		}
	}

	// Insert the node in the planarization copyCopy.
	m_inserter.call(copyCopy, dual, vCopy, pCostCopy);
}

void PlanarizerChordlessCycle::transferToPlanRep(
	PlanRep &pr,
	const GraphCopy &graphCopy,
	const GraphCopy &copyCopy)
{
	EdgeArray<SListPure<int>> crossings {pr.original()};

	// Get number of crossing in the planarization copyCopy.
	int numCrossings {0};
	NodeArray<int> index {copyCopy, -1};
	for (node v : copyCopy.nodes) {
		if (copyCopy.isDummy(v)) {
			index[v] = numCrossings++;
		}
	}

	// For each edge in original graph: Get crossings.
	for (edge eCopyCopy : copyCopy.edges) {
		if (copyCopy.original(eCopyCopy->source()) != nullptr) {
			edge eCopy {copyCopy.original(eCopyCopy)};
			edge eOrig {graphCopy.original(eCopy)};
			ListConstIterator<edge> it = copyCopy.chain(eCopy).begin();
			for (++it; it.valid(); ++it) {
				crossings[eOrig].pushBack(index[(*it)->source()]);
			}
		}
	}

	// For all edges ePr in pr:
	// Create crossings in ePr as remembered in the crossings edge array.
	Array<node> id2Node(0, numCrossings-1, nullptr);
	SListPure<edge> edges;
	pr.allEdges(edges);

	for (edge ePr : edges) {
		edge eOrig {pr.original(ePr)};

		for (int i : crossings[eOrig]) {
			node x {id2Node[i]};
			edge ePrOld {ePr};
			ePr = pr.split(ePr);
			node y {ePr->source()};

			if (x == nullptr) {
				id2Node[i] = y;
			} else {
				pr.moveTarget(ePrOld, x);
				pr.moveSource(ePr, x);
				pr.delNode(y);
			}
		}
	}
}

Module::ReturnType PlanarizerChordlessCycle::doCall(
	PlanRep &pr,
	int cc,
	const EdgeArray<int> *pCostOrig,
	const EdgeArray<bool> *pForbiddenOrig,
	const EdgeArray<uint32_t> *pEdgeSubGraphs,
	int &crossingNumber)
{
	OGDF_ASSERT(isSimpleUndirected(pr));
	crossingNumber = 0;
	pr.initCC(cc);

	// The graph copies used here are as follows:
	// G -copy-> graphCopy (building up copy) -copy-> copyCopy (planarization)
	// G -copy-> pr (final planarization, assigned at the end)
	const Graph &G {pr.original()};
	GraphCopy graphCopy;
	graphCopy.createEmpty(G);

	// Find a chordless cycle in G. If none could be found, G is planar.
	List<node> cycle;
	bool cycleFound {findChordlessCycle(G, cycle)};
	if (!cycleFound) {
#ifdef OGDF_DEBUG
		bool planar =
#endif
			planarEmbed(pr);
		OGDF_ASSERT(planar);
		crossingNumber = 0;
		return ReturnType::Optimal;
	}

	// Start with graphCopy only containing a chordless cycle.
	NodeArray<bool> activeNodes {G, false};
	for (node cycleNodeOrig : cycle) {
		activeNodes[cycleNodeOrig] = true;
	}
	EdgeArray<edge> edgeCopies {G, nullptr};
	graphCopy.initByActiveNodes(cycle, activeNodes, edgeCopies);

	// Must hold for a chordless cycle:
	OGDF_ASSERT(graphCopy.numberOfNodes() == cycle.size());
	OGDF_ASSERT(graphCopy.numberOfEdges() == cycle.size());

	// Create array for edge costs in graphCopy.
	EdgeArray<int> *pCostCopy {nullptr};
	if (pCostOrig) {
		pCostCopy = new EdgeArray<int> {graphCopy};
		for (edge eCopy : graphCopy.edges) {
			(*pCostCopy)[eCopy] = (*pCostOrig)[graphCopy.original(eCopy)];
		}
	}

	// Initialize planarization copyCopy and corresponding dual graph.
	GraphCopy copyCopy {dynamic_cast<const Graph &>(graphCopy)};
	CombinatorialEmbedding emb {copyCopy};
	DynamicDualGraph copyCopyDual {emb};

	// For each node in the chordless cycle:
	// Start a dfs to find a node that is not yet in graphCopy.
	for (node cycleNodeOrig : cycle) {
		ArrayBuffer<node> stack;
		stack.push(cycleNodeOrig);
		while (!stack.empty()) {
			node vOrig {stack.popRet()};

			for (adjEntry adj : vOrig->adjEntries) {
				node wOrig {adj->twinNode()};

				// If a neighbor is not yet in graphCopy/copyCopy, add it to
				// them and push it on the stack.
				if (graphCopy.copy(wOrig) == nullptr) {
					addToGraphCopy(graphCopy, copyCopy, copyCopyDual,
							wOrig, pCostOrig, pCostCopy);
					stack.push(wOrig);
				}
			}
		}
	}
	delete pCostCopy;

	OGDF_ASSERT(G.numberOfNodes() == graphCopy.numberOfNodes());
	OGDF_ASSERT(G.numberOfEdges() == graphCopy.numberOfEdges());
	OGDF_ASSERT(G.numberOfNodes() == pr.numberOfNodes());
	OGDF_ASSERT(G.numberOfEdges() == pr.numberOfEdges());

	// Move the resulting planarization from copyCopy to pr.
	transferToPlanRep(pr, graphCopy, copyCopy);

	// Remove pseudo crossings and recompute crossing number.
#ifdef OGDF_DEBUG
	bool planar =
#endif
		planarEmbed(pr);
	OGDF_ASSERT(planar);
	pr.removePseudoCrossings();
	crossingNumber = computeCrossingNumber(pr, pCostOrig, pEdgeSubGraphs);

	OGDF_ASSERT(isPlanar(pr));
	OGDF_ASSERT(!pr.hasNonSimpleCrossings());
	return ReturnType::Feasible;
}

}
