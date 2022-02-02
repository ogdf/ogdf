/** \file
 * \brief Implementation of class StarInserter.
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

#include <ogdf/planarity/StarInserter.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/graphalg/Dijkstra.h>
#include <set>

namespace ogdf {

StarInserter &StarInserter::operator=(const StarInserter &inserter)
{
	return *this;
}

node StarInserter::getOptimalDualNode(node origNode,
		const EdgeArray<int> *pCostOrig,
		PredecessorMap &predecessors)
{
	const Graph &dualGraph {m_dual->getGraph()};

	// Setup parameters for shortest path computations.
	EdgeArray<int> weight {dualGraph}; // weights of edges in dualGraph
	for (edge e : dualGraph.edges) {
		weight[e] = pCostOrig == nullptr ? 1 :
			(*pCostOrig)[m_graphCopy->original(m_dual->primalEdge(e))];
	}
	NodeArray<int> distance {dualGraph}; // tmp result for each Dijkstra run
	NodeArray<int> distanceSum {dualGraph, 0}; // sums over all Dijkstra runs
	Dijkstra<int> dijkstra;

	// For every (copy) neighbor w of origNode:
	for (adjEntry origAdj : origNode->adjEntries) {
		edge origEdge {origAdj->theEdge()};
		node w {m_graphCopy->copy(origEdge->opposite(origNode))};

		// Collect dual nodes of faces incident to w as sources.
		std::set<node> sourceSet;
		for (adjEntry adj : w->adjEntries) {
			sourceSet.insert(m_dual->dualNode(m_combEmbedding->rightFace(adj)));
		}

		List<node> sources;
		for (node v : sourceSet) {
			sources.pushBack(v);
		}

		// Get shortest paths starting at collected sources.
		predecessors[w].reset(new NodeArray<edge>{dualGraph});
		dijkstra.call(dualGraph, weight, sources,
			*(predecessors[w]), distance, false);

		// Update the distance sums for every dual node.
		for (auto dualNode : dualGraph.nodes) {
			distanceSum[dualNode] += distance[dualNode] *
				(pCostOrig == nullptr ? 1 : (*pCostOrig)[origEdge]);
		}
	}

	// Get the optimal dual node, i.e. the one with the minimum distance over
	// all shortest path searches.
	node optimalDualNode {dualGraph.firstNode()};
	int minimumDistanceSum {std::numeric_limits<int>::max()};
	for (auto dualNode : dualGraph.nodes) {
		if (distanceSum[dualNode] < minimumDistanceSum) {
			optimalDualNode = dualNode;
			minimumDistanceSum = distanceSum[dualNode];
		}
	}
	return optimalDualNode;
}

void StarInserter::makePredsConsistent(node origNode,
		node optimalDualNode,
		PredecessorMap &predecessors)
{
	std::unordered_map<node, node> dualNodeNeighbor;
	NodeArray<edge> successor {*m_dual, nullptr};

	// Create a single shortest path tree in the form of a successor array.
	for (adjEntry origAdj : origNode->adjEntries) {
		edge origEdge {origAdj->theEdge()};
		node w {m_graphCopy->copy(origEdge->opposite(origNode))};

		edge dualPredEdge {(*predecessors[w])[optimalDualNode]};
		node lastDualNode {optimalDualNode};
		while (dualPredEdge != nullptr) {
			node opposite {dualPredEdge->opposite(lastDualNode)};
			if (successor[opposite] == nullptr) {
				// Only set successor if there does not exist one already
				// because of a different predecessor path.
				successor[opposite] = dualPredEdge;
			}
			lastDualNode = opposite;
			dualPredEdge = (*predecessors[w])[lastDualNode];
		}

		dualNodeNeighbor[w] = lastDualNode;
	}

	// Update predecessor arrays to conform to the shortest path tree given by
	// successor: Start at the dual node neighbor of w, follow successors up to
	// the inserted node (the copy of origNode).
	for (adjEntry origAdj : origNode->adjEntries) {
		edge origEdge {origAdj->theEdge()};
		node w {m_graphCopy->copy(origEdge->opposite(origNode))};

		node lastDualNode {dualNodeNeighbor[w]};
		edge dualSuccEdge {successor[lastDualNode]};
		while (dualSuccEdge != nullptr) {
			lastDualNode = dualSuccEdge->opposite(lastDualNode);
			(*predecessors[w])[lastDualNode] = dualSuccEdge;
			dualSuccEdge = successor[lastDualNode];
		}
	}
}

adjEntry StarInserter::getCrossedAdjEntry(edge primalEdgeToSplit, node leftDualNode)
{
	OGDF_ASSERT((*m_edgeInChainToSplit)[primalEdgeToSplit] == nullptr ||
		(*m_edgeInChainToSplit)[primalEdgeToSplit] == primalEdgeToSplit);
	face leftFace {m_dual->primalFace(leftDualNode)};

	adjEntry adjSrc {primalEdgeToSplit->adjSource()};
	adjEntry adjTgt {primalEdgeToSplit->adjTarget()};
	if (m_combEmbedding->leftFace(adjSrc) == leftFace) {
		return adjSrc;
	} else {
		OGDF_ASSERT(m_combEmbedding->leftFace(adjTgt) == leftFace);
		return adjTgt;
	}
}

adjEntry StarInserter::getAdjEntry(node primalNode,
		node rightDualNode,
		node otherPrimalNode)
{
	if (primalNode->degree() <= 1) {
		// If there only exists one adjEntry, return it.
		// Return nullptr if no adjEntry exists.
		return primalNode->firstAdj();
	}

	// If otherPrimalNode is isolated, return the adjEntry whose right face
	// corresponds to rightDualNode.
	face rightFace {oldPrimalFace(rightDualNode)};
	if (otherPrimalNode->degree() == 0) {
		for (adjEntry adj : primalNode->adjEntries) {
			if ((*m_newToOldFace)[m_combEmbedding->rightFace(adj)] == rightFace) {
				return adj;
			}
		}
		OGDF_ASSERT(false);
	}

	// Otherwise return the common face of primalNode and otherPrimalNode.
	adjEntry adj {
		m_combEmbedding->findCommonFace(primalNode, otherPrimalNode, false)
	};
	OGDF_ASSERT((*m_newToOldFace)[m_combEmbedding->rightFace(adj)] == rightFace);
	return adj;
}

adjEntry StarInserter::getAdjEntry(node primalNode, node rightDualNode, edge primalEdge, bool first)
{
	OGDF_ASSERT((*m_edgeInChainToSplit)[primalEdge] == nullptr ||
		(*m_edgeInChainToSplit)[primalEdge] == primalEdge);
	if (primalNode->degree() <= 1) {
		// If there only exists one adjEntry, return it.
		// Return nullptr if no adjEntry exists.
		return primalNode->firstAdj();
	}

	adjEntry adjSrc {primalEdge->adjSource()};
	adjEntry adjTgt {primalEdge->adjTarget()};
	face rightFace {oldPrimalFace(rightDualNode)};
	adjEntry adjStart {nullptr};
	if ((*m_newToOldFace)[m_combEmbedding->rightFace(adjSrc)] == rightFace) {
		adjStart = adjSrc;
	} else {
		OGDF_ASSERT((*m_newToOldFace)[m_combEmbedding->rightFace(adjTgt)] == rightFace);
		adjStart = adjTgt;
	}
	// We should never cross through an edge that is incident to an endpoint of
	// the inserted edge.
	OGDF_ASSERT(adjStart->theNode() != primalNode);

	// Cycle through adjEntries of the face in clockwise order, starting at adj
	// of primalEdge. Return the first adj at the correct primal node.
	adjEntry adjNext {first ? adjStart->clockwiseFacePred() :
		adjStart->clockwiseFaceSucc()};
	while (adjNext != adjStart) {
		if (adjNext->theNode() == primalNode) {
			return adjNext;
		}
		adjNext = first ? adjNext->clockwiseFacePred() :
			adjNext->clockwiseFaceSucc();
	}

	// primalNode and primalEdge have no common face. This should not happen.
	OGDF_ASSERT(false);
	return nullptr;
}

edge StarInserter::collectAdjEntries(node w,
		node insertedNode,
		node optimalDualNode,
		const PredecessorMap &predecessors,
		List<adjEntry> &crossedEdges)
{
	// Follow predecessors of optimal dual node in the insertion path and
	// collect crossed edges for use in GraphCopy::insertEdgePathEmbedded().
	edge dualPredEdge {(*predecessors.at(w))[optimalDualNode]};
	node lastDualNode {optimalDualNode};

	// First, the adjEntry of insertedNode with optimalDualNode to its right.
	edge edgeToSplit {nullptr};
	if (dualPredEdge) {
		// In case the edge has been split already by some insertion path,
		// determine the correct edge in the chain to split.
		edgeToSplit = (*m_edgeInChainToSplit)[
			(*m_originalEdge)[m_dual->primalEdge(dualPredEdge)]
		];
		dualPredEdge = m_dual->dualEdge(edgeToSplit);
		crossedEdges.pushBack(getAdjEntry(insertedNode, lastDualNode, edgeToSplit, true));

		// For the case that optimalDualNode was already split previously, make
		// sure that lastDualNode is actually incident to dualPredEdge.
		if (oldPrimalFace(dualPredEdge->source()) == oldPrimalFace(lastDualNode)) {
			lastDualNode = dualPredEdge->source();
		} else {
			OGDF_ASSERT(oldPrimalFace(dualPredEdge->target())
				== oldPrimalFace(lastDualNode));
			lastDualNode = dualPredEdge->target();
		}
	} else {
		crossedEdges.pushBack(getAdjEntry(insertedNode, optimalDualNode, w));
	}

	// Then, the adjEntries of crossed edges (crossed from left to right).
	edge lastDualPredEdge {dualPredEdge};
	while (dualPredEdge != nullptr) {
		crossedEdges.pushBack(getCrossedAdjEntry(edgeToSplit, lastDualNode));

		// Get the next dual node.
		if (oldPrimalFace(dualPredEdge->source()) == oldPrimalFace(lastDualNode)) {
			lastDualNode = dualPredEdge->target();
		} else {
			OGDF_ASSERT(oldPrimalFace(dualPredEdge->target())
				== oldPrimalFace(lastDualNode));
			lastDualNode = dualPredEdge->source();
		}

		lastDualPredEdge = dualPredEdge;

		// Get the next predecessor edge. In case the edge has been split
		// already by some insertion path, determine the correct edge in the
		// chain to split.
		dualPredEdge = (*predecessors.at(w))[
			m_dual->dualNode(oldPrimalFace(lastDualNode))
		];
		if (dualPredEdge) {
			edgeToSplit = (*m_edgeInChainToSplit)[
				(*m_originalEdge)[m_dual->primalEdge(dualPredEdge)]
			];
			dualPredEdge = m_dual->dualEdge(edgeToSplit);
		}
	}

	// Lastly, the adjEntry of w with lastDualNode to its right.
	crossedEdges.pushBack(
		lastDualPredEdge ?
		getAdjEntry(w, lastDualNode, m_dual->primalEdge(lastDualPredEdge), false) :
		getAdjEntry(w, lastDualNode, insertedNode)
	);

	// If the first adjEntry is nullptr, we are in the first iteration of the
	// main loop and the inserted node has no incident edges yet (and hence no
	// adjEntries). Create a dummy edge such that the first edge can be embedded
	// easily (delete the dummy edge at the end of the main loop).
	edge dummyEdge {nullptr};
	ListConstIterator<adjEntry> itStart {crossedEdges.begin()};
	if (*itStart == nullptr) {
		++itStart;
		dummyEdge = m_dual->addEdgeToIsolatedNodePrimal(insertedNode,
			lastDualPredEdge ? (*itStart)->twin() : *itStart);
		crossedEdges.popFront();
		crossedEdges.pushFront(insertedNode->firstAdj());
	}
	return dummyEdge;
}

void StarInserter::transferCrossedEdges(
	const List<adjEntry> &crossedEdges,
	SList<adjEntry> &finalCrossedEdges,
	bool startAtSource)
{
	// GraphCopy::insertEdgePathEmbedded() expects an SList, not a List.
	// If the inserted edge is directed from insertedNode to w, just copy the
	// list elements.
	if (startAtSource) {
		for (adjEntry adj : crossedEdges) {
			finalCrossedEdges.pushBack(adj);
		}
	} else {
		// If the edge is directed from w to insertedNode, the list of
		// adjEntries has to be reversed for insertEdgePathEmbedded. All
		// adjEntries (except the first and last) have to be reversed, too.
		ListConstReverseIterator<adjEntry> it {crossedEdges.rbegin()};
		finalCrossedEdges.pushBack(*it);
		for (it++; it.valid() && it.succ().valid(); it++) {
			finalCrossedEdges.pushBack((*it)->twin());
		}
		finalCrossedEdges.pushBack(*it);
	}

#ifdef OGDF_HEAVY_DEBUG
	// The final adjEntries can be used for insertEdgePathEmbedded().
	SListConstIterator<adjEntry> it {finalCrossedEdges.begin()};
	for (; it.valid() && it.succ().valid() && it.succ().succ().valid(); it++) {
		OGDF_ASSERT(m_combEmbedding->rightFace(*it) ==
		            m_combEmbedding->leftFace(*(it.succ())));
	}
	OGDF_ASSERT(m_combEmbedding->rightFace(*it) ==
				m_combEmbedding->rightFace(*(it.succ())));
#endif
}

void StarInserter::initMemberData(GraphCopy &graphCopy, DynamicDualGraph &dualGraph)
{
	m_graphCopy = &graphCopy;
	m_combEmbedding = &dualGraph.getPrimalEmbedding();
	m_dual = &dualGraph;

	// Remember which faces existed before star insertion paths were inserted.
	m_newToOldFace = new FaceArray<face>{*m_combEmbedding, nullptr};
	for (face f : m_combEmbedding->faces) {
		(*m_newToOldFace)[f] = f;
	}

	// Remember the edges in the GraphCopy before insertion paths were inserted.
	// For each such edge, specify the part of its chain that should be split.
	m_edgeInChainToSplit = new EdgeArray<edge>{*m_combEmbedding, nullptr};
	m_originalEdge = new EdgeArray<edge>{*m_combEmbedding, nullptr};
	for (edge e : m_combEmbedding->getGraph().edges) {
		(*m_edgeInChainToSplit)[e] = e;
		(*m_originalEdge)[e] = e;
	}

}

void StarInserter::updateMemberData(edge origEdge, bool startAtSource)
{
	// Traverse the chain of the newly inserted edge.
	edge prevCopyEdge {nullptr};
	for (edge copyEdge : m_graphCopy->chain(origEdge)) {
		// For each new face, remember the corresponding old face prior to the
		// edge insertion. This is needed since predecessors only refers to
		// these old faces.
		face newFace {m_combEmbedding->leftFace(copyEdge->adjSource())};
		face oldFace {m_combEmbedding->rightFace(copyEdge->adjSource())};
		if ((*m_newToOldFace)[oldFace] == nullptr) {
			std::swap(oldFace, newFace);
		}
		(*m_newToOldFace)[newFace] = (*m_newToOldFace)[oldFace];

		// For each split edge, remember the part of its chain that can be split
		// next (its the next one in a clockwise order such that the later
		// inserted edges do not cross old inserted ones).
		if (prevCopyEdge != nullptr) {
			bool nextIsEdgeToSplit {true};
			edge newEdge {nullptr};
			edge oldEdge {nullptr};

			// Iterate over the incident edges e of the crossing:
			node crossing {copyEdge->commonNode(prevCopyEdge)};
			edge edgeToSplit {crossing->firstAdj()->theEdge()};
			for (adjEntry adj : crossing->adjEntries) {
				edge e {adj->theEdge()};

				// If e is not part of the inserted edge:
				if (e != copyEdge && e != prevCopyEdge) {
					if ((*m_originalEdge)[e] == nullptr) {
						// e is the new edge created by a split.
						newEdge = e;
					} else {
						// e is the old edge involved in a split.
						oldEdge = e;
					}
					if (nextIsEdgeToSplit) {
						// e is the next edge to split.
						edgeToSplit = e;
					}
					nextIsEdgeToSplit = false;
				} else {
					// If copyEdge is closer to w than prevCopyEdge,
					// use the edge after copyEdge as the one to split.
					// Otherwise, use the edge after prevCopyEdge.
					nextIsEdgeToSplit = (e == copyEdge) == startAtSource;
				}
			}

			edge originalEdge {(*m_originalEdge)[oldEdge]};
			(*m_originalEdge)[newEdge] = originalEdge;
			(*m_edgeInChainToSplit)[originalEdge] = edgeToSplit;
		}
		prevCopyEdge = copyEdge;
	}
}

void StarInserter::call(
	GraphCopy &graphCopy,
	DynamicDualGraph &dualGraph,
	node origNode,
	const EdgeArray<int> *pCostOrig)
{
	OGDF_ASSERT(graphCopy.representsCombEmbedding());

	initMemberData(graphCopy, dualGraph);

	// Compute the dual node that is optimal for the insertion of a star.
	// Also get the insertion paths in form of a predecessor relation.
	// Since they may cross each other, ensure that the entirety of the
	// insertion paths forms a tree with optimalDualNode as the root.
	PredecessorMap predecessors;
	node optimalDualNode {getOptimalDualNode(origNode, pCostOrig, predecessors)};
	makePredsConsistent(origNode, optimalDualNode, predecessors);

	// If multiple edges e1...ek cross the same edge a, we have to determine the
	// order of these crossings. Otherwise e1...ek might cross each other.
	ArrayBuffer<edge> sortedEdges;
	for (adjEntry adj : origNode->adjEntries) {
		sortedEdges.push(adj->theEdge());
	}
	EdgeOrderComparer comp(origNode, optimalDualNode, predecessors,
		*m_graphCopy, m_dual);
	sortedEdges.quicksort(comp);

	// Insert origNode (in optimal face given by optimalDualNode).
	node insertedNode {m_graphCopy->newNode(origNode)};

	// Determine optimal insertion paths from the optimal dual node to the dual
	// nodes incident to insertedNode's neighbors.
	// For every copy neighbor w of origNode:
	edge dummyEdge {nullptr};
	for (edge origEdge : sortedEdges) {
		bool startAtSource {origEdge->source() == origNode};
		node w {m_graphCopy->copy(origEdge->opposite(origNode))};

		// Collect the adjEntries for all crossed adjEntries in the insertion
		// path. Transfer them to an SList and pass them to
		// GraphCopy::insertEdgePathEmbedded() to embed origEdge.
		List<adjEntry> crossedEdges;
		edge createdEdge {collectAdjEntries(w, insertedNode, optimalDualNode,
			predecessors, crossedEdges)};
		if (dummyEdge == nullptr && createdEdge != nullptr) {
			dummyEdge = createdEdge;
		}

		SList<adjEntry> finalCrossedEdges;
		transferCrossedEdges(crossedEdges, finalCrossedEdges, startAtSource);
		m_graphCopy->insertEdgePathEmbedded(origEdge, *m_combEmbedding, *m_dual,
			finalCrossedEdges);
		updateMemberData(origEdge, startAtSource);
	}

	// Delete the dummy edge that was created.
	OGDF_ASSERT(dummyEdge != nullptr);
	m_dual->joinFacesPrimal(dummyEdge);

	// Remove all non-simple crossings.
	graphCopy.removeNonSimpleCrossings(origNode, &dualGraph);

	// Free memory.
	delete m_newToOldFace;
	delete m_edgeInChainToSplit;
	delete m_originalEdge;

	// Verify computed planarization and free memory.
	OGDF_ASSERT(m_graphCopy->representsCombEmbedding());
}

}
