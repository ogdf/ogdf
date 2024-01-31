/** \file
 * \brief Implementation of an Alternating Tree helper structure for the Blossom algorithm.
 *
 * \author Joshua Sangmeister
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
#include <ogdf/basic/List.h>
#include <ogdf/graphalg/matching_blossom/BlossomHelper.h>
#include <ogdf/graphalg/matching_blossom/Cycle.h>
#include <ogdf/graphalg/matching_blossom/Pseudonode.h>
#include <ogdf/graphalg/matching_blossom/utils.h>

#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace ogdf {
namespace matching_blossom {

template<class TWeight>
class AlternatingTree {
	//! Type of callback function when iterating the tree.
	using IteratorCallback = std::function<void(edge, bool)>;

	//! Reference to the helper class.
	BlossomHelper<TWeight>& m_helper;

	//! The root of the tree. May be empty to signalize an invalid tree which needs to be rebuild.
	node m_root;

	//! All even nodes, mapped to their respective back edge in the tree (or nullptr for the root).
	std::unordered_map<node, edge> m_evenNodes;

	//! All odd nodes, mapped to their respective back edge in the tree.
	std::unordered_map<node, edge> m_oddNodes;

	//! Epsilon test for floating point numbers.
	EpsilonTest m_eps;

	//! Find the node where the paths from both endpoints of \p e to the root first meet.
	node findCycleStartNode(edge e) {
		std::unordered_set<node> potentialIntersections = {e->source(), e->target()};
		std::vector<node> nodes = {e->source(), e->target()};
		while (true) {
			// iterate both nodes in the vector by reference and update them to the next even node
			for (node& n : nodes) {
				if (n != m_root) {
					n = getNextEvenNode(n);
					if (potentialIntersections.find(n) != potentialIntersections.end()) {
						return n;
					}
					potentialIntersections.insert(n);
				}
			}
		}
	}

	//! Return the next even node in root direction. \p v must be even. \p callback is executed for
	//! both edges on the path, if given.
	node getNextEvenNode(node v, IteratorCallback callback = nullptr) {
		OGDF_ASSERT(isEven(v));
		if (v == m_root) {
			return nullptr;
		}
		edge e = evenBackEdge(v);
		if (callback) {
			callback(e, true);
		}
		v = e->opposite(v);
		e = oddBackEdge(v);
		if (callback) {
			callback(e, false);
		}
		return e->opposite(v);
	}

	//! Iterate the tree from \p start to \p end. \p callback is executed for each edge on the path.
	//! Both \p start and \p end must be even.
	void iterateTree(node start, node end, IteratorCallback callback) {
		OGDF_ASSERT(isEven(start) && isEven(end));
		while (start != end) {
			start = getNextEvenNode(start, callback);
		}
	}

	//! Exchange the end node \p oldNode of edge \p e in graph \p graph for node \p newNode.
	void moveEdge(edge e, node oldNode, node newNode) {
		if (oldNode == e->source()) {
			m_helper.graph().moveSource(e, newNode);
		} else {
			m_helper.graph().moveTarget(e, newNode);
		}
	}

public:
	KeyIteratorContainer<node, edge> evenNodes;

	KeyIteratorContainer<node, edge> oddNodes;

	AlternatingTree(BlossomHelper<TWeight>& helper, node root = nullptr)
		: m_helper(helper), evenNodes(m_evenNodes), oddNodes(m_oddNodes) {
		reset(root);
	}

	node root() { return m_root; }

	bool hasRoot() { return m_root != nullptr; }

	edge evenBackEdge(node v) { return tryGetPointerFromMap(m_evenNodes, v); }

	bool isEven(node v) { return m_evenNodes.find(v) != m_evenNodes.end(); }

	edge oddBackEdge(node v) { return tryGetPointerFromMap(m_oddNodes, v); }

	bool isOdd(node v) { return m_oddNodes.find(v) != m_oddNodes.end(); }

	bool contains(node v) { return isEven(v) || isOdd(v); }

	//! Finds the common node between \p e and the tree. If both nodes are in the tree, only the
	//! source node is returned.
	node commonNode(edge e) {
		if (contains(e->source())) {
			return e->source();
		} else {
			OGDF_ASSERT(contains(e->target()));
			return e->target();
		}
	}

	//! Reset the tree with new the new \p root. If \p root is nullptr, the tree will be invalid.
	void reset(node root = nullptr) {
		m_root = root;
		m_evenNodes.clear();
		m_oddNodes.clear();
		if (m_root != nullptr) {
			OGDF_ASSERT(m_root->graphOf() == &m_helper.graph());
			m_evenNodes[m_root] = nullptr;
		}
	}

	/**
	 * @brief Grow the tree by adding e and f.
	 *
	 * @param u the node from which to grow, must be an even node in the tree.
	 * @param e the edge connecting u to a node incident to f.
	 * @param f a matching edge which should be added to the tree.
	 */
	void grow(node u, edge e, edge f) {
		node v = e->opposite(u);
		m_oddNodes[v] = e;
		node w = f->opposite(v);
		m_evenNodes[w] = f;
	}

	//! Return the odd-length cycle which is closed in this tree with \p cycleEdge.
	Cycle* getCycle(edge cycleEdge) {
		// Follow the path in the direction of the root node given by the tree for both
		// end nodes of the found edge until an intersection is found to form the cycle.
		Cycle* cycle = new Cycle(cycleEdge);
		node startNode = findCycleStartNode(cycleEdge);
		iterateTree(cycleEdge->source(), startNode, [&](edge e, bool isEven) { cycle->addEdge(e); });
		std::vector<edge> stack;
		iterateTree(cycleEdge->target(), startNode, [&](edge e, bool isEven) { stack.push_back(e); });
		for (auto it = stack.rbegin(); it != stack.rend(); ++it) {
			cycle->addEdge(*it);
		}
		return cycle;
	}

	//! Augment the matching along the path from \p startingEdge to the root. \p startingEdge must
	//! be incident to a free node.
	void augmentMatching(edge startingEdge) {
		m_helper.addToMatching(startingEdge);
		iterateTree(commonNode(startingEdge), m_root, [&](edge e, bool isEven) {
			if (!isEven) {
				m_helper.addToMatching(e);
			}
		});
	}

	//! Shrink the \p cycle in this tree and return the generated pseudonode.
	Pseudonode* shrink(Cycle* cycle, std::vector<std::tuple<edge, bool>>& _selfLoops) {
		node newNode = m_helper.graph().newNode();
		Pseudonode* pseudonode = new Pseudonode(newNode, cycle);
		std::unordered_map<node, edge> edgeMap; // map non-cycle nodes to their best edge
		std::unordered_map<node, std::vector<edge>> selfLoops; // map cycle nodes to all edges that will become self loops
		node matchingNode = nullptr, backNode = nullptr;
		for (node u : cycle->nodes()) {
			if (!matchingNode) {
				edge matchingEdge = m_helper.matching(u);
				if (matchingEdge != nullptr && !cycle->contains(matchingEdge->opposite(u))) {
					matchingNode = matchingEdge->opposite(u);
				}
			}
			if (!backNode) {
				edge backEdge = evenBackEdge(u);
				if (backEdge != nullptr && !cycle->contains(backEdge->opposite(u))) {
					backNode = backEdge->opposite(u);
				}
			}
			m_helper.matching(u) = nullptr;
			// since the adjacency list is modified, we have to copy it
			List<edge> edges;
			u->adjEdges(edges);
			for (edge e : edges) {
				node v = e->opposite(u);
				if (!cycle->contains(v)) {
					auto it = edgeMap.find(v);
					bool firstEdge = it == edgeMap.end();
					bool betterEdge = firstEdge
							|| m_eps.less(m_helper.getRealReducedWeight(e),
									m_helper.getRealReducedWeight(it->second));
					if (!firstEdge) {
						edge selfLoop = betterEdge ? it->second : e;
						node loopNode = selfLoop->opposite(v);
						selfLoops[loopNode].push_back(selfLoop);
					}
					if (betterEdge) {
						edgeMap[v] = e;
					}
				}
			}
			if (u == m_root) {
				m_root = newNode;
			}
		}
		// move self loops
		std::unordered_map<node, edge> selfLoopMap;
		for (node u : cycle->nodes()) {
			bool even = isEven(u);
			for (edge e : selfLoops[u]) {
				node v = e->opposite(u);
				edge bestEdge = edgeMap[v];
				node w = bestEdge->opposite(v);
				// set edge weight to difference between reduced weights of this edge to used edge
				// in edge map
				// this means the reduced weight of the edge becomes invalid, but self loops are
				// simply ignored in all calculations
				m_helper.c(e) -= m_helper.c(bestEdge) + m_helper.y(u) - m_helper.y(w);
				pseudonode->addReference(bestEdge, e, m_helper.pseudonode(v));
				moveEdge(e, v, u);
				_selfLoops.push_back(std::make_tuple(e, even));
				if (evenBackEdge(v) == e) {
					m_evenNodes[v] = bestEdge;
				} else if (oddBackEdge(v) == e) {
					m_oddNodes[v] = bestEdge;
				}
			}
			if (even) {
				m_evenNodes.erase(u);
			} else {
				m_oddNodes.erase(u);
			}
		}

		for (auto entry : edgeMap) {
			node v = entry.first;
			edge e = entry.second;
			node u = e->opposite(v);
			// edge between a cycle node and a non-cycle node: bend it to new node
			moveEdge(e, u, newNode);
			m_helper.c(e) -= m_helper.y(u);
		}
		m_evenNodes[newNode] = backNode ? edgeMap[backNode] : nullptr;
		if (matchingNode) {
			m_helper.addToMatching(edgeMap[matchingNode]);
		}
		m_helper.addPseudonode(pseudonode);
		return pseudonode;
	}

	//! Expand the given \p pseudonode in this tree.
	void expand(Pseudonode* pseudonode) {
		node graphNode = pseudonode->graphNode;
		auto cycle = pseudonode->cycle;

		// move edges back to old nodes
		List<edge> refEdges;
		graphNode->adjEdges(refEdges);
		for (edge e : refEdges) {
			node uOrig = m_helper.getBaseNode(e, graphNode);
			node u = m_helper.reprChild(uOrig);
			OGDF_ASSERT(m_helper.repr(u) == graphNode);
			m_helper.c(e) += m_helper.y(u);
			moveEdge(e, graphNode, u);
		}
		// move back self loops
		while (!refEdges.empty()) {
			edge refEdge = refEdges.popFrontRet();
			for (edge e : pseudonode->referenceEdges.selfLoops(refEdge)) {
				refEdges.pushBack(e);
				node source, target;
				std::tie(source, target) = m_helper.getBaseNodes(e);
				node currentSource = m_helper.repr(source);
				node u, v;
				if (currentSource == graphNode) {
					u = m_helper.reprChild(source);
					v = m_helper.repr(target);
					m_helper.graph().moveSource(e, u);
					m_helper.graph().moveTarget(e, v);
				} else {
					u = m_helper.reprChild(target);
					v = currentSource;
					m_helper.graph().moveSource(e, v);
					m_helper.graph().moveTarget(e, u);
				}
				OGDF_ASSERT(m_helper.repr(u) == graphNode);
				node w = refEdge->opposite(v);
				m_helper.c(e) += m_helper.c(refEdge) + m_helper.y(u) - m_helper.y(w);
				pseudonode->referenceEdges.removeFromOther(e);
			}
		}
#ifdef OGDF_DEBUG
		// Check that all self loops have been moved
		for (node u : cycle->nodes()) {
			List<edge> inEdges;
			u->inEdges(inEdges);
			for (edge e : inEdges) {
				OGDF_ASSERT(!e->isSelfLoop());
			}
		}
#endif

		// find the matching edges which go in and out of the pseudonode in the current tree and
		// mark the respective nodes in the cycle as start and end node
		// the tree enters the pseudonode/cycle at startCycleNode with a non-matching edge since
		// the pseudonode was odd and leaves it at endCycleNode with a matching edge
		node startCycleNode = cycle->startNode();
		auto it = m_oddNodes.find(graphNode);
		if (it != m_oddNodes.end()) {
			edge e = it->second;
			node origStart = m_helper.getBaseNode(e, graphNode);
			startCycleNode = m_helper.reprChild(origStart);
			OGDF_ASSERT(cycle->contains(startCycleNode));
			m_oddNodes[startCycleNode] = e;
		}
		edge matchingEdge = m_helper.matching(graphNode);
		node origEnd = m_helper.getBaseNode(matchingEdge, graphNode);
		node endCycleNode = m_helper.reprChild(origEnd);
		OGDF_ASSERT(cycle->contains(endCycleNode));
		m_helper.addToMatching(matchingEdge);

		size_t startNodeEdgeIndex, endNodeEdgeIndex;
		std::tie(startNodeEdgeIndex, endNodeEdgeIndex) = cycle->indexOf(startCycleNode, endCycleNode);

		// We start at the endCycleNode and iterate the even-length path to the startCycleNode and
		// then back via the odd-length path and add every second edge to the matching.
		// Simultanously, we add the nodes and edges of the even length path to the tree.
		// the direction of the iteration depends the direction of the cycle and the order of the
		// start/end nodes therein.
		auto edges = cycle->edgeOrder();
		bool moveForward = (endNodeEdgeIndex - startNodeEdgeIndex) % 2
				== (startNodeEdgeIndex < endNodeEdgeIndex);
		node currentNode = endCycleNode;
		for (size_t i = 0; i < edges.size() - 1; ++i) {
			int index = endNodeEdgeIndex;
			if (moveForward) {
				index += i + 1;
			} else {
				index += edges.size() - i;
			}
			edge e = edges[index % edges.size()];
			if (i % 2 == 1) {
				m_helper.addToMatching(e);
			}
			if (currentNode != startCycleNode) {
				if (i % 2 == 0) {
					m_oddNodes[currentNode] = e;
				} else {
					m_evenNodes[currentNode] = e;
				}
				currentNode = e->opposite(currentNode);
			}
		}

		// remove pseudonode from all data structures
		m_helper.removePseudonode(pseudonode);
		m_oddNodes.erase(graphNode);
	}
};

}
}
