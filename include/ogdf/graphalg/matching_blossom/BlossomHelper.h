/** \file
 * \brief Utility class for the Blossom algorithm, providiing
 * access to all important data structures.
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

#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/EpsilonTest.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/Logger.h>
#include <ogdf/basic/NodeArray.h>
#include <ogdf/graphalg/matching_blossom/Pseudonode.h>
#include <ogdf/graphalg/matching_blossom/utils.h>

#include <algorithm>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <unordered_set>

namespace ogdf {
namespace matching_blossom {

//! Helper class for the blossom matching algorithms.
template<class TWeight>
class BlossomHelper : private Logger {
protected:
	//! Whether or not to use the greedy initialization.
	bool m_greedyInit;

	//! A copy of the graph to work on
	GraphCopySimple m_graph;

	//! LP-induced data used by the algorithm
	EdgeArray<TWeight> m_c;
	NodeArray<TWeight> m_y;

	//! The current matching, mapping both endpoints to the corresponding matching edge
	NodeArray<edge> m_matching;

	//! A mapping of all pseudonodes in the graph
	std::unordered_map<node, Pseudonode*> m_pseudonodes;

	//! The tree induced by the pseudonodes, mapping each node index to its parent
	std::vector<node> m_repr;

	//! A shortcut array for the tree, mapping each node index to the penultimate node
	std::vector<node> m_shortcuts;

	//! The epsilon test used for floating point comparisons
	EpsilonTest m_eps;

	// Multiply all integer weights by 4 to ensure that all primal and dual values stay integral.
	static const int WEIGHT_FACTOR = std::numeric_limits<TWeight>::is_integer ? 4 : 1;

	//! Initializes the dual solution with the given y values.
	bool initDualSolution(NodeArray<TWeight>& minY) {
		// for _some_ reason, the array has to be initialized with 0, otherwise the algorithm fails
		m_y.init(m_graph, 0);
		for (node v : m_graph.nodes) {
			if (v->degree() == 0) {
				return false;
			}
			m_y[v] = minY[v];
		}
		if (!m_greedyInit) {
			return true;
		}
		for (node v : m_graph.nodes) {
			if (matching(v) != nullptr) {
				continue;
			}
			edge matchingEdge = nullptr;
			TWeight maxChange = infinity<TWeight>();
			for (auto adj : v->adjEntries) {
				edge e = adj->theEdge();
				node w = adj->twinNode();
				TWeight reducedWeight = getReducedWeight(e);
				if (m_eps.leq(reducedWeight, maxChange)) {
					maxChange = reducedWeight;
					if (!matching(w)) {
						matchingEdge = e;
					}
				}
			}
			y(v) += maxChange;
			if (matchingEdge != nullptr && isEqualityEdge(matchingEdge)) {
				addToMatching(matchingEdge);
			}
		}
#ifdef OGDF_DEBUG
		for (node v : m_graph.nodes) {
			lout() << "initial y for node " << v << ": " << y(v) << std::endl;
			edge e = matching(v);
			if (e && e->source() == v) {
				OGDF_ASSERT(e == matching(e->target()));
				OGDF_ASSERT(isEqualityEdge(e));
				lout() << "initial matching edge: " << e << std::endl;
			}
		}
		for (edge e : m_graph.edges) {
			OGDF_ASSERT(getReducedWeight(e) >= 0);
			if (isEqualityEdge(e)) {
				OGDF_ASSERT(matching(e->source()) || matching(e->target()));
			}
		}
#endif
		return true;
	}

	//! Finds the parent of \p v in the tree induced by the pseudonodes.
	node findParentInRepr(node v, node child = nullptr) {
		node parent = m_repr[v->index()];
		OGDF_ASSERT(parent != nullptr);
		if (parent == v) {
			return child;
		} else {
			node grandparent = findParentInRepr(parent, v);
			m_shortcuts[v->index()] = grandparent;
			return grandparent;
		}
	}

	//! Deletes all pseudonodes.
	void deletePseudonodes() {
		for (auto entry : m_pseudonodes) {
			delete entry.second;
		}
		m_pseudonodes.clear();
	}

public:
#ifdef OGDF_HEAVY_DEBUG
	void getReducedEdgeWeights(std::unordered_map<edge, TWeight>& edges) {
		for (edge e : m_graph.edges) {
			if (repr(e->source()) == e->source() && repr(e->target()) == e->target()) {
				edges[e] = getRealReducedWeight(e);
			}
		}
	}

	void assertConsistentEdgeWeights(std::unordered_map<edge, TWeight>& edges) {
		for (edge e : m_graph.edges) {
			if (repr(e->source()) == e->source() && repr(e->target()) == e->target()
					&& edges.find(e) != edges.end()) {
				OGDF_ASSERT(edges[e] == getRealReducedWeight(e));
			}
		}
	}
#endif
	/**
	 * @brief Construct a new Blossom V Helper object
	 *
	 * @param greedyInit whether or not to use the greedy initialization
	 */
	BlossomHelper(bool greedyInit) : m_greedyInit(greedyInit) { }

	~BlossomHelper() { deletePseudonodes(); }

	//! Reinitialize the helper class with a new graph and edge weights. Resets all helper members.
	//! Returns false if the graph cannot have a perfect matching.
	template<class WeightContainer>
	bool init(const Graph& graph, const WeightContainer& weights) {
		if (graph.numberOfNodes() % 2 == 1) {
			return false;
		}
		m_graph.clear();
		m_graph.init(graph);
		m_c.init(m_graph);
		NodeArray<TWeight> minY(m_graph, -1);
		for (edge e : m_graph.edges) {
			if (e->isSelfLoop()) {
				throw std::runtime_error("Self-loops are not supported.");
			}
			// copy initial edge costs from original edge, multiplied by WEIGHT_FACTOR
			TWeight weight = WEIGHT_FACTOR * getWeight<TWeight>(m_graph.original(e), weights);
			m_c[e] = weight;
			// init dual variables for all nodes
			auto halfWeight = weight / 2;
			for (node v : e->nodes()) {
				if (minY[v] == -1 || m_eps.less(halfWeight, minY[v])) {
					minY[v] = halfWeight;
				}
			}
		}
		m_matching.init(m_graph, nullptr);
		deletePseudonodes();
		m_repr.clear();
		m_shortcuts.clear();
		for (node v : m_graph.nodes) {
			m_repr.push_back(v);
			m_shortcuts.push_back(v);
		}
		return initDualSolution(minY);
	}

	/* Getters */

	GraphCopySimple& graph() { return m_graph; }

	//! Returns the base edge weight of \p e.
	TWeight& c(edge e) { return m_c[e]; }

	//! Returns the base y value of \p v.
	TWeight& y(node v) { return m_y[v]; }

	NodeArray<edge>& matching() { return m_matching; }

	//! Returns the matching edge of \p v, or nullptr if \p v is not matched.
	edge& matching(node v) { return m_matching[v]; }

	std::unordered_map<node, Pseudonode*>& pseudonodes() { return m_pseudonodes; }

	//! Returns the pseudonode corresponding to \p v, or nullptr if \p v is not a pseudonode.
	Pseudonode* pseudonode(node v) { return tryGetPointerFromMap(m_pseudonodes, v); }

	//! Checks whether \p v is a pseudonode.
	bool isPseudonode(node v) { return pseudonode(v) != nullptr; }

	/* End of getters */

	//! Returns the representative of \p v in the tree induced by the pseudonodes.
	node repr(node v) {
		node parent = m_repr[v->index()];
		if (parent == v || parent == nullptr) {
			return v;
		}
		node child = reprChild(v);
		parent = m_repr[child->index()];
		if (parent == nullptr) {
			return child;
		} else {
			return parent;
		}
	}

	//! Returns the child of the representative of \p v in the tree induced by the pseudonodes.
	node reprChild(node v) {
		node shortcut = m_shortcuts[v->index()];
		node parent = m_repr[shortcut->index()];
		if (parent == shortcut || parent == nullptr) {
			return findParentInRepr(v);
		} else {
			parent = findParentInRepr(shortcut);
			m_shortcuts[v->index()] = parent;
			return parent;
		}
	}

	//! Adds a pseudonode to the current representation structure.
	void addPseudonode(Pseudonode* pseudonode) {
		m_pseudonodes[pseudonode->graphNode] = pseudonode;
		m_repr.push_back(pseudonode->graphNode);
		m_shortcuts.push_back(pseudonode->graphNode);
		OGDF_ASSERT(m_repr[pseudonode->graphNode->index()] == pseudonode->graphNode);
		for (node v : pseudonode->cycle->nodes()) {
			m_repr[v->index()] = pseudonode->graphNode;
		}
	}

	//! Removes a pseudonode from the current representation structure.
	void removePseudonode(Pseudonode* pseudonode) {
		node theNode = pseudonode->graphNode;
		m_pseudonodes.erase(theNode);
		m_matching[theNode] = nullptr;
		m_repr[theNode->index()] = nullptr;
		m_graph.delNode(theNode);
		for (node v : pseudonode->cycle->nodes()) {
			m_repr[v->index()] = v;
		}
	}

	//! Expands a pseudonode in the current representation structure.
	void expandRepr(Pseudonode* pseudonode) {
		m_repr[pseudonode->graphNode->index()] = nullptr;
		for (node v : pseudonode->cycle->nodes()) {
			m_repr[v->index()] = v;
		}
	}

	//! Returns the original end point of \p e which is currently represented by \p v.
	node getBaseNode(edge e, node v) { return getBaseNodes(e, v).first; }

	//! Returns the original end point of \p e where the other end point is currently represented by \p v.
	node getOppositeBaseNode(edge e, node v) { return getBaseNodes(e, v).second; }

	//! Return both original end points of \p e where the first end point is currently represented by \p v.
	std::pair<node, node> getBaseNodes(edge e, node v = nullptr) {
		edge eOrig = m_graph.original(e);
		node source = m_graph.copy(eOrig->source());
		node target = m_graph.copy(eOrig->target());
		if (v != nullptr && repr(target) == v) {
			return {target, source};
		} else {
			OGDF_ASSERT(v == nullptr || repr(source) == v);
			return {source, target};
		}
	}

	//! Helper function to determine whether a floating point value is 0.
	bool isZero(TWeight x) { return m_eps.equal(x, (TWeight)0); }

	//! Returns the reduced weight of \p e, taking into account the y values of the endpoints.
	TWeight getReducedWeight(edge e) { return c(e) - y(e->source()) - y(e->target()); }

	//! Returns the real reduced weight. Can be overridden by subclasses to consider additional factors.
	virtual TWeight getRealReducedWeight(edge e) { return getReducedWeight(e); }

	//! Checks whether \p e is an equality edge, i.e. the reduced weight is 0.
	virtual bool isEqualityEdge(edge e) { return isZero(getRealReducedWeight(e)); }

	//! Checks whether \p v is a zero-cost node, i.e. the y value is 0.
	bool isZeroCostNode(node v) { return isZero(y(v)); }

	//! Fills \p matching with the original edges which correspond to the edges in m_matching after
	//! expanding it with the correct edges currently contracted in blossoms.
	void getOriginalMatching(std::unordered_set<edge>& matching) {
		// Extend matching in the copy of the graph with the edges in the contracted cycles
		std::vector<Pseudonode*> stack;
		for (auto entry : m_pseudonodes) {
			if (repr(entry.first) == entry.first) {
				stack.push_back(entry.second);
			}
		}

		while (!stack.empty()) {
			auto pseudonode = stack.back();
			stack.pop_back();
			node graphNode = pseudonode->graphNode;
			auto cycle = pseudonode->cycle;
			// Find the edge entering the cycle
			edge matchingEdge = m_matching[graphNode];
			// Find the node in the cycle where the matching edge enters
			node matchedNodeInner = getBaseNode(matchingEdge, graphNode);
			node matchedNode = reprChild(matchedNodeInner);
			// Set the matching edge of the node in case it is a pseudonode and also needs to be
			// expanded later. This breaks the invariant that m_matching always contains exactly the two
			// endpoints pointing to the same edge, but since this is the last step of the algorithm, it
			// doesn't matter.
			m_matching[matchedNode] = matchingEdge;
			auto startIndex = cycle->indexOf(matchedNode);
			auto edges = cycle->edgeOrder();
			// Iterate the edges and add every second one from the perspective of the matchedNode to the
			// matching. Since startIndex is the index of the edge _before_ the matchedNode, we start at
			// index 2.
			for (size_t i = 2; i < edges.size(); i += 2) {
				edge e = edges[(startIndex + i) % edges.size()];
				addToMatching(e);
			}
			// Also expand all children
			for (node v : cycle->nodes()) {
				if (auto child = this->pseudonode(v)) {
					stack.push_back(child);
				}
			}
			expandRepr(pseudonode);
		}
		// Add all respective edges of the original graph to the output parameter
		for (edge e : m_matching) {
			matching.insert(m_graph.original(e));
		}
	}

	//! Adds the edge \p e to the matching.
	void addToMatching(edge e) { m_matching[e->source()] = m_matching[e->target()] = e; }
};

}
}
