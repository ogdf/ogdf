/** \file
 * \brief Implementation of the original Blossom algorithm by Edmonds (1965).
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
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/Logger.h>
#include <ogdf/graphalg/MatchingModule.h>
#include <ogdf/graphalg/matching_blossom/AlternatingTree.h>
#include <ogdf/graphalg/matching_blossom/BlossomHelper.h>
#include <ogdf/graphalg/matching_blossom/Cycle.h>
#include <ogdf/graphalg/matching_blossom/Pseudonode.h>
#include <ogdf/graphalg/matching_blossom/utils.h>

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace ogdf {

/**
 * Implementation of the Blossom-I algorithm by Edmonds (1965) for Minimum
 * Weight Perfect Matching.
 * Use MatchingBlossomV for a more sophisticated and faster version of this
 * algorithm.
 */
template<typename TWeight>
class MatchingBlossom : public MatchingModule<TWeight> {
	using Logger::lout;

	//! Helper class to store the current state of the algorithm.
	BlossomHelper<TWeight> m_helper;

	//! Pointer to the HiddenEdgeSet containing all non-equality edges.
	std::unique_ptr<Graph::HiddenEdgeSet> m_nonEqualityEdgesHiddenSet;

	//! Set of all non-equality edges.
	std::unordered_set<edge> m_nonEqualityEdges;

	//! All nodes currently present in the graph
	std::unordered_set<node> m_graphNodes;

	//! All nodes which are not part of the matching yet
	std::unordered_set<node> m_unmatchedNodes;

	//! All top-level pseudonodes
	std::unordered_set<node> m_pseudonodes;

	//! The alternating tree
	std::unique_ptr<AlternatingTree<TWeight>> m_tree;

#ifdef OGDF_HEAVY_DEBUG
	//! Debug function to assert that all helper variables are in a consistent state.
	void assertConsistency() {
		if (m_unmatchedNodes.empty()) {
			return;
		}
		// graph nodes & euqality edges
		int matchedNodes = 0;
		for (node v : m_graphNodes) {
			if (m_helper.matching(v)) {
				matchedNodes++;
			}
			OGDF_ASSERT(v->graphOf() == &m_helper.graph());
		}
		OGDF_ASSERT((size_t)m_nonEqualityEdgesHiddenSet->size() == m_nonEqualityEdges.size());
		for (edge e : m_helper.graph().edges) {
			bool isGraphNode = m_graphNodes.find(e->source()) != m_graphNodes.end();
			OGDF_ASSERT(isGraphNode == (m_graphNodes.find(e->target()) != m_graphNodes.end()));
			if (isGraphNode) {
				OGDF_ASSERT((m_nonEqualityEdges.find(e) == m_nonEqualityEdges.end())
						== m_helper.isEqualityEdge(e));
			}
		}
		// tree structure
		for (node v : m_tree->evenNodes) {
			// all even nodes except the root have a corresponding matching edge
			if (v != m_tree->root()) {
				OGDF_ASSERT(m_helper.matching(v) == m_tree->evenBackEdge(v));
			}
			// and no back edge in the tree
			OGDF_ASSERT(!m_tree->isOdd(v));
		}
		for (node v : m_tree->oddNodes) {
			OGDF_ASSERT(!m_tree->isEven(v));
			// all odd nodes have both a back edge in the tree and a forward matching
			// edge, since the matching edges are always defined for both ends
			node w = m_helper.matching(v)->opposite(v);
			OGDF_ASSERT(m_tree->isEven(w));
		}
		// assert that all matching edges are correctly defined on both ends
		for (node v : m_helper.graph().nodes) {
			if (edge matchingEdge = m_helper.matching(v)) {
				OGDF_ASSERT(matchingEdge->isIncident(v));
				OGDF_ASSERT(m_helper.matching(matchingEdge->opposite(v)) == matchingEdge);
			}
		}
		// matching
		long unsigned int hiddenNodes = 0;
		for (auto entry : m_helper.pseudonodes()) {
			auto pseudonode = entry.second;
			hiddenNodes += pseudonode->cycle->nodes().size();
			for (node v : pseudonode->cycle->nodes()) {
				OGDF_ASSERT(m_helper.matching(v) == nullptr);
			}
		}
		OGDF_ASSERT(matchedNodes + m_unmatchedNodes.size() + hiddenNodes
				== (size_t)m_helper.graph().numberOfNodes());
		if (m_tree->hasRoot()) {
			OGDF_ASSERT(m_unmatchedNodes.find(m_tree->root()) != m_unmatchedNodes.end());
		}
	}
#endif

	//! Helper function to hide all non-equality edges
	void hideNonEqualityEdges() {
		for (edge e : m_nonEqualityEdges) {
			m_nonEqualityEdgesHiddenSet->hide(e);
		}
	}

	//! Helper function to restore all non-equality edges
	void restoreNonEqualityEdges() { m_nonEqualityEdgesHiddenSet->restore(); }

	/**
	 * @brief Helper function to get a new root for the tree.
	 *
	 * \pre There must still be unmatched nodes left.
	 *
	 * @return the new root
	 */
	node getNewRoot() { return *m_unmatchedNodes.begin(); }

public:
	/**
	 * @brief Construct a MatchingBlossom instance.
	 *
	 * @param greedyInit whether or not to use the greedy initialization
	 */
	MatchingBlossom(bool greedyInit = true) : m_helper(greedyInit) { }

private:
	bool doCall(const Graph& G, const EdgeArray<TWeight>& weights,
			std::unordered_set<edge>& matching) {
		return _doCall(G, weights, matching);
	}

	bool doCall(const GraphAttributes& GA, std::unordered_set<edge>& matching) {
		return _doCall(GA.constGraph(), GA, matching);
	}

	//! Helper for the main call function since abstract functions cannot be templated.
	template<class WeightContainer>
	bool _doCall(const Graph& G, const WeightContainer& weights, std::unordered_set<edge>& matching) {
		// init
		lout() << std::endl;
		if (!m_helper.init(G, weights)) {
			return false;
		}
		m_nonEqualityEdgesHiddenSet.reset(new Graph::HiddenEdgeSet(m_helper.graph()));
		m_nonEqualityEdges.clear();
		for (edge e : m_helper.graph().edges) {
			// hide all non-equality edges
			if (!m_helper.isEqualityEdge(e)) {
				m_nonEqualityEdges.insert(e);
			}
		}
		hideNonEqualityEdges();
		m_graphNodes.clear();
		m_unmatchedNodes.clear();
		for (node v : m_helper.graph().nodes) {
			m_graphNodes.insert(v);
			if (!m_helper.matching(v)) {
				m_unmatchedNodes.insert(v);
			}
		}
		m_tree.reset(new AlternatingTree<TWeight>(m_helper));

		return findMatching(matching);
	}

	/**
	 * @brief Main function of the algorithm. Finds a minimum weight perfect matching in the graph
	 * if one exists.
	 *
	 * @returns whether a matching was found
	 */
	bool findMatching(std::unordered_set<edge>& matching) {
		// main loop
		while (!m_unmatchedNodes.empty()) {
#ifdef OGDF_HEAVY_DEBUG
			assertConsistency();
#endif
			if (!m_tree->hasRoot()) {
				// we need to generate a new root for the tree
				m_tree->reset(getNewRoot());
				lout() << "new root: " << m_tree->root() << std::endl;
			}
			if (!primalChange() && !dualChange()) {
				return false;
			}
		}
		// matching found
		m_helper.getOriginalMatching(matching);
		return true;
	}

	//! Executes a primal change step.
	bool primalChange() {
#ifdef OGDF_HEAVY_DEBUG
		std::unordered_map<edge, TWeight> edgeWeights;
		m_helper.getReducedEdgeWeights(edgeWeights);
#endif
		bool result = findMatchingAugmentation() || findTreeAugmentation() || findShrinkableCycle()
				|| findExpandablePseudonode();
#ifdef OGDF_HEAVY_DEBUG
		m_helper.assertConsistentEdgeWeights(edgeWeights);
#endif
		return result;
	}

	/**
	 * @brief Finds and executes a matching augmentation, if possible.
	 *
	 * @return whether the matching was augmented
	 */
	bool findMatchingAugmentation() {
		for (node u : m_tree->evenNodes) {
			for (adjEntry adj : u->adjEntries) {
				node v = adj->twinNode();
				if (!m_helper.matching(v) && v != m_tree->root()) {
					// found free edge, augment matching
					m_tree->augmentMatching(adj->theEdge());
					m_unmatchedNodes.erase(m_tree->root());
					m_unmatchedNodes.erase(v);
					m_tree->reset();
					lout() << "augmented matching with " << adj->theEdge() << std::endl;
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * @brief Finds and executes a tree augmentation, if possible.
	 *
	 * @return whether the tree was augmented
	 */
	bool findTreeAugmentation() {
		for (node u : m_tree->evenNodes) {
			for (adjEntry adj : u->adjEntries) {
				node v = adj->twinNode();
				edge matchingEdge = m_helper.matching(v);
				if (matchingEdge && !m_tree->contains(v)) {
					m_tree->grow(u, adj->theEdge(), matchingEdge);
					lout() << "augmented tree with " << adj->theEdge() << std::endl;
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * @brief Finds and shrinks an odd cycle, if possible.
	 *
	 * @return whether a cycle was shrunken
	 */
	bool findShrinkableCycle() {
		for (node u : m_tree->evenNodes) {
			for (adjEntry adj : u->adjEntries) {
				node v = adj->twinNode();
				if (m_tree->isEven(v)) {
					// found an odd cycle: shrink it
					shrink(adj->theEdge());
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * @brief Finds and expands an odd pseudonode, if possible.
	 *
	 * @return whether a pseudonode was expanded
	 */
	bool findExpandablePseudonode() {
		for (node v : m_pseudonodes) {
			if (m_tree->isOdd(v) && m_helper.isZeroCostNode(v)) {
				expand(m_helper.pseudonode(v));
				return true;
			}
		}
		return false;
	}

	//! Executes a dual change step.
	bool dualChange() {
		TWeight delta = infinity<TWeight>();
		// restore all edges to find potential candidates
		restoreNonEqualityEdges();
		// find ++/+0 edges
		for (node u : m_tree->evenNodes) {
			for (adjEntry adj : u->adjEntries) {
				node v = adj->twinNode();
				TWeight potentialChange = m_helper.getReducedWeight(adj->theEdge());
				if (m_tree->isEven(v)) {
					// ++ edge: max dual change has to be halfed
					potentialChange *= 0.5;
				} else if (m_tree->isOdd(v)) {
					// +- edge: ignore
					continue;
				}
				// (modified) reduced weight gives an upper bound to the dual change
				delta = std::min(delta, potentialChange);
			}
		}
		// find - pseudonodes
		for (node v : m_pseudonodes) {
			if (m_tree->isOdd(v)) {
				delta = std::min(delta, m_helper.y(v));
			}
		}
		// no dual change possible, so no matching can be found
		if (delta == infinity<TWeight>()) {
			return false;
		}
		OGDF_ASSERT(delta > 0);

		// do update of dual variables
		for (node v : m_tree->evenNodes) {
			m_helper.y(v) += delta;
		}
		for (node v : m_tree->oddNodes) {
			m_helper.y(v) -= delta;
		}
		// update non-equality edges
		for (node v : m_tree->evenNodes) {
			for (adjEntry adj : v->adjEntries) {
				edge e = adj->theEdge();
				if (m_helper.isEqualityEdge(e)) {
					auto it = m_nonEqualityEdges.find(e);
					if (it != m_nonEqualityEdges.end()) {
						m_nonEqualityEdges.erase(it);
					}
				}
			}
		}
		for (node u : m_tree->oddNodes) {
			// no adjacent edge to a node outside of the current tree can be an equality edge anymore
			for (adjEntry adj : u->adjEntries) {
				node v = adj->twinNode();
				edge e = adj->theEdge();
				if (!m_tree->isEven(v)) {
					m_nonEqualityEdges.insert(e);
				}
			}
		}
		hideNonEqualityEdges();
		lout() << "dual change: " << delta << std::endl;
		return true;
	}

	//! Shrink the odd cycle induced by the current tree together with \p cycleEdge.
	void shrink(edge cycleEdge) {
		restoreNonEqualityEdges();
		Cycle* cycle = m_tree->getCycle(cycleEdge);
		std::vector<std::tuple<edge, bool>> selfLoops;
		Pseudonode* pseudonode = m_tree->shrink(cycle, selfLoops);
		node newNode = pseudonode->graphNode;

		for (node u : pseudonode->cycle->nodes()) {
			m_unmatchedNodes.erase(u);
			m_graphNodes.erase(u);
		}
		for (node u : pseudonode->cycle->nodes()) {
			for (auto adj : u->adjEntries) {
				edge e = adj->theEdge();
				auto it = m_nonEqualityEdges.find(e);
				if (it != m_nonEqualityEdges.end()) {
					m_nonEqualityEdges.erase(it);
				}
			}
		}
		m_graphNodes.insert(newNode);
		m_pseudonodes.insert(newNode);
		if (m_helper.matching(newNode) == nullptr) {
			m_unmatchedNodes.insert(newNode);
		}
		lout() << "shrunk odd cycle of length " << cycle->nodes().size() << " at " << cycleEdge
			   << " into pseudonode " << pseudonode->graphNode << std::endl;
		hideNonEqualityEdges();
	}

	//! Expand the pseudonode \p pseudonode.
	void expand(Pseudonode* pseudonode) {
		node graphNode = pseudonode->graphNode;
		m_pseudonodes.erase(graphNode);
		m_graphNodes.erase(graphNode);
		OGDF_ASSERT(m_unmatchedNodes.empty() || m_helper.isZeroCostNode(graphNode));
		restoreNonEqualityEdges();
		int pseudonodeIndex = graphNode->index();
		m_tree->expand(pseudonode);
		for (node u : pseudonode->cycle->nodes()) {
			m_graphNodes.insert(u);
		}
		for (node u : pseudonode->cycle->nodes()) {
			for (auto adj : u->adjEntries) {
				edge e = adj->theEdge();
				if (!m_helper.isEqualityEdge(e)) {
					m_nonEqualityEdges.insert(e);
				}
			}
		}
		delete pseudonode;
		hideNonEqualityEdges();
		lout() << "expanded pseudonode " << pseudonodeIndex << std::endl;
	}
};

}
