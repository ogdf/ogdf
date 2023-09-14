/** \file
 * \brief Declaration and implementation of ogdf::MinimumCutStoerWagner
 *
 * \author Mathias Jansen, Stephan Beyer
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

#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/Math.h>
#include <ogdf/basic/PriorityQueue.h>
#include <ogdf/graphalg/MinimumCutModule.h>

namespace ogdf {

/**
 * Computes a minimum cut in a graph.
 *
 * The minimum-cut algorithm according to an approach of Stoer and Wagner 1997
 * is used.
 *
 * @ingroup ga-cut
 */
template<typename T = double>
struct MinimumCutStoerWagner : public MinimumCutModule<T> {
	//! Computes a minimum cut on graph \p G.
	virtual T call(const Graph& G) override {
		EdgeArray<T> weights(G, 1);
		return call(G, weights);
	}

	//! Computes a minimum cut on graph \p G with non-negative \p weights on edges.
	virtual T call(const Graph& G, const EdgeArray<T>& weights) override {
		m_minCut = std::numeric_limits<T>::max();
		m_GC.init(G);
		m_w.init(m_GC);
		m_contractedNodes.init(m_GC);

		for (edge e : m_GC.edges) {
			m_w[e] = weights[m_GC.original(e)];
			OGDF_ASSERT(m_w[e] >= T {});
		}
		for (node v : m_GC.nodes) {
			m_contractedNodes[v].pushBack(m_GC.original(v));
		}

		mainLoop();

		return value();
	}

	/**
	 * Computes the edges defining the computed mincut and returns them.
	 * When calling this method multiple times, the cut edges are only
	 * recomputed if the main min cut algorithm has been called in between.
	 */
	virtual const ArrayBuffer<edge>& edges() override {
		if (!m_cutEdges.empty()) {
			return m_cutEdges;
		}

		NodeArray<bool> inPartition(m_GC.original(), false);

		for (node v : m_partition) {
			inPartition[v] = true;
		}

		for (node v : m_partition) {
			for (adjEntry adj : v->adjEntries) {
				if (!inPartition[adj->twinNode()]) {
					m_cutEdges.push(adj->theEdge());
				}
			}
		}
		return m_cutEdges;
	}

	//! Returns a const-reference to the list of nodes belonging to one side of the bipartition.
	virtual const ArrayBuffer<node>& nodes() override { return m_partition; }

	virtual T value() const override { return m_minCut; }

protected:
	//! Stores the value of the minimum cut
	T m_minCut;

	//! The modifiable version of the input graph (used for contractions)
	GraphCopy m_GC;

	//! An EdgeArray containing the corresponding edge weights.
	EdgeArray<T> m_w;

	//! Store one side of the computed bipartition.
	ArrayBuffer<node> m_partition;

	//! Store cut edges if computed.
	ArrayBuffer<edge> m_cutEdges;

	//! Each node has a list containing the nodes with which it has been contracted.
	//! Because the GraphCopy #m_GC is destroyed during the algorithm, this is
	//! necessary to be able to determine the original nodes in the end.
	NodeArray<List<node>> m_contractedNodes;

	//! Computes minimum cut by invoking minimumCutPhase() O(|V|) times.
	void mainLoop() {
		m_partition.clear();
		m_cutEdges.clear();
		while (m_GC.numberOfNodes() > 1) {
			Math::updateMin(m_minCut, minimumCutPhase());
			if (m_minCut == T {}) { // cannot get better so return early
				return;
			}
		}
	}

	//! Computes and returns the value of the minimum cut of the current phase (iteration).
	T minimumCutPhase();

	//! Contracts the nodes \p s and \p t, i.e., \p s is collapsed to \p t.
	//! The edge (if existing) between \p s and \p t is deleted.
	//! Edges incident to \p s are redirected to \p t.
	//! If parallel edges occur, one of them is deleted and its weight is added to the other one.
	void contraction(node t, node s);
};

template<typename T>
void MinimumCutStoerWagner<T>::contraction(node t, node s) {
	OGDF_ASSERT(t != s);

	// Since nodes in #m_GC correspond to sets of nodes (due to the node contraction),
	// the NodeArray #m_contractedNodes has to be updated.
	// Hence append the list of contracted nodes of s to the list of t.
	m_contractedNodes[t].conc(m_contractedNodes(s));

	// Redirect edges and delete node s
	safeForEach(s->adjEntries, [&](adjEntry adj) {
		edge e = adj->theEdge();
		if (e->source() == t || e->target() == t) {
			m_GC.delEdge(e);
		} else if (e->source() == s) {
			m_GC.moveSource(e, t);
		} else {
			m_GC.moveTarget(e, t);
		}
	});
	m_GC.delNode(s);

	// NodeArray containing the first edge of parallel edges
	NodeArray<edge> incident {m_GC, nullptr};

	// Delete parallel edges and add their weights
	safeForEach(t->adjEntries, [&](adjEntry adj) {
		node v {adj->twinNode()};
		if (v != t) {
			edge e {adj->theEdge()};
			edge& f {incident[v]};
			if (f == nullptr) {
				f = e;
			} else {
				m_w[f] += m_w[e];
				m_GC.delEdge(e);
			}
		}
	});
}

template<typename T>
T MinimumCutStoerWagner<T>::minimumCutPhase() {
	/*
	 * This function computes the mincut of the current phase.
	 * First, nodes are explored successively in descending order of the sum of
	 * their incident edge weights; \a s and \a t are always the two last added nodes.
	 * Afterwards, the current mincut value \a cutOfThePhase is computed, which corresponds to the
	 * sum of the weights of those edges incident to node \a t.
	 * At the end, \a s and \a t are contracted and the \a cutOfThePhase is returned.
	 */

	// A priority queue containing the unexplored nodes.
	// Priorities are the (negative) sum of edge weights incident to the explored nodes.
	PrioritizedMapQueue<node, T> pq(m_GC);
	for (node v : m_GC.nodes) {
		pq.push(v, T {});
	}
	std::function<void(node)> updatePQ {[&](node currentNode) {
		OGDF_ASSERT(!pq.contains(currentNode));
		for (adjEntry adj : currentNode->adjEntries) {
			node v {adj->twinNode()};
			// The code below is at it is due to numerical issues with T=double.
			if (pq.contains(v)) {
				T newPriority {pq.priority(v) - m_w[adj->theEdge()]};
				if (newPriority < pq.priority(v)) {
					pq.decrease(adj->twinNode(), newPriority);
				}
			}
		}
	}};

	// The start node can be chosen arbitrarily. It has no effect on the correctness of the algorithm.
	node s {nullptr};
	node t {pq.topElement()};
	pq.pop();

	updatePQ(t);

	// Successively adding the most tightly connected node.
	while (!pq.empty()) {
		// Find the most tightly connected node and remove it from the priority queue
		node currentNode {pq.topElement()};
		pq.pop();

		// Push the current node to the 2-element list consisting of s and t
		s = currentNode;
		std::swap(s, t);

		// Update the priorities
		updatePQ(currentNode);
	}

	// Contains the mincut value according to the current phase.
	T phaseCut {};
	for (adjEntry adj : t->adjEntries) {
		edge e {adj->theEdge()};
		if (!e->isSelfLoop()) {
			phaseCut += m_w[adj->theEdge()];
		}
	}

	// If the current \a phaseCut is strictly smaller than the global mincut value,
	// the partition defining the mincut has to be updated.
	if (phaseCut < m_minCut) {
		m_partition.clear();
		for (node v : m_contractedNodes[t]) {
			m_partition.push(v);
		}
	}

	// Performing the node contraction of nodes \a s and \a t.
	contraction(t, s);

	return phaseCut;
}

}
