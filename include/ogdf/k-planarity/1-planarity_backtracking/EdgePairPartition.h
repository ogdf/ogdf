/** \file
 * \brief Partial solutions for backtracking 1-Planarity.
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

#pragma once

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/GraphSets.h>
#include <ogdf/basic/basic.h>

#include <list>
#include <memory>
#include <set>
#include <stack>
#include <vector>

namespace ogdf::oneplan_backtracking {

//! A pair of distinct graph elements, ordered by their index.
template<typename GraphElementPointer>
class OrderedPair {
private:
	GraphElementPointer m_first, m_second;

public:
	//! Creates a new pair that consists of \p a and \p b.
	OrderedPair(GraphElementPointer a, GraphElementPointer b) {
		OGDF_ASSERT(a && b);
		OGDF_ASSERT(a->graphOf() == b->graphOf());
		OGDF_ASSERT(a != b);
		if (a->index() < b->index()) {
			m_first = a;
			m_second = b;
		} else {
			m_first = b;
			m_second = a;
		}
	}

	//! Returns the element with smaller index.
	GraphElementPointer first() const { return m_first; }

	//! Returns the element with larger index.
	GraphElementPointer second() const { return m_second; }

	//! Returns whether \p el is contained in the pair.
	bool contains(GraphElementPointer el) const { return m_first == el || m_second == el; }

	bool operator==(const OrderedPair& rhs) const {
		OGDF_ASSERT(graphOf() == rhs.graphOf());
		OGDF_ASSERT(m_first->index() <= m_second->index());
		OGDF_ASSERT(rhs.m_first->index() <= rhs.m_second->index());

		return m_first == rhs.m_first && m_second == rhs.m_second;
	}

	bool operator!=(const OrderedPair& rhs) const { return !operator==(rhs); }

	bool operator<(const OrderedPair& rhs) const {
		OGDF_ASSERT(graphOf() == rhs.graphOf());
		OGDF_ASSERT(m_first->index() <= m_second->index());
		OGDF_ASSERT(rhs.m_first->index() <= rhs.m_second->index());

		return m_first->index() < rhs.m_first->index()
				|| ((m_first->index() == rhs.m_first->index())
						&& (m_second->index() < rhs.m_second->index()));
	}

#ifdef OGDF_DEBUG
	const Graph* graphOf() const { return m_first->graphOf(); }
#endif
};

using VertexPair = OrderedPair<node>;
using EdgePair = OrderedPair<edge>;

//! The different modes for 1-Planarity.
enum class OGDF_EXPORT OneplanMode {
	Normal, //!< 1-Planarity
	NIC, //!< 1-Planarity, where two crossed edge pairs may share at most one vertex (NIC-Planarity)
	IC //!< 1-Planarity, where no two crossed edges may share a vertex (IC-Planarity)
};

//! A partial solution for a 1-Planarity instance, representing a node in the backtracking tree.
/**
 * All independent pairs of the input graph are partitioned into crossing, uncrossable, and free
 * (i.e., still allowed to cross). Maintains a stack of transactions that can be later undone.
 */
class OGDF_EXPORT EdgePairPartition {
private:
	struct UndoInformation {
		std::set<EdgePair> m_crossings;
		std::set<VertexPair> m_kiteEdges;
		std::set<EdgePair> m_nonCrossingPairs;
		std::list<EdgePair> m_todoCrossings;
		std::unique_ptr<EdgePair> m_previousCrossing = nullptr;
	};

	OneplanMode m_mode = OneplanMode::Normal;
	const Graph& m_graph;
	EdgeSet m_crossedEdges;
	std::set<VertexPair> m_kiteEdges;
	std::set<EdgePair> m_crossingEdgePairs;
	std::set<EdgePair> m_freeEdgePairs;

	std::stack<UndoInformation> m_undoInformation;

public:
	//! Creates a new EdgePairPartition, where all pairs of independent edges are crossable.
	/**
	 * @param g The input graph.
	 * @param m The kind of 1-Planarity represented
	 */
	explicit EdgePairPartition(const Graph& g, OneplanMode m = OneplanMode::Normal)
		: m_mode(m), m_graph(g), m_crossedEdges(m_graph) {
		for (edge e1 : m_graph.edges) {
			for (edge e2 : m_graph.edges) {
				if (!e1->isAdjacent(e2)) {
					m_freeEdgePairs.emplace(e1, e2);
				}
			}
		}
	}

	//! Creates a new EdgePairPartition as a copy of \p other, but does not copy its transaction stack.
	EdgePairPartition(const EdgePairPartition& other)
		: m_mode(other.m_mode)
		, m_graph(other.m_graph)
		, m_crossedEdges(other.m_crossedEdges)
		, m_kiteEdges(other.m_kiteEdges)
		, m_crossingEdgePairs(other.m_crossingEdgePairs)
		, m_freeEdgePairs(other.m_freeEdgePairs) { }

	//! Returns the associated graph.
	const Graph& graph() const { return m_graph; }

	//! Returns the set of kite edges.
	const std::set<VertexPair>& kiteEdges() const { return m_kiteEdges; }

	//! Returns the set of all crossed edge pairs.
	const std::set<EdgePair>& crossingEdgePairs() const { return m_crossingEdgePairs; }

	//! Starts a new transaction that can later be undone.
	void startTransaction() { m_undoInformation.emplace(); }

	//! Starts a new transaction and associates an exhaustive set of crossings that will be branched over.
	void startTransaction(std::vector<EdgePair>& todoCrossings) {
		UndoInformation& ui = m_undoInformation.emplace();
		ui.m_todoCrossings.assign(todoCrossings.begin(), todoCrossings.end());
	}

	//! Reverts the previous transaction.
	void undoTransaction();

	//! Returns whether there are crossings that should be branched over.
	bool hasTodoCrossing() {
		if (m_undoInformation.empty()) {
			return false;
		}
		return !m_undoInformation.top().m_todoCrossings.empty();
	}

	//! Removes and returns the next stored crossing that should be considered.
	EdgePair getNextTodoCrossing() {
		OGDF_ASSERT(hasTodoCrossing());
		EdgePair ep = m_undoInformation.top().m_todoCrossings.front();
		m_undoInformation.top().m_todoCrossings.pop_front();
		m_undoInformation.top().m_previousCrossing = std::make_unique<EdgePair>(ep);
		return ep;
	}

	//! Crosses \p pair and computes the corresponding kite edges.
	void crossEdgePair(const EdgePair& pair);

	//! Marks \p pair as non-crossing.
	void setNonCrossing(const EdgePair& pair) {
		if (m_freeEdgePairs.count(pair) > 0) {
			m_freeEdgePairs.erase(pair);
			m_undoInformation.top().m_nonCrossingPairs.insert(pair);
		}
	}

	//! Returns whether \p pair is still allowed to cross.
	bool isFree(const EdgePair& pair) const {
		OGDF_ASSERT(pair.graphOf() == &m_graph);
		if (m_crossedEdges.contains(pair.first()) || m_crossedEdges.contains(pair.second())) {
			return false;
		}
		if (m_kiteEdges.count({pair.first()->source(), pair.first()->target()}) > 0
				|| m_kiteEdges.count({pair.second()->source(), pair.second()->target()}) > 0) {
			return false;
		}
		return m_freeEdgePairs.count(pair) > 0;
	}

	//! Returns whether there exists an edge that \p e may cross.
	bool isFree(edge e) const {
		if (m_kiteEdges.count({e->source(), e->target()}) > 0 || m_crossedEdges.contains(e)) {
			return false;
		}
		for (EdgePair ep : m_freeEdgePairs) {
			if (ep.contains(e)) {
				return true;
			}
		}
		return false;
	}

	//! Returns whether \p e is part of a crossing.
	bool isCrossed(edge e) const {
		for (EdgePair ep : m_crossingEdgePairs) {
			if (ep.contains(e)) {
				return true;
			}
		}
		return false;
	}

	//! Returns an edge between \p a and \p b, or null
	static edge getEdgeBetween(node a, node b) {
		for (adjEntry adj : a->adjEntries) {
			if (adj->twinNode() == b) {
				return adj->theEdge();
			}
		}
		return nullptr;
	}

private:
	void computeKiteEdges(const EdgePair& crossedPair);

	void addKiteEdge(edge e) {
		OGDF_ASSERT(e->graphOf() == &m_graph);
		if (!m_crossedEdges.contains(e)) {
			addKiteEdge({e->source(), e->target()});
		}
	}

	void addKiteEdge(const VertexPair& vp) {
		if (m_kiteEdges.count(vp) == 0) {
			m_kiteEdges.insert(vp);
			m_undoInformation.top().m_kiteEdges.insert(vp);
		}
	}
};
}
