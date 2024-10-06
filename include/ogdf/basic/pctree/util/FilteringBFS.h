/** \file
 * \brief An iterator-based BFS through a Graph. TODO should be moved to a central location; add DFS?
 *
 * \author Simon D. Fink <ogdf@niko.fink.bayern>
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
#include <ogdf/basic/Queue.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/internal/copy_move.h>
#include <ogdf/basic/internal/graph_iterators.h>

#include <functional>
#include <initializer_list>
#include <iosfwd>
#include <iterator>

namespace ogdf {

class FilteringBFSIterator;

/**
 * An iterator-based BFS through a Graph.
 *
 * Allows specifying filters to not visit or descend from certain nodes.
 */
// see also FilteringPCTreeDFS/BFS
class OGDF_EXPORT FilteringBFS {
	Queue<node> m_pending;
	NodeArray<bool> m_visited;
	std::function<bool(adjEntry)> m_visit;
	std::function<bool(node)> m_descend;

public:
	template<typename T>
	static bool return_true(T t) {
		return true;
	}

	explicit FilteringBFS() = default;

	OGDF_DEFAULT_COPY(FilteringBFS)
	OGDF_DEFAULT_MOVE(FilteringBFS)

	template<typename Container>
	explicit FilteringBFS(const Graph& G, Container& nodes,
			const std::function<bool(adjEntry)>& visit = return_true<adjEntry>,
			const std::function<bool(node)>& descend_from = return_true<node>)
		: m_pending(), m_visited(G, false), m_visit(visit), m_descend(descend_from) {
		for (node n : nodes) {
			m_pending.append(n);
		}
	}

	explicit FilteringBFS(const Graph& G, std::initializer_list<node> nodes,
			const std::function<bool(adjEntry)>& visit = return_true<adjEntry>,
			const std::function<bool(node)>& descend_from = return_true<node>)
		: m_pending(nodes), m_visited(G, false), m_visit(visit), m_descend(descend_from) { }

	bool operator==(const FilteringBFS& rhs) const {
		return m_pending.getList() == rhs.m_pending.getList();
	}

	bool operator!=(const FilteringBFS& rhs) const {
		return m_pending.getList() != rhs.m_pending.getList();
	}

	FilteringBFSIterator begin();

	FilteringBFSIterator end();

	void next() {
		OGDF_ASSERT(!m_pending.empty());
		node n = m_pending.pop();
		OGDF_ASSERT(!m_visited[n]);
		m_visited[n] = true;
		if (m_descend(n)) {
			for (adjEntry adj : n->adjEntries) {
				node twin = adj->twinNode();
				if (!m_visited[twin] && m_visit(adj)) {
					m_pending.append(twin);
				}
			}
		}
		while (!m_pending.empty() && m_visited[m_pending.top()]) {
			m_pending.pop();
		}
	}

	node current() {
		OGDF_ASSERT(!m_pending.empty());
		return m_pending.top();
	}

	operator bool() const { return valid(); }

	bool valid() const { return !m_pending.empty(); }

	void append(node n) {
		m_visited[n] = false;
		m_pending.append(n);
	}

	bool hasVisited(node n) const { return m_visited[n]; }

	bool willVisitTarget(adjEntry adj) const { return m_visit(adj); }

	bool willDescendFrom(node n) const { return m_descend(n); }

	void setVisitFilter(const std::function<bool(adjEntry)>& mVisit) { m_visit = mVisit; }

	void setDescendFilter(const std::function<bool(node)>& mDescend) { m_descend = mDescend; }

	int pendingCount() const { return m_pending.size(); }
};

class FilteringBFSIterator {
	FilteringBFS* m_bfs;

public:
	// iterator traits
	using iterator_category = std::input_iterator_tag;
	using value_type = node;
	using difference_type = std::ptrdiff_t;
	using pointer = node*;
	using reference = node&;

	explicit FilteringBFSIterator() : m_bfs(nullptr) { }

	explicit FilteringBFSIterator(FilteringBFS* bfs) : m_bfs(bfs) { }

	bool operator==(const FilteringBFSIterator& rhs) const {
		if (m_bfs) {
			if (rhs.m_bfs) {
				return m_bfs == rhs.m_bfs;
			} else {
				return !m_bfs->valid();
			}
		} else {
			if (rhs.m_bfs) {
				return !rhs.m_bfs->valid();
			} else {
				return true;
			}
		}
	}

	bool operator!=(const FilteringBFSIterator& rhs) const { return !(*this == rhs); }

	node operator*() {
		OGDF_ASSERT(m_bfs != nullptr);
		return m_bfs->current();
	}

	FilteringBFSIterator& operator++() {
		OGDF_ASSERT(m_bfs != nullptr);
		m_bfs->next();
		return *this;
	}
};

inline FilteringBFSIterator FilteringBFS::begin() { return FilteringBFSIterator(this); }

inline FilteringBFSIterator FilteringBFS::end() { return FilteringBFSIterator(nullptr); }

}
