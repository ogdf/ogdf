/** \file
 * \brief An iterator-based BFD through a Graph.
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

// TOOD move to a more central location, add DFS?

#pragma once

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/Queue.h>

namespace ogdf {

/**
 * An iterator-based BFD through a Graph.
 *
 * Allows specifying filters to not visit or descend from certain nodes.
 */
class OGDF_EXPORT FilteringBFS {
	Queue<node> m_pending;
	NodeArray<bool> m_visited;
	std::function<bool(adjEntry)> m_visit;
	std::function<bool(node)> m_descend;

public:
	// iterator traits
	using iterator_category = std::input_iterator_tag;
	using value_type = node;
	using difference_type = std::ptrdiff_t;
	using pointer = node*;
	using reference = node&;

	template<typename T>
	static bool return_true(T t) {
		return true;
	}

	explicit FilteringBFS() { }

	template<typename Container>
	explicit FilteringBFS(const Graph& G, Container& nodes,
			std::function<bool(adjEntry)> visit = return_true<adjEntry>,
			std::function<bool(node)> descend_from = return_true<node>)
		: m_pending(), m_visited(G, false), m_visit(visit), m_descend(descend_from) {
		for (node n : nodes) {
			m_pending.append(n);
		}
	}

	explicit FilteringBFS(const Graph& G, std::initializer_list<node> nodes,
			std::function<bool(adjEntry)> visit = return_true<adjEntry>,
			std::function<bool(node)> descend_from = return_true<node>)
		: m_pending(nodes), m_visited(G, false), m_visit(visit), m_descend(descend_from) { }

	bool operator==(const FilteringBFS& rhs) const {
		return m_pending.getList() == rhs.m_pending.getList();
	}

	bool operator!=(const FilteringBFS& rhs) const {
		return m_pending.getList() != rhs.m_pending.getList();
	}

	FilteringBFS& begin() { return *this; }

	FilteringBFS end() const { return FilteringBFS(); }

	node operator*() {
		OGDF_ASSERT(!m_pending.empty());
		return m_pending.top();
	}

	//! Increment operator (prefix, returns result).
	FilteringBFS& operator++() {
		next();
		return *this;
	}

	//! Increment operator (postfix, returns previous value).
	OGDF_DEPRECATED("Calling FilteringBFS++ will copy the array of visited nodes")

	FilteringBFS operator++(int) {
		FilteringBFS before = *this;
		next();
		return before;
	}

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

	operator bool() const { return valid(); }

	bool valid() const { return !m_pending.empty(); }

	void append(node n) { m_pending.append(n); }

	bool hasVisited(node n) const { return m_visited[n]; }

	int pendingCount() const { return m_pending.size(); }
};

}
