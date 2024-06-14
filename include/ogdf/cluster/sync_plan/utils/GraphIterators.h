#pragma once

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/Queue.h>

using namespace ogdf;

class FilteringBFS {
	Queue<node> m_pending;
	NodeArray<bool> m_visited;
	std::function<bool(adjEntry)> m_visit;
	std::function<bool(node)> m_descend;

public:
	template<typename T>
	static bool return_true(T t) {
		return true;
	}

	explicit FilteringBFS() { }

	explicit FilteringBFS(const FilteringBFS& copy) = default;

	explicit FilteringBFS(FilteringBFS&& move) = default;

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
	FilteringBFS* bfs;

public:
	using iterator_category = std::input_iterator_tag;
	using value_type = node;
	using difference_type = void;

	explicit FilteringBFSIterator() : bfs(nullptr) { }

	explicit FilteringBFSIterator(FilteringBFS* bfs) : bfs(bfs) { }

	bool operator==(const FilteringBFSIterator& rhs) const {
		if (bfs) {
			if (rhs.bfs) {
				return bfs == rhs.bfs;
			} else {
				return !bfs->valid();
			}
		} else {
			if (rhs.bfs) {
				return !rhs.bfs->valid();
			} else {
				return true;
			}
		}
	}

	bool operator!=(const FilteringBFSIterator& rhs) const { return !(*this == rhs); }

	node operator*() {
		OGDF_ASSERT(bfs != nullptr);
		return bfs->current();
	}

	FilteringBFSIterator& operator++() {
		OGDF_ASSERT(bfs != nullptr);
		bfs->next();
		return *this;
	}
};

FilteringBFSIterator begin(FilteringBFS& bfs);

FilteringBFSIterator end(FilteringBFS& bfs);
