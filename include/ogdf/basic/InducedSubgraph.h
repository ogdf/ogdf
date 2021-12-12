#pragma once


#include <ogdf/basic/GraphObserver.h>

namespace ogdf {

template<typename BaseIterator>
struct filtered_iterator { // TODO IteratorTraits
	typedef std::function<bool(const typename BaseIterator::value_type&)> filter_type;

	filtered_iterator(filter_type filter, BaseIterator base, BaseIterator end = {})
		: _cur(base), _end(end), _filter(filter) {
		while (_cur != _end && !_filter(*_cur)) {
			++_cur;
		}
	}

	filtered_iterator begin() const { return *this; }

	filtered_iterator end() const { return {_filter, _end, _end}; }

	typename BaseIterator::value_type operator*() const { return *_cur; }

	filtered_iterator& operator++() {
		do {
			++_cur;
		} while (_cur != _end && !_filter(*_cur));
		return *this;
	}

	filtered_iterator operator++(int) {
		filtered_iterator copy = *this;
		++*this;
		return copy;
	}

	bool operator==(const filtered_iterator& rhs) const {
		return _cur == rhs._cur && _end == rhs._end
				&& _filter.target_type() == rhs._filter.target_type();
	}

	bool operator!=(const filtered_iterator& rhs) const { return !(rhs == *this); }

private:
	BaseIterator _cur;
	BaseIterator _end;
	filter_type _filter;
};

template<typename BaseIterator>
filtered_iterator<BaseIterator> make_filtered_iterator(
		typename filtered_iterator<BaseIterator>::filter_type filter, BaseIterator base,
		BaseIterator end = {}) {
	return {filter, base, end};
}

template<NodeIter NI, EdgeIter EI, bool copyEmbedding, bool copyIDs, bool notifyObservers>
std::pair<int, int> Graph::insert(const NI& nodesBegin, const NI& nodesEnd, const EI& edgesBegin,
		const EI& edgesEnd, NodeArray<node>& nodeMap,
		EdgeArray<edge>& edgeMap) { // TODO use template magic to switch to a faster implementation, if possible
	// TODO reserve size in m_regNodeArrays and m_regEdgeArrays
	OGDF_ASSERT(nodeMap.registeredAt());
	OGDF_ASSERT(edgeMap.registeredAt());
	if (nodesBegin == nodesEnd) {
		return {0, 0};
	}
	int newNodes = 0, newEdges = 0;

	for (auto it = nodesBegin; it != nodesEnd; ++it) {
		node vG = *it;
		if (copyIDs) {
			m_nodeIdCount = max(m_nodeIdCount, vG->index() + 1);
		}
		node v = nodeMap[vG] = pureNewNode(copyIDs ? vG->index() : m_nodeIdCount++);
		newNodes++;
		if (notifyObservers) {
			m_regNodeArrays.keyAdded(v);
			for (GraphObserver* obs : m_regObservers) {
				obs->nodeAdded(v);
			}
		}
	}

	if (edgesBegin == edgesEnd) {
		return {newNodes, newEdges};
	}

	for (auto it = edgesBegin; it != edgesEnd; ++it) {
		edge eG = *it;
		node src = nodeMap[eG->source()];
		node tgt = nodeMap[eG->target()];
		if (src == nullptr || tgt == nullptr) {
			continue;
		}
		if (copyIDs) {
			m_edgeIdCount = max(m_edgeIdCount, eG->index() + 1);
		}
		edge e = edgeMap[eG] = pureNewEdge(src, tgt, copyIDs ? eG->index() : m_edgeIdCount++);
		newEdges++;
		if (!copyEmbedding) {
			src->adjEntries.pushBack(e->m_adjSrc);
			tgt->adjEntries.pushBack(e->m_adjTgt);
			if (notifyObservers) {
				m_regEdgeArrays.keyAdded(e);
				m_regAdjArrays.keyAdded(e->adjSource());
				for (GraphObserver* obs : m_regObservers) {
					obs->edgeAdded(e);
				}
			}
		}
	}

	if (!copyEmbedding) {
#ifdef OGDF_HEAVY_DEBUG
		consistencyCheck();
#endif
		return {newNodes, newEdges};
	}

	for (auto it = nodesBegin; it != nodesEnd; ++it) {
		node vG = *it;
		node v = nodeMap[vG];
		for (adjEntry adjG : vG->adjEntries) {
			edge e = adjG->m_edge;
			edge eC = edgeMap[e];
			if (eC == nullptr) {
				continue;
			}
			adjEntry adj = adjG->isSource() ? eC->adjSource() : eC->adjTarget();
			v->adjEntries.pushBack(adj);
		}
	}

	// notify observers of added edges after adjEntries are initialized
	if (notifyObservers) { // FIXME call registry.reserve(), maybe use && !m_regObservers.empty()
		for (auto it = edgesBegin; it != edgesEnd; ++it) {
			edge eG = *it;
			edge e = edgeMap[eG];
			if (e == nullptr) {
				continue;
			}
			m_regEdgeArrays.keyAdded(e);
			m_regAdjArrays.keyAdded(e->adjSource());
			for (GraphObserver* obs : m_regObservers) {
				obs->edgeAdded(e);
			}
		}
	}

#ifdef OGDF_HEAVY_DEBUG
	consistencyCheck();
#endif

	return {newNodes, newEdges};
}

template<NodeFilter NF, EdgeFilter EF, bool copyEmbedding, bool copyIDs, bool notifyObservers>
std::pair<int, int> Graph::insert(const Graph& G, const NF& nodeFilter, const EF& edgeFilter,
		NodeArray<node>& nodeMap, EdgeArray<edge>& edgeMap) {
	if (!nodeMap.registeredAt()) {
		nodeMap.init(G);
	}
	OGDF_ASSERT(nodeMap.registeredAt()->graphOf() == &G);
	if (!edgeMap.registeredAt()) {
		edgeMap.init(G);
	}
	OGDF_ASSERT(edgeMap.registeredAt()->graphOf() == &G);
	filtered_iterator<node_iterator> nodes_it {nodeFilter, G.nodes.begin(), G.nodes.end()};
	filtered_iterator<edge_iterator> edges_it {edgeFilter, G.edges.begin(), G.edges.end()};
	return insert<filtered_iterator<node_iterator>, filtered_iterator<edge_iterator>, copyEmbedding,
			copyIDs, notifyObservers>(nodes_it, nodes_it.end(), edges_it, edges_it.end(), nodeMap,
			edgeMap);
}
}