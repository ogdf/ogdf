#pragma once


#include <ogdf/basic/GraphObserver.h>

namespace ogdf {

namespace internal {
template<typename...>
using void_t = void;

template<class T, class = void>
struct is_iterator : std::false_type { };

template<class T>
struct is_iterator<T, void_t<typename std::iterator_traits<T>::iterator_category>> : std::true_type {
};

template<class It>
typename std::iterator_traits<It>::difference_type do_guess_dist(It first, It last,
		std::input_iterator_tag) {
	return 0;
}

template<class It>
typename std::iterator_traits<It>::difference_type do_guess_dist(It first, It last,
		std::random_access_iterator_tag) {
	return last - first;
}

template<class It>
typename std::enable_if<!is_iterator<It>::value, int>::type guess_dist(It first, It last) {
	return 0;
}

template<class It>
typename std::enable_if<is_iterator<It>::value, typename std::iterator_traits<It>::difference_type>::type
guess_dist(It first, It last) {
	return do_guess_dist(first, last, typename std::iterator_traits<It>::iterator_category());
}
}

template<OGDF_NODE_ITER NI, OGDF_EDGE_ITER EI, bool copyEmbedding, bool copyIDs, bool notifyObservers>
std::pair<int, int> Graph::insert(const NI& nodesBegin, const NI& nodesEnd, const EI& edgesBegin,
		const EI& edgesEnd, NodeArray<node>& nodeMap, EdgeArray<edge>& edgeMap) {
	OGDF_ASSERT(nodeMap.valid());
	OGDF_ASSERT(edgeMap.valid());
	OGDF_ASSERT(nodeMap.graphOf() == edgeMap.graphOf());
	int newNodes = 0, newEdges = 0;
	void* cbData = preInsert(copyEmbedding, copyIDs, notifyObservers, false, nodeMap, edgeMap,
			&newNodes, &newEdges);
	if (nodesBegin == nodesEnd) {
		postInsert(cbData, newNodes, newEdges);
		return {newNodes, newEdges};
	}

	int guessedNodes = internal::guess_dist(nodesBegin, nodesEnd);
	if (guessedNodes > 0 && notifyObservers) {
		m_regNodeArrays.reserveSpace(guessedNodes);
	}

	for (auto it = nodesBegin; it != nodesEnd; ++it) {
		node vG = *it;
		if (copyIDs) {
			m_nodeIdCount = max(m_nodeIdCount, vG->index() + 1);
		}
		// nodeMap[vG] is overwritten if it is != nullptr
		node v = nodeMap[vG] = pureNewNode(copyIDs ? vG->index() : m_nodeIdCount++);
		newNodes++;
		if (notifyObservers) {
			m_regNodeArrays.keyAdded(v);
			nodeInserted(cbData, vG, v);
			for (GraphObserver* obs : getObservers()) {
				obs->nodeAdded(v);
			}
		}
	}

	if (edgesBegin == edgesEnd) {
		postInsert(cbData, newNodes, newEdges);
		return {newNodes, newEdges};
	}

	if (!copyEmbedding && notifyObservers) {
		int guessedEdges = internal::guess_dist(edgesBegin, edgesEnd);
		if (guessedEdges > 0) {
			m_regEdgeArrays.reserveSpace(guessedEdges);
			m_regAdjArrays.reserveSpace(guessedEdges); // registry adds factor 2 in calculateArraySize
		}
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
		// edgeMap[eG] is overwritten if it is != nullptr
		edge e = edgeMap[eG] = pureNewEdge(src, tgt, copyIDs ? eG->index() : m_edgeIdCount++);
		newEdges++;
		if (!copyEmbedding) {
			src->adjEntries.pushBack(e->m_adjSrc);
			tgt->adjEntries.pushBack(e->m_adjTgt);
			if (notifyObservers) {
				m_regEdgeArrays.keyAdded(e);
				m_regAdjArrays.keyAdded(e->adjSource());
				edgeInserted(cbData, eG, e);
				for (GraphObserver* obs : getObservers()) {
					obs->edgeAdded(e);
				}
			}
		}
	}

	if (!copyEmbedding) {
#ifdef OGDF_HEAVY_DEBUG
		consistencyCheck();
#endif
		postInsert(cbData, newNodes, newEdges);
		return {newNodes, newEdges};
	}

	for (auto it = nodesBegin; it != nodesEnd; ++it) {
		node vG = *it;
		node v = nodeMap[vG];
		for (adjEntry adjG : vG->adjEntries) {
			edge eG = adjG->m_edge;
			edge e = edgeMap[eG];
			if (e == nullptr) {
				continue;
			}
			adjEntry adj = adjG->isSource() ? e->adjSource() : e->adjTarget();
			// edgeMap[eG] might be an old entry that was already inserted into the list
			// so check whether adj is already in the list, indicated by having succ or pred,
			// or being the only entry in the list
			if (adj->succ() != nullptr || adj->pred() != nullptr || v->adjEntries.head() == adj) {
				continue;
			}
			v->adjEntries.pushBack(adj);
		}
	}

	// notify observers of added edges after adjEntries are initialized
	if (notifyObservers) {
		m_regEdgeArrays.reserveSpace(newEdges);
		m_regAdjArrays.reserveSpace(newEdges); // registry adds factor 2 in calculateArraySize

		for (auto it = edgesBegin; it != edgesEnd; ++it) {
			edge eG = *it;
			edge e = edgeMap[eG];
			if (nodeMap[eG->source()] == nullptr || nodeMap[eG->target()] == nullptr) {
				continue;
			}
			OGDF_ASSERT(e != nullptr);
			m_regEdgeArrays.keyAdded(e);
			m_regAdjArrays.keyAdded(e->adjSource());
			edgeInserted(cbData, eG, e);
			for (GraphObserver* obs : getObservers()) {
				obs->edgeAdded(e);
			}
		}
	}

#ifdef OGDF_HEAVY_DEBUG
	consistencyCheck();
#endif

	postInsert(cbData, newNodes, newEdges);
	return {newNodes, newEdges};
}

template<OGDF_NODE_ITER NI, OGDF_EDGE_FILTER EF, bool copyIDs, bool notifyObservers>
std::pair<int, int> Graph::insert(const NI& nodesBegin, const NI& nodesEnd, const EF& edgeFilter,
		NodeArray<node>& nodeMap, EdgeArray<edge>& edgeMap) {
	OGDF_ASSERT(nodeMap.valid());
	OGDF_ASSERT(edgeMap.valid());
	OGDF_ASSERT(nodeMap.graphOf() == edgeMap.graphOf());
	int newNodes = 0, newEdges = 0;
	void* cbData =
			preInsert(true, copyIDs, notifyObservers, true, nodeMap, edgeMap, &newNodes, &newEdges);
	if (nodesBegin == nodesEnd) {
		postInsert(cbData, newNodes, newEdges);
		return {newNodes, newEdges};
	}

	int guessedNodes = internal::guess_dist(nodesBegin, nodesEnd);
	if (guessedNodes > 0 && notifyObservers) {
		m_regNodeArrays.reserveSpace(guessedNodes);
	}

	for (auto it = nodesBegin; it != nodesEnd; ++it) {
		node vG = *it;
		if (copyIDs) {
			m_nodeIdCount = max(m_nodeIdCount, vG->index() + 1);
		}
		// nodeMap[vG] is overwritten if it is != nullptr
		node v = nodeMap[vG] = pureNewNode(copyIDs ? vG->index() : m_nodeIdCount++);
		newNodes++;
		if (notifyObservers) {
			m_regNodeArrays.keyAdded(v);
			nodeInserted(cbData, vG, v);
			for (GraphObserver* obs : getObservers()) {
				obs->nodeAdded(v);
			}
		}
	}

	for (auto it = nodesBegin; it != nodesEnd; ++it) {
		node vG = *it;
		node v = nodeMap[vG];
		for (adjEntry adjG : vG->adjEntries) {
			edge eG = adjG->m_edge;
			if (!edgeFilter(eG)) {
				continue;
			}
			// edgeMap[eG] is *not* overwritten if it is != nullptr
			edge e = edgeMap[eG];
			if (e == nullptr) {
				// add the first endpoint of the edge
				node twin = nodeMap[adjG->twinNode()];
				if (twin == nullptr) {
					// we can be sure that the other endpoint wasn't selected and
					// we thus cannot add this (selected) edge
					continue;
				}
				if (copyIDs) {
					m_edgeIdCount = max(m_edgeIdCount, eG->index() + 1);
				}
				if (adjG->isSource()) {
					e = edgeMap[eG] = pureNewEdge(v, twin, copyIDs ? eG->index() : m_edgeIdCount++);
					v->adjEntries.pushBack(e->m_adjSrc);
				} else {
					e = edgeMap[eG] = pureNewEdge(twin, v, copyIDs ? eG->index() : m_edgeIdCount++);
					v->adjEntries.pushBack(e->m_adjTgt);
				}
				newEdges++;
			} else {
				// complete the edge with its second endpoint
				adjEntry adj = adjG->isSource() ? e->adjSource() : e->adjTarget();
				// edgeMap[eG] might be an old entry that was already inserted into the list
				// so check whether adj is already in the list, indicated by having succ or pred,
				// or being the only entry in the list
				if (adj->succ() == nullptr && adj->pred() == nullptr && v->adjEntries.head() != adj) {
					v->adjEntries.pushBack(adj);
					// at this point, other edges might still be incomplete, so we cannot call observers
				}
			}
		}
	}

	// notify observers of added edges after all adjEntries are initialized
	if (notifyObservers) {
		m_regEdgeArrays.reserveSpace(newEdges);
		m_regAdjArrays.reserveSpace(newEdges); // registry adds factor 2 in calculateArraySize

		for (auto it = nodesBegin; it != nodesEnd; ++it) {
			node vG = *it;
			for (adjEntry adjG : vG->adjEntries) {
				if (!adjG->isSource()) {
					continue;
				}
				edge eG = adjG->m_edge;
				edge e = edgeMap[eG];
				// we will call Observers for *all* edgeMap entries
				if (e == nullptr) {
					continue;
				}
				m_regEdgeArrays.keyAdded(e);
				m_regAdjArrays.keyAdded(e->adjSource());
				edgeInserted(cbData, eG, e);
				for (GraphObserver* obs : getObservers()) {
					obs->edgeAdded(e);
				}
			}
		}
	}

#ifdef OGDF_HEAVY_DEBUG
	consistencyCheck();
#endif

	postInsert(cbData, newNodes, newEdges);
	return {newNodes, newEdges};
}

}
