/** \file
 * \brief Implementation of Graph class
 *
 * \author Carsten Gutwenger
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

#include <ogdf/basic/Array.h>
#include <ogdf/basic/ArrayBuffer.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/simple_graph_alg.h>

#include <functional>
#include <ostream>
#include <utility>

namespace ogdf {

Graph::Graph()
	: m_regNodeArrays(this, &m_nodeIdCount)
	, m_regEdgeArrays(this, &m_edgeIdCount)
	, m_regAdjArrays(this, &m_edgeIdCount) {
	m_nodeIdCount = m_edgeIdCount = 0;
}

Graph::Graph(const Graph& G)
	: m_regNodeArrays(this, &m_nodeIdCount)
	, m_regEdgeArrays(this, &m_edgeIdCount)
	, m_regAdjArrays(this, &m_edgeIdCount) {
	m_nodeIdCount = m_edgeIdCount = 0;
	insert(G);
}

Graph::~Graph() {
	clearObservers();
	restoreAllEdges();

	// this is only necessary because GraphObjectContainer simply deallocs its memory without calling destructors
	for (node v = nodes.head(); v; v = v->succ()) {
		v->adjEntries.~GraphObjectContainer<AdjElement>();
	}
}

void Graph::clear() {
	restoreAllEdges();

	// tell all structures to clear their graph-initialized data
	for (GraphObserver* obs : getObservers()) {
		obs->cleared();
	}

	m_regNodeArrays.keysCleared();
	m_regEdgeArrays.keysCleared();
	m_regAdjArrays.keysCleared();

	for (node v = nodes.head(); v; v = v->succ()) {
		v->adjEntries.~GraphObjectContainer<AdjElement>();
	}
	nodes.clear();
	edges.clear();

	m_nodeIdCount = m_edgeIdCount = 0;

#ifdef OGDF_HEAVY_DEBUG
	consistencyCheck();
#endif
}

Graph& Graph::operator=(const Graph& G) {
	clear();
	insert(G);
	return *this;
}

void Graph::move(edge e, adjEntry adjSrc, Direction dirSrc, adjEntry adjTgt, Direction dirTgt) {
	OGDF_ASSERT(e->graphOf() == this);
	OGDF_ASSERT(adjSrc->graphOf() == this);
	OGDF_ASSERT(adjTgt->graphOf() == this);
	OGDF_ASSERT(adjSrc != e->m_adjSrc);
	OGDF_ASSERT(adjSrc != e->m_adjTgt);
	OGDF_ASSERT(adjTgt != e->m_adjSrc);
	OGDF_ASSERT(adjTgt != e->m_adjTgt);

	node v = adjSrc->m_node, w = adjTgt->m_node;
	adjEntry adj1 = e->m_adjSrc, adj2 = e->m_adjTgt;
	e->m_src->adjEntries.move(adj1, v->adjEntries, adjSrc, dirSrc);
	e->m_tgt->adjEntries.move(adj2, w->adjEntries, adjTgt, dirTgt);

	e->m_src->m_outdeg--;
	e->m_tgt->m_indeg--;

	adj1->m_node = e->m_src = v;
	adj2->m_node = e->m_tgt = w;

	v->m_outdeg++;
	w->m_indeg++;
}

void Graph::moveTarget(edge e, node v) {
	OGDF_ASSERT(e->graphOf() == this);
	OGDF_ASSERT(v->graphOf() == this);

	adjEntry adj = e->m_adjTgt;
	e->m_tgt->adjEntries.move(adj, v->adjEntries);

	e->m_tgt->m_indeg--;
	adj->m_node = e->m_tgt = v;
	v->m_indeg++;
}

void Graph::moveTarget(edge e, adjEntry adjTgt, Direction dir) {
	node v = adjTgt->theNode();

	OGDF_ASSERT(e->graphOf() == this);
	OGDF_ASSERT(v->graphOf() == this);

	adjEntry adj = e->m_adjTgt;
	e->m_tgt->adjEntries.move(adj, v->adjEntries, adjTgt, dir);

	e->m_tgt->m_indeg--;
	adj->m_node = e->m_tgt = v;
	v->m_indeg++;
}

// By Leipert
void Graph::moveSource(edge e, node v) {
	OGDF_ASSERT(e->graphOf() == this);
	OGDF_ASSERT(v->graphOf() == this);

	adjEntry adj = e->m_adjSrc;
	e->m_src->adjEntries.move(adj, v->adjEntries);

	e->m_src->m_outdeg--;
	adj->m_node = e->m_src = v;
	v->m_outdeg++;
}

void Graph::moveSource(edge e, adjEntry adjSrc, Direction dir) {
	node v = adjSrc->theNode();

	OGDF_ASSERT(e->graphOf() == this);
	OGDF_ASSERT(v->graphOf() == this);

	adjEntry adj = e->m_adjSrc;
	e->m_src->adjEntries.move(adj, v->adjEntries, adjSrc, dir);

	e->m_src->m_outdeg--;
	adj->m_node = e->m_src = v;
	v->m_outdeg++;
}

edge Graph::split(edge e) {
	OGDF_ASSERT(e != nullptr);
	OGDF_ASSERT(e->graphOf() == this);

	node /* src = e->source(), */ tgt = e->target();
	adjEntry eadjSrc = e->adjSource(), eadjTgt = e->adjTarget();

	// before:
	// src (eadjSrc) -- e -> (eadjTgt) tgt
	// after:
	// src (eadjSrc) -- e -> (uadjTgt) u (uadjSrc) -- e2 -> (eadjTgt) tgt

	node u = newNode();
	u->m_indeg = u->m_outdeg = 1;

	adjEntry uadjTgt = new AdjElement(u);
	uadjTgt->m_twin = eadjSrc;
	eadjSrc->m_twin = uadjTgt;
	u->adjEntries.pushBack(uadjTgt);

	adjEntry uadjSrc = new AdjElement(u);
	uadjSrc->m_twin = eadjTgt;
	eadjTgt->m_twin = uadjSrc;
	u->adjEntries.pushBack(uadjSrc);

	edge e2 = new EdgeElement(u, tgt, uadjSrc, eadjTgt, m_edgeIdCount);
	edges.pushBack(e2);

	// update indices
	uadjSrc->m_id = m_edgeIdCount << 1;
	uadjTgt->m_id = eadjTgt->m_id; // must be read before being overwritten
	eadjTgt->m_id = (uadjSrc->m_id) | 1;
	m_edgeIdCount++;

	// update edge
	uadjTgt->m_edge = e;
	uadjSrc->m_edge = e2;
	eadjTgt->m_edge = e2;

	e->m_tgt = u;
	e->m_adjTgt = uadjTgt;

	m_regEdgeArrays.keyAdded(e2); // note: registry observers won't see copied entry
	m_regAdjArrays.keyAdded(e2->adjSource());
	m_regAdjArrays.keyAdded(e2->adjTarget());

	// copy array entries from the original adjEntries to the new ones
	m_regAdjArrays.copyArrayEntries(eadjTgt->m_id, uadjTgt->m_id);
	m_regAdjArrays.copyArrayEntries(uadjSrc->m_id, eadjSrc->m_id);

	for (GraphObserver* obs : getObservers()) {
		obs->edgeAdded(e2);
	}
	return e2;
}

void Graph::unsplit(node u) {
	edge eIn = u->firstAdj()->theEdge();
	edge eOut = u->lastAdj()->theEdge();

	if (eIn->target() != u) {
		std::swap(eIn, eOut);
	}

	unsplit(eIn, eOut);
}

void Graph::unsplit(edge eIn, edge eOut) {
	node u = eIn->target();

	// u must be a node with exactly one incoming edge eIn and one outgoing
	// edge eOut
	OGDF_ASSERT(u->graphOf() == this);
	OGDF_ASSERT(u->indeg() == 1);
	OGDF_ASSERT(u->outdeg() == 1);
	OGDF_ASSERT(eOut->source() == u);

	// none of them is a self-loop!
	OGDF_ASSERT(!eIn->isSelfLoop());
	OGDF_ASSERT(!eOut->isSelfLoop());

	// notify all registered observers while everything is still valid
	m_regAdjArrays.keyRemoved(eOut->adjSource());
	m_regAdjArrays.keyRemoved(eOut->adjTarget());
	m_regEdgeArrays.keyRemoved(eOut);
	for (GraphObserver* obs : getObservers()) {
		obs->edgeDeleted(eOut);
	}
	m_regNodeArrays.keyRemoved(u);
	for (GraphObserver* obs : getObservers()) {
		obs->nodeDeleted(u);
	}

	// before:
	// src (adjSrc) -- eIn -> u -- eOut -> (adjTgt) tgt
	// after:
	// src (adjSrc) -- eIn -> (adjTgt) tgt

	// we reuse these adjacency entries
	adjEntry adjSrc = eIn->m_adjSrc;
	adjEntry adjTgt = eOut->m_adjTgt;

	// adapt adjacency entry index to hold invariant
	m_regAdjArrays.swapArrayEntries(eIn->m_adjTgt->m_id, eOut->m_adjTgt->m_id);
	eOut->m_adjTgt->m_id = eIn->m_adjTgt->m_id;
	eIn->m_adjTgt = eOut->m_adjTgt;
	eIn->m_tgt = eOut->m_tgt;

	adjSrc->m_twin = adjTgt;
	adjTgt->m_twin = adjSrc;
	adjTgt->m_edge = eIn;

	// remove structures that are no longer used
	edges.del(eOut);
	nodes.del(u);
}

void Graph::delNode(node v) {
	OGDF_ASSERT(v != nullptr);
	OGDF_ASSERT(v->graphOf() == this);

	// delete all edges first, notifying the respective observers
	for (AdjElement* adj = v->adjEntries.head(); adj != nullptr; adj = v->adjEntries.head()) {
		delEdge(adj->m_edge);
	}

	m_regNodeArrays.keyRemoved(v);
	for (GraphObserver* obs : getObservers()) {
		obs->nodeDeleted(v);
	}

	nodes.del(v);
}

void Graph::delEdge(edge e) {
	OGDF_ASSERT(e != nullptr);
	OGDF_ASSERT(e->graphOf() == this);

	m_regAdjArrays.keyRemoved(e->adjSource());
	m_regAdjArrays.keyRemoved(e->adjTarget());
	m_regEdgeArrays.keyRemoved(e);
	for (GraphObserver* obs : getObservers()) {
		obs->edgeDeleted(e);
	}

	node src = e->m_src, tgt = e->m_tgt;

	src->adjEntries.del(e->m_adjSrc);
	src->m_outdeg--;
	tgt->adjEntries.del(e->m_adjTgt);
	tgt->m_indeg--;

	edges.del(e);
}

void Graph::reverseEdge(edge e) {
	OGDF_ASSERT(e != nullptr);
	OGDF_ASSERT(e->graphOf() == this);
	node &src = e->m_src, &tgt = e->m_tgt;

	std::swap(src, tgt);
	std::swap(e->m_adjSrc, e->m_adjTgt);
	src->m_outdeg++;
	src->m_indeg--;
	tgt->m_outdeg--;
	tgt->m_indeg++;
}

void Graph::reverseAllEdges() {
	for (edge e = edges.head(); e; e = e->succ()) {
		reverseEdge(e);
	}

#ifdef OGDF_HEAVY_DEBUG
	consistencyCheck();
#endif
}

void Graph::reverseAdjEdges() {
	for (node v : nodes) {
		reverseAdjEdges(v);
	}
}

node Graph::chooseNode(std::function<bool(node)> includeNode, bool isFastTest) const {
	return *chooseIteratorFrom<internal::GraphObjectContainer<NodeElement>, node>(
			const_cast<internal::GraphObjectContainer<NodeElement>&>(nodes),
			[&](const node& v) { return includeNode(v); }, isFastTest);
}

edge Graph::chooseEdge(std::function<bool(edge)> includeEdge, bool isFastTest) const {
	return *chooseIteratorFrom<internal::GraphObjectContainer<EdgeElement>, edge>(
			const_cast<internal::GraphObjectContainer<EdgeElement>&>(edges),
			[&](const edge& e) { return includeEdge(e); }, isFastTest);
}

edge Graph::searchEdge(node v, node w, bool directed) const {
	OGDF_ASSERT(v != nullptr);
	OGDF_ASSERT(v->graphOf() == this);
	OGDF_ASSERT(w != nullptr);
	OGDF_ASSERT(w->graphOf() == this);

	bool swapped = false;
	if (w->degree() < v->degree()) {
		std::swap(v, w);
		swapped = true;
	}

	for (adjEntry adj : v->adjEntries) {
		if (adj->twinNode() == w && (!directed || swapped != adj->isSource())) {
			return adj->theEdge();
		}
	}
	return nullptr;
}

void Graph::restoreAllEdges() {
	while (!m_hiddenEdgeSets.empty()) {
		HiddenEdgeSet* set = m_hiddenEdgeSets.popFrontRet();
		set->restore();
		set->m_graph = nullptr;
	}
}

int Graph::genus() const {
	if (empty()) {
		return 0;
	}

	int nIsolated = 0;
	for (node v : nodes) {
		if (v->degree() == 0) {
			++nIsolated;
		}
	}

	NodeArray<int> component(*this);
	int nCC = connectedComponents(*this, component);

	AdjEntryArray<bool> visited(*this, false);
	int nFaceCycles = 0;

	for (node v : nodes) {
		for (adjEntry adj1 : v->adjEntries) {
			if (visited[adj1]) {
				continue;
			}

			adjEntry adj = adj1;
			do {
				visited[adj] = true;
				adj = adj->faceCycleSucc();
			} while (adj != adj1);

			++nFaceCycles;
		}
	}

	return (numberOfEdges() - numberOfNodes() - nIsolated - nFaceCycles + 2 * nCC) / 2;
}


#ifdef OGDF_DEBUG
void Graph::consistencyCheck() const {
	int n = 0;
	for (node v : nodes) {
		OGDF_ASSERT(v->graphOf() == this);

		n++;
		int in = 0, out = 0;

		for (adjEntry adj : v->adjEntries) {
			edge e = adj->m_edge;
			OGDF_ASSERT(adj->m_twin->m_edge == e);

			if (e->m_adjSrc == adj) {
				out++;
			} else {
				OGDF_ASSERT(e->m_adjTgt == adj);
				in++;
			}

			OGDF_ASSERT(adj->m_node == v);
			OGDF_ASSERT(adj->graphOf() == this);
		}

		OGDF_ASSERT(v->m_indeg == in);
		OGDF_ASSERT(v->m_outdeg == out);
	}

	OGDF_ASSERT(n == nodes.size());

	int m = 0;
	for (edge e : edges) {
		m++;
		OGDF_ASSERT(e->graphOf() == this);
		OGDF_ASSERT(e->m_adjSrc != e->m_adjTgt);
		OGDF_ASSERT(e->m_adjSrc->m_edge == e);
		OGDF_ASSERT(e->m_adjTgt->m_edge == e);
		OGDF_ASSERT(e->m_adjSrc->m_node == e->m_src);
		OGDF_ASSERT(e->m_adjTgt->m_node == e->m_tgt);
	}

	OGDF_ASSERT(m == edges.size());
}
#endif


void Graph::resetEdgeIdCount(int maxId) {
	m_edgeIdCount = maxId + 1;

#ifdef OGDF_HEAVY_DEBUG
	for (edge e : edges) {
		// if there is an edge with higer index than maxId, we cannot
		// set the edge id count to maxId+1
		OGDF_ASSERT(e->index() <= maxId);
	}
#endif
}

void Graph::resetNodeIdCount(int maxId) {
	m_nodeIdCount = maxId + 1;

#ifdef OGDF_HEAVY_DEBUG
	for (node n : nodes) {
		// if there is a node with higer index than maxId, we cannot
		// set the node id count to maxId+1
		OGDF_ASSERT(n->index() <= maxId);
	}
#endif
}

node Graph::splitNode(adjEntry adjStartLeft, adjEntry adjStartRight) {
	OGDF_ASSERT(adjStartLeft != nullptr);
	OGDF_ASSERT(adjStartRight != nullptr);
	OGDF_ASSERT(adjStartLeft->graphOf() == this);
	OGDF_ASSERT(adjStartRight->graphOf() == this);
	OGDF_ASSERT(adjStartLeft->theNode() == adjStartRight->theNode());

	node w = newNode();

	adjEntry adj, adjSucc;
	for (adj = adjStartRight; adj != adjStartLeft; adj = adjSucc) {
		adjSucc = adj->cyclicSucc();
		moveAdj(adj, w);
	}

	if (adjStartLeft == adjStartRight) {
		newEdge(adjStartLeft->cyclicPred(), w);
	} else {
		newEdge(adjStartLeft, adjStartRight, Direction::before);
	}

	return w;
}

node Graph::contract(edge e, bool keepSelfLoops) {
	adjEntry adjSrc = e->adjSource();
	adjEntry adjTgt = e->adjTarget();
	node v = e->source();

	adjEntry adjNext;
	for (adjEntry adj = adjTgt->cyclicSucc(); adj != adjTgt; adj = adjNext) {
		adjNext = adj->cyclicSucc();
		if (keepSelfLoops || adj->twinNode() != v) {
			if (adj->isSource()) {
				moveSource(adj->theEdge(), adjSrc, Direction::before);
			} else {
				moveTarget(adj->theEdge(), adjSrc, Direction::before);
			}
		}
	}

	delNode(adjTgt->theNode());

	return v;
}

void Graph::moveAdj(adjEntry adj, node w) {
	node v = adj->m_node;

	v->adjEntries.move(adj, w->adjEntries);
	adj->m_node = w;

	edge e = adj->m_edge;
	if (adj->isSource()) {
		--v->m_outdeg;
		e->m_src = w;
		++w->m_outdeg;
	} else {
		--v->m_indeg;
		e->m_tgt = w;
		++w->m_indeg;
	}
}

std::ostream& operator<<(std::ostream& os, ogdf::node v) {
	if (v) {
		os << v->index();
	} else {
		os << "nil";
	}
	return os;
}

std::ostream& operator<<(std::ostream& os, ogdf::edge e) {
	if (e) {
		os << "(" << e->source() << "," << e->target() << ")";
	} else {
		os << "nil";
	}
	return os;
}

std::ostream& operator<<(std::ostream& os, ogdf::adjEntry adj) {
	if (adj) {
		ogdf::edge e = adj->theEdge();
		if (adj == e->adjSource()) {
			os << e->source() << "->" << e->target();
		} else {
			os << e->target() << "->" << e->source();
		}
	} else {
		os << "nil";
	}
	return os;
}

Graph::CCsInfo::CCsInfo(const Graph& G)
	: m_graph(&G), m_nodes(G.numberOfNodes()), m_edges(G.numberOfEdges()) {
	NodeArray<int> component(G, -1);

	ArrayBuffer<node> S;
	SListPure<int> startNodes;
	SListPure<int> startEdges;
	int nComponent = 0, n = 0, m = 0;

	for (node v : G.nodes) {
		if (component[v] != -1) {
			continue;
		}

		S.push(v);
		component[v] = nComponent;

		while (!S.empty()) {
			node w = S.popRet();
			m_nodes[n++] = w;

			for (adjEntry adj : w->adjEntries) {
				if ((adj->index() & 1) == 0) {
					m_edges[m++] = adj->theEdge();
				}
				node x = adj->twinNode();
				if (component[x] == -1) {
					component[x] = nComponent;
					S.push(x);
				}
			}
		}

		++nComponent;
		startNodes.pushBack(n);
		startEdges.pushBack(m);
	}

	m_startNode.init(nComponent + 1);
	m_startNode[0] = 0;

	int i = 1;
	for (int j : startNodes) {
		m_startNode[i++] = j;
	}


	m_startEdge.init(nComponent + 1);
	m_startEdge[0] = 0;

	i = 1;
	for (int j : startEdges) {
		m_startEdge[i++] = j;
	}

	m_numCC = nComponent;
}

void Graph::HiddenEdgeSet::hide(edge e) {
	OGDF_ASSERT(m_graph == e->graphOf());
	OGDF_ASSERT(!e->m_hidden);

	node src = e->m_src, tgt = e->m_tgt;

	src->adjEntries.delPure(e->m_adjSrc);
	src->m_outdeg--;
	tgt->adjEntries.delPure(e->m_adjTgt);
	tgt->m_indeg--;

	m_graph->edges.move(e, m_edges);
#ifdef OGDF_DEBUG
	e->m_hidden = true;
#endif
}

void Graph::HiddenEdgeSet::restore(edge e) {
	OGDF_ASSERT(m_graph == e->graphOf());
	OGDF_ASSERT(e->m_hidden);
	OGDF_ASSERT(!m_edges.empty());

	node v = e->m_src;
	v->adjEntries.pushBack(e->m_adjSrc);
	++v->m_outdeg;

	node w = e->m_tgt;
	w->adjEntries.pushBack(e->m_adjTgt);
	++w->m_indeg;

	m_edges.move(e, m_graph->edges);
#ifdef OGDF_DEBUG
	e->m_hidden = false;
#endif
}

void Graph::HiddenEdgeSet::restore() {
	OGDF_ASSERT(m_graph != nullptr);

	while (!m_edges.empty()) {
		restore(m_edges.head());
	}
}

int Graph::HiddenEdgeSet::size() { return m_edges.size(); }

bool Graph::HiddenEdgeSet::empty() { return m_edges.empty(); }

internal::GraphList<EdgeElement>::iterator Graph::HiddenEdgeSet::begin() { return m_edges.begin(); }

internal::GraphList<EdgeElement>::iterator Graph::HiddenEdgeSet::end() { return m_edges.end(); }

std::ostream& operator<<(std::ostream& os, const Graph::EdgeType& et) {
	switch (et) {
	case Graph::EdgeType::association:
		os << "association";
		break;
	case Graph::EdgeType::generalization:
		os << "generalization";
		break;
	case Graph::EdgeType::dependency:
		os << "dependency";
		break;
	}
	return os;
}

#ifndef DOXYGEN_IGNORE
namespace internal {
GraphNodeRegistry::iterator begin(const GraphNodeRegistry& reg) {
	return reg.graphOf()->nodes.begin();
}

GraphNodeRegistry::iterator end(const GraphNodeRegistry& reg) { return reg.graphOf()->nodes.end(); }

GraphEdgeRegistry::iterator begin(const GraphEdgeRegistry& reg) {
	return reg.graphOf()->edges.begin();
}

GraphEdgeRegistry::iterator end(const GraphEdgeRegistry& reg) { return reg.graphOf()->edges.end(); }

GraphAdjRegistry::iterator begin(const GraphAdjRegistry& reg) {
	return GraphAdjIterator(reg.graphOf()).begin();
}

GraphAdjRegistry::iterator end(const GraphAdjRegistry& reg) {
	return GraphAdjIterator(reg.graphOf());
}
}
#endif

GraphAdjIterator::GraphAdjIterator(Graph* graph, adjEntry entry)
	: m_pGraph(graph), m_entry(entry) { }

GraphAdjIterator GraphAdjIterator::begin() {
	node v = m_pGraph->firstNode();
	while (v != nullptr && v->firstAdj() == nullptr) {
		v = v->succ();
	}
	return {m_pGraph, (v != nullptr) ? v->firstAdj() : nullptr};
}

void GraphAdjIterator::next() {
	OGDF_ASSERT(m_entry != nullptr);
	if (m_entry->succ() != nullptr) {
		m_entry = m_entry->succ();
	} else {
		node v = m_entry->theNode()->succ();
		while (v != nullptr && v->firstAdj() == nullptr) {
			v = v->succ();
		}
		m_entry = (v != nullptr) ? v->firstAdj() : nullptr;
	}
}

void GraphAdjIterator::prev() {
	OGDF_ASSERT(m_entry != nullptr);
	if (m_entry->pred() != nullptr) {
		m_entry = m_entry->pred();
	} else {
		node v = m_entry->theNode()->pred();
		while (v != nullptr && v->lastAdj() == nullptr) {
			v = v->pred();
		}
		m_entry = (v != nullptr) ? v->lastAdj() : nullptr;
	}
}

bool filter_any_edge(edge e) { return true; }

bool filter_any_node(node n) { return true; }

}
