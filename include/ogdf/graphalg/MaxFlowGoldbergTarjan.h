/** \file
 * \brief Declaration and implementation of Goldberg-Tarjan max-flow algorithm
 *        with global relabeling and gap relabeling heuristic
 *
 * \author Stephan Beyer, Hennes Hoffmann, Tilo Wiedera
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.txt in the root directory of the OGDF installation for details.
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
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/

#ifdef _MSC_VER
#pragma once
#endif

#ifndef OGDF_MAX_FLOW_GOLDBERG_TARJAN_H
#define OGDF_MAX_FLOW_GOLDBERG_TARJAN_H

#include <ogdf/module/MaxFlowModule.h>

//#define USE_GAP_RELABEL_HEURISTIC
#define USE_MAX_ACTIVE_LABEL
#ifdef USE_GAP_RELABEL_HEURISTIC
#  define GRH_STEPS	1	// gap relabel frequency: call gapRelabel() after GRH_STEPS relabel() operations (1 == off)
#endif
#define USE_PUSH_RELABEL_SECOND_STAGE
// world666 is much better without USE_PUSH_RELABEL_SECOND_STAGE

namespace ogdf {

//! Computes a max flow via Preflow-Push (global relabeling and gap relabeling heuristic).
/**
 * @ingroup ga-flow
 */
template<typename TCap>
class OGDF_EXPORT MaxFlowGoldbergTarjan : public MaxFlowModule<TCap>
{
	NodeArray<int> m_label;
	NodeArray<TCap> m_ex; // ex_f(v) values will be saved here to save runtime
#ifdef USE_MAX_ACTIVE_LABEL
	NodeArray< ListIterator<node> > m_activeLabelListPosition; // holds the iterator of every active node in the corresp. list of m_labeList
	Array< List<node> > m_activeLabelList; // array indexed by label, contains list of active nodes with that label
	int m_maxLabel; // the maximum label among all active nodes
#endif
#ifdef USE_GAP_RELABEL_HEURISTIC
	NodeArray< ListIterator<node> > m_labelListPosition; // holds the iterator of every node in the corresp. list of m_labeList
	Array< List<node> > m_labelList; // array indexed by label, contains list of nodes with that label
#endif

	mutable List<node> m_cutNodes;
	mutable List<edge> m_cutEdges;

	inline TCap getCap(const edge e) const {
		return e->target() == *this->m_s ? 0 : (*this->m_cap)[e];
	}

	inline bool isResidualEdge(const adjEntry adj) const
	{
		const edge e = adj->theEdge();
		if (adj->theNode() == e->source()) {
			return this->m_et->less((*this->m_flow)[e], getCap(e));
		}
		if (this->m_et->greater((*this->m_flow)[e], (TCap) 0)) {
			return true;
		}
		return false;
	}

	inline bool isAdmissible(const adjEntry adj) const
	{
		OGDF_ASSERT(adj);
		return (isResidualEdge(adj)
		     && m_label[adj->theNode()] == m_label[adj->twinNode()] + 1);
	}

	inline bool isActive(const node v) const
	{
		OGDF_ASSERT((v != *this->m_s && v != *this->m_t)
		  || (m_label[*this->m_s] == this->m_G->numberOfNodes() && m_label[*this->m_t] == 0));
		return (this->m_et->greater(m_ex[v], (TCap) 0)
		     && this->m_G->numberOfNodes() > m_label[v]
		     && m_label[v] > 0);
	}

#ifdef USE_MAX_ACTIVE_LABEL
	inline void setActive(const node v)
	{
		const int label = m_label[v];
		OGDF_ASSERT(0 < label);
		OGDF_ASSERT(label < this->m_G->numberOfNodes());
		OGDF_ASSERT(!m_activeLabelListPosition[v].valid());
		m_activeLabelListPosition[v] = m_activeLabelList[label].pushBack(v);
		if (label > m_maxLabel) {
			m_maxLabel = label;
		}
	}

	inline void findNewMaxLabel()
	{
		while (m_maxLabel > 0
		    && m_activeLabelList[m_maxLabel].empty()) {
			--m_maxLabel;
		}
	}

	inline void setInactive(const node v)
	{
		OGDF_ASSERT(m_activeLabelListPosition[v].valid());
		m_activeLabelList[m_label[v]].del(m_activeLabelListPosition[v]);
		m_activeLabelListPosition[v] = nullptr;
		findNewMaxLabel();
	}
#endif // USE_MAX_ACTIVE_LABEL

	// sets label of v, maintaining m_labelList (moves node v to the correct list in the array)
	inline void setLabel(const node v, int label)
	{
#ifdef USE_GAP_RELABEL_HEURISTIC
		if (m_labelListPosition[v].valid()) {
			m_labelList[m_label[v]].del(m_labelListPosition[v]); // delete node from old list using noted position
		}
		m_labelListPosition[v] = m_labelList[label].pushBack(v); // push node to new list and update iterator
#endif
#ifdef USE_MAX_ACTIVE_LABEL
		if (m_activeLabelListPosition[v].valid()) {
			OGDF_ASSERT(0 < m_label[v]);
			OGDF_ASSERT(m_label[v] < this->m_G->numberOfNodes());
			setInactive(v);
		}
		m_label[v] = label; // update label
		if (v != *this->m_s
		 && v != *this->m_t
		 && isActive(v)) {
			setActive(v);
		}
#else
		m_label[v] = label; // update label
#endif
	}

#ifdef USE_GAP_RELABEL_HEURISTIC
	void gapRelabel()
	{
#  ifdef USE_MAX_ACTIVE_LABEL
		// XXX: this is a test but it seems to work and it seems to be fast!
		const int n = m_maxLabel + 1;
#  else
		const int n = this->m_G->numberOfNodes();
#  endif
		for (int i = 1; i < n - 1; ++i) {
			if (m_labelList[i].empty()) {
				for (int j = i + 1; j < n; ++j) {
					while (!m_labelList[j].empty()) {
						setLabel(m_labelList[j].front(), this->m_G->numberOfNodes());
					}
				}
				break;
			}
		}
	}
#endif

	void push(const adjEntry adj)
	{
		const edge e = adj->theEdge();
		const node v = adj->theNode();
		if (v == e->source()) {
			const TCap value = min(m_ex[v], getCap(e) - (*this->m_flow)[e]);
			OGDF_ASSERT(this->m_et->geq(value, (TCap) 0));
			(*this->m_flow)[e] += value;
			m_ex[v] -= value;
			m_ex[adj->twinNode()] += value;
		} else {
			const TCap value = min(m_ex[v], (*this->m_flow)[adj]);
			OGDF_ASSERT(this->m_et->geq(value, (TCap) 0));
			(*this->m_flow)[adj] -= value;
			m_ex[v] -= value;
			m_ex[adj->twinNode()] += value;
		}
	}

	void globalRelabel()
	{
		// breadth-first search to relabel nodes with their respective distance to the sink in the residual graph
		const int n = this->m_G->numberOfNodes();
		NodeArray<int> dist(*this->m_G, n); // distance array
		List<node> queue; // reachable, not processed nodes
		dist[*this->m_t] = 0;
		queue.pushBack(*this->m_t);
		while (!queue.empty()) { // is there a node to check?
			node w = queue.popFrontRet();
			adjEntry adj;
			forall_adj(adj, w) {
				node x = adj->twinNode();
				if (isResidualEdge(adj->twin())
				 && dist[x] == n) { // not already seen
					dist[x] = dist[w] + 1; // set distance of node to sink
					queue.pushBack(x);
				}
			}
		}
		// set distance of unreachable nodes to "number of nodes" thus making them inactive
		node w;
		forall_nodes(w, *this->m_G) {
			setLabel(w, dist[w]);
		}
	}

	void relabel(const node v)
	{
		int minLabel = this->m_G->numberOfNodes() - 1;
		adjEntry adj;
		forall_adj(adj, v) {
			if (isResidualEdge(adj)) {
				const int label = m_label[adj->twinNode()];
				if (label < minLabel) {
					minLabel = label;
				}
			}
		}
		if (minLabel + 1 != m_label[v]) { // == can happen after global relabel
			setLabel(v, minLabel + 1);
		}
	}

	void relabelStage2(const node v)
	{
		int minLabel = this->m_G->numberOfNodes() - 1;
		adjEntry adj;
		forall_adj(adj, v) {
			if (isResidualEdge(adj)) {
				const int label = m_label[adj->twinNode()];
				if (label < minLabel) {
					minLabel = label;
				}
			}
		}
		OGDF_ASSERT(minLabel + 1 != m_label[v]);
		m_label[v] = minLabel + 1;
	}

public:
	// first stage: push excess towards sink
	TCap computeValue(const EdgeArray<TCap> &cap, const node &s, const node &t)
	{
		// TODO: init this stuff in the module?
		this->m_s = &s;
		this->m_t = &t;
		this->m_cap = &cap;
		this->m_flow->init(*this->m_G, (TCap) 0);

		m_label.init(*this->m_G);
		m_ex.init(*this->m_G, 0);
#ifdef USE_MAX_ACTIVE_LABEL
		m_activeLabelListPosition.init(*this->m_G, nullptr);
		m_activeLabelList.init(1, this->m_G->numberOfNodes() - 1);
		m_maxLabel = 0;
#endif
#ifdef USE_GAP_RELABEL_HEURISTIC
		m_labelListPosition.init(*this->m_G, nullptr);
		m_labelList.init(this->m_G->numberOfNodes() + 1);
#endif
		m_cutNodes.clear();

		// initialize residual graph for first preflow
		edge e;
		forall_edges(e, *this->m_G) {
			if (e->source() == *this->m_s
			&& e->target() != *this->m_s) { // ignore loops
				(*this->m_flow)[e] = getCap(e);
				m_ex[e->target()] += getCap(e); // "+" needed for the case of multigraphs
			}
		}

		if(*this->m_t == *this->m_s) {
			return (TCap) 0;
		}

		NodeArray<adjEntry> curr(*this->m_G);
		for (node v = this->m_G->firstNode(); v; v = v->succ()) {
			curr[v] = v->firstAdj();
		}

		globalRelabel(); // initialize distance labels

		int relCount = 0; // counts the relabel operations for the global relabeling heuristic
#ifdef USE_MAX_ACTIVE_LABEL
		while (m_maxLabel != 0) {
			OGDF_ASSERT(!m_activeLabelList[m_maxLabel].empty());
			const node v = m_activeLabelList[m_maxLabel].front();
			OGDF_ASSERT(m_maxLabel == m_label[v]);
			OGDF_ASSERT(m_activeLabelListPosition[v] == m_activeLabelList[m_maxLabel].begin());
#else
		List<node> active;
		forall_adj_edges(e, *this->m_s) {
			if (e->target() != *this->m_s) {
				active.pushBack(e->target());
			}
		}
		while (!active.empty()) {
			const node v = active.front();
#endif
			adjEntry &adj = curr[v];
			if (v == *this->m_s
			 || v == *this->m_t
			 || !isActive(v)) {
				// source, sink or not active: remove activity status
#ifdef USE_MAX_ACTIVE_LABEL
				setInactive(v);
#else
				active.popFront();
#endif
			} else {
				while (this->m_et->greater(m_ex[v], (TCap) 0)) {
					if (isAdmissible(adj)) {
						// push and adjacent node becomes active
#ifdef USE_MAX_ACTIVE_LABEL
						const node w = adj->twinNode();
						if (w != *this->m_s
						 && w != *this->m_t
						 && !isActive(w)) {
							// w will become active after push
							setActive(w);
						}
						push(adj);
						if (v != *this->m_s
						 && !isActive(v)) {
							setInactive(v);
						}
#else
						push(adj);
						active.pushBack(adj->twinNode());
#endif
					} else {
						if (adj != v->lastAdj()) {
							adj = adj->succ();
						} else { // end of adjacency list
							adj = v->firstAdj();
							relabel(v);
							++relCount;
#ifdef USE_GAP_RELABEL_HEURISTIC
							// only gapRelabel if we do not do a globalRelabel directly afterwards
							if (relCount != this->m_G->numberOfNodes()
#  if (GRH_STEPS > 1)
							 && relCount % GRH_STEPS == 0 // obey frequency of gap relabel heuristic
#  endif
							  ) {
								gapRelabel();
							}
#endif
							break;
						}
					}
				}
				if (relCount == this->m_G->numberOfNodes()) {
					relCount = 0;
					globalRelabel();
				}
			}
		}

		TCap result = 0;
		forall_adj_edges(e, *this->m_t) {
			if(e->target() == *this->m_t) {
				result += (*this->m_flow)[e];
			} else {
				result -= (*this->m_flow)[e];
			}
		}
		return result;
	}

	// second stage: push excess that has not reached the sink back towards source
	void computeFlowAfterValue()
	{
		List<node> active;
#ifdef USE_PUSH_RELABEL_SECOND_STAGE
		NodeArray<adjEntry> curr(*this->m_G);
		for (node v = this->m_G->firstNode(); v; v = v->succ()) {
			curr[v] = v->firstAdj();
			m_label[v] = 1;
			if (this->m_et->greater(m_ex[v], (TCap) 0) && v != *this->m_s && v != *this->m_t) {
				active.pushBack(v);
			}
		}
		if (active.empty()) {
			return;
		}

		m_label[*this->m_s] = 0;
		while (!active.empty()) {
			node v = active.front();
			if (v == *this->m_s
			 || v == *this->m_t
			 || !isActive(v)) {
				active.popFront();
			} else {
				adjEntry &adj = curr[v];
				if (isAdmissible(adj)) {
					push(adj);
					active.pushBack(adj->twinNode());
				} else {
					if (adj == v->lastAdj()) {
						// no admissible outgoing edge found -> relabel node!
						relabelStage2(v);
						adj = v->firstAdj();
#if 0
						// node is still active but move it to the end of the queue
						// (don't know if this is really necessary)
						active.popFront();
						active.pushBack(v);
#endif
					} else {
						adj = adj->succ();
					}
				}
			}
		}
#else // USE_PUSH_RELABEL_SECOND_STAGE
		m_ex[*this->m_s] = m_ex[*this->m_t] = 0;
		for (node v = this->m_G->firstNode(); v; v = v->succ()) {
			if (this->m_et->greater(m_ex[v], (T) 0)) {
				active.pushBack(v);
			}
		}
		while (!active.empty()) {
			const node v = active.popFrontRet();
			if (this->m_et->greater(m_ex[v], (T) 0) && v != *this->m_s && v != *this->m_t) {
				for (adjEntry adj = v->firstAdj(); adj; adj = adj->succ()) {
					const edge e = adj->theEdge();
					const node u = e->source();
					if (u != v) { // e is incoming edge
						if (this->m_et->greater(m_ex[v], (T) 0)
						 && isResidualEdge(adj)) {
							push(adj);
							if (u != *this->m_s) {
								active.pushFront(u);
							}
						}
					}
				}
			}
		}
#endif // USE_PUSH_RELABEL_SECOND_STAGE
	}
	using MaxFlowModule<TCap>::useEpsilonTest;
	using MaxFlowModule<TCap>::init;
	using MaxFlowModule<TCap>::computeFlow;
	using MaxFlowModule<TCap>::computeFlowAfterValue;
	using MaxFlowModule<TCap>::MaxFlowModule;
};

} // namespace ogdf

#endif // OGDF_MAX_FLOW_GOLDBERG_TARJAN_H
