/** \file
 * \brief Implementation of EdgePairPartition methods.
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

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/GraphSets.h>
#include <ogdf/basic/basic.h>
#include <ogdf/k-planarity/1-planarity_backtracking/EdgePairPartition.h>

#include <deque>
#include <initializer_list>
#include <list>
#include <memory>
#include <set>
#include <stack>
#include <vector>

using namespace std;
using namespace ogdf;
using namespace oneplan_backtracking;

void EdgePairPartition::undoTransaction() {
	for (const EdgePair& ep : m_undoInformation.top().m_crossings) {
		m_freeEdgePairs.insert(ep);
		m_crossingEdgePairs.erase(ep);
		m_crossedEdges.remove(ep.first());
		m_crossedEdges.remove(ep.second());
	}

	for (const VertexPair& vp : m_undoInformation.top().m_kiteEdges) {
		m_kiteEdges.erase(vp);
	}

	for (const EdgePair& ep : m_undoInformation.top().m_nonCrossingPairs) {
		m_freeEdgePairs.insert(ep);
	}
	m_undoInformation.pop();

	if (!m_undoInformation.empty() && !m_undoInformation.top().m_todoCrossings.empty()
			&& m_undoInformation.top().m_previousCrossing) {
		setNonCrossing(*m_undoInformation.top().m_previousCrossing);
	}
}

void EdgePairPartition::crossEdgePair(const EdgePair& pair) {
	OGDF_ASSERT(m_freeEdgePairs.count(pair) > 0);
	OGDF_ASSERT(!m_crossedEdges.contains(pair.first()));
	OGDF_ASSERT(!m_crossedEdges.contains(pair.second()));
	OGDF_ASSERT(m_kiteEdges.count({pair.first()->source(), pair.first()->target()}) == 0);
	OGDF_ASSERT(m_kiteEdges.count({pair.second()->source(), pair.second()->target()}) == 0);
	m_freeEdgePairs.erase(pair);
	m_crossingEdgePairs.insert(pair);
	m_crossedEdges.insert(pair.first());
	m_crossedEdges.insert(pair.second());
	m_undoInformation.top().m_crossings.insert(pair);
	if (m_mode == OneplanMode::NIC) {
		vector<node> nl = {pair.first()->source(), pair.first()->target(), pair.second()->source(),
				pair.second()->target()};
		set<VertexPair> distinctPairs;
		for (node v1 : nl) {
			for (node v2 : nl) {
				if (v1 != v2) {
					distinctPairs.insert({v1, v2});
				}
			}
		}
		for (const VertexPair& vp : distinctPairs) {
			for (adjEntry adj1 : vp.first()->adjEntries) {
				for (adjEntry adj2 : vp.second()->adjEntries) {
					if (adj1->theEdge() == adj2->theEdge()) {
						continue;
					}
					EdgePair ep = {adj1->theEdge(), adj2->theEdge()};
					if (isFree(ep)) {
						setNonCrossing(ep);
					}
				}
			}
		}
	} else if (m_mode == OneplanMode::IC) {
		vector<node> nl = {pair.first()->source(), pair.first()->target(), pair.second()->source(),
				pair.second()->target()};
		for (node v : nl) {
			for (adjEntry adj : v->adjEntries) {
				if (isFree(adj->theEdge())) {
					addKiteEdge(adj->theEdge());
				}
			}
		}
	}
	computeKiteEdges(pair);

#ifdef OGDF_DEBUG
	for (edge e : m_graph.edges) {
		int nCrossings = 0;
		for (edge e2 : m_graph.edges) {
			if (e == e2) {
				continue;
			}
			if (m_crossingEdgePairs.count({e, e2}) > 0) {
				++nCrossings;
			}
		}
		OGDF_ASSERT(nCrossings <= 1);
	}
#endif
}

void EdgePairPartition::computeKiteEdges(const EdgePair& crossedPair) {
	VertexPair p1(crossedPair.first()->source(), crossedPair.second()->source()),
			p2(crossedPair.first()->source(), crossedPair.second()->target()),
			p3(crossedPair.first()->target(), crossedPair.second()->source()),
			p4(crossedPair.first()->target(), crossedPair.second()->target());

	for (VertexPair p : {p1, p2, p3, p4}) {
		ogdf::edge e = getEdgeBetween(p.first(), p.second());
		if (e) {
			addKiteEdge(e);
		} else {
			addKiteEdge(p);
		}

		// degree-2 paths between the endpoints can alo be realized without crossings -> mark as uncrossable
		std::stack<ogdf::node> dfsStack({p.first()});
		ogdf::NodeArray<bool> visited(m_graph, false);
		visited[p.first()] = true;

		std::vector<ogdf::edge> pathEnds;
		while (!dfsStack.empty()) {
			ogdf::node n = dfsStack.top();
			dfsStack.pop();

			visited[n] = true;
			if (n != p.first() && n->degree() != 2) {
				continue;
			}
			for (ogdf::adjEntry adj : n->adjEntries) {
				ogdf::node t = adj->twinNode();
				if (t == p.second()) {
					pathEnds.push_back(adj->theEdge());
				} else if (!visited[t]) {
					dfsStack.push(t);
				}
				visited[t] = true;
			}
		}

		for (ogdf::edge s : pathEnds) {
			addKiteEdge(s);
			ogdf::node n = s->opposite(p.second());
			while (n != p.first()) {
				OGDF_ASSERT(n->degree() == 2);
				s = n->firstAdj()->theEdge() == s ? n->lastAdj()->theEdge()
												  : n->firstAdj()->theEdge();
				addKiteEdge(s);
				n = s->opposite(n);
			}
		}
	}
}
