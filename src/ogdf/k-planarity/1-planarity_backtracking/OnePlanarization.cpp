/** \file
 * \brief Implementation of OnePlanarization methods.
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
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/GraphSets.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/basic.h>
#include <ogdf/k-planarity/1-planarity_backtracking/EdgePairPartition.h>
#include <ogdf/k-planarity/1-planarity_backtracking/OnePlanarization.h>

#include <set>

using namespace ogdf;
using namespace oneplan_backtracking;

void OnePlanarization::init(const EdgePairPartition* ep) {
	m_epp = ep;
	GraphCopy::init(m_epp->graph());
	m_saturatedVertices.init(*this);
	m_crossingVertices.init(*this);
	m_crossingEdges.init(*this);
	m_kiteEdges.init(*this);
	m_freeEdges.init(*this);
	m_remainingEdges.init(*this);

	for (const EdgePair& crossPair : m_epp->crossingEdgePairs()) {
		edge c1 = copy(crossPair.second());
		node crossingVertex = insertCrossing(c1, copy(crossPair.first()), true)->source();
		m_crossingVertices.insert(crossingVertex);
		for (adjEntry adj : crossingVertex->adjEntries) {
			OGDF_ASSERT(original(adj->theEdge()) == crossPair.first()
					|| original(adj->theEdge()) == crossPair.second());
			m_crossingEdges.insert(adj->theEdge());
		}
	}

	for (const VertexPair& vp : m_epp->kiteEdges()) {
		edge e = EdgePairPartition::getEdgeBetween(vp.first(), vp.second());
		if (e == nullptr) {
			m_kiteEdges.insert(newEdge(copy(vp.first()), copy(vp.second())));
		}
	}

	for (edge e : m_epp->graph().edges) {
		if (m_epp->isFree(e)) {
			m_freeEdges.insert(copy(e));
		} else if (!m_epp->isCrossed(e)) {
			m_remainingEdges.insert(copy(e));
		}
	}

	for (edge e : edges) {
		if (!m_freeEdges.isMember(e)) {
			m_saturatedVertices.insert(e->source());
			m_saturatedVertices.insert(e->target());
		}
	}

#ifdef OGDF_DEBUG
	OGDF_ASSERT(numberOfEdges()
			== m_crossingEdges.size() + m_kiteEdges.size() + m_freeEdges.size()
					+ m_remainingEdges.size());
	OGDF_ASSERT(m_crossingVertices.size() == (int)m_epp->crossingEdgePairs().size());
	OGDF_ASSERT(numberOfEdges()
			== m_epp->graph().numberOfEdges() + 2 * m_crossingVertices.size() + m_kiteEdges.size());
	OGDF_ASSERT(numberOfNodes() == m_epp->graph().numberOfNodes() + m_crossingVertices.size());
	for (node cv : m_crossingVertices) {
		for (adjEntry adj : cv->adjEntries) {
			OGDF_ASSERT(m_crossingEdges.isMember(adj->theEdge()));
		}
	}
#endif
}
