/** \file
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
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/basic.h>
#include <ogdf/cluster/sync_plan/basic/OverlappingGraphCopies.h>

namespace ogdf {

const Graph& OverlappingGraphCopy::original() const { return *m_pOGC->m_G; }

node OverlappingGraphCopy::copy(node v) const {
	if (m_pOGC->m_node_copies.contains(v, this)) {
		node cv = m_pOGC->m_node_copies(v, this);
		OGDF_ASSERT(cv == nullptr || cv->graphOf() == this);
		return cv;
	} else {
		return nullptr;
	}
}

edge OverlappingGraphCopy::copy(edge e) const {
	if (m_pOGC->m_edge_copies.contains(e, this)) {
		edge ce = m_pOGC->m_edge_copies(e, this);
		OGDF_ASSERT(ce == nullptr || ce->graphOf() == this);
		return ce;
	} else {
		return nullptr;
	}
}

node OverlappingGraphCopy::newNode(node vOrig) {
	OGDF_ASSERT(vOrig != nullptr);
	OGDF_ASSERT(vOrig->graphOf() == &original());
	node v = Graph::newNode();
	m_vOrig[v] = vOrig;
	OGDF_ASSERT(m_pOGC->m_node_copies(vOrig, this) == nullptr);
	m_pOGC->m_node_copies(vOrig, this) = v;
	return v;
}

edge OverlappingGraphCopy::newEdge(edge eOrig) {
	OGDF_ASSERT(eOrig != nullptr);
	OGDF_ASSERT(eOrig->graphOf() == &original());
	node s = m_pOGC->m_node_copies(eOrig->source(), this);
	node t = m_pOGC->m_node_copies(eOrig->target(), this);
	edge e = Graph::newEdge(s, t);
	m_eOrig[e] = eOrig;
	OGDF_ASSERT(m_pOGC->m_edge_copies(eOrig, this) == nullptr);
	m_pOGC->m_edge_copies(eOrig, this) = e;
	return e;
}

void OverlappingGraphCopy::delNode(node v) {
	if (!isDummy(v)) {
		m_pOGC->m_node_copies.remove(original(v), this);
	}
	Graph::delNode(v);
}

void OverlappingGraphCopy::delEdge(edge e) {
	if (!isDummy(e)) {
		m_pOGC->m_edge_copies.remove(original(e), this);
	}
	Graph::delEdge(e);
}

void OverlappingGraphCopy::unmap() {
	if (m_pOGC == nullptr) {
		return;
	}
	for (edge e : edges) {
		if (!isDummy(e)) {
			m_pOGC->m_edge_copies.remove(original(e), this);
		}
	}
	for (node n : nodes) {
		if (!isDummy(n)) {
			m_pOGC->m_node_copies.remove(original(n), this);
		}
	}
}

}
