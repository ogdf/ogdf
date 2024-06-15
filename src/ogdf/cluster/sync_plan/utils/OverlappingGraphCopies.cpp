#include <ogdf/cluster/sync_plan/basic/OverlappingGraphCopies.h>

const Graph& OverlappingGraphCopy::original() const { return *m_pOGC->G; }

node OverlappingGraphCopy::copy(node v) const {
	if (m_pOGC->node_copies.contains(v, this)) {
		node cv = m_pOGC->node_copies(v, this);
		OGDF_ASSERT(cv == nullptr || cv->graphOf() == this);
		return cv;
	} else {
		return nullptr;
	}
}

edge OverlappingGraphCopy::copy(edge e) const {
	if (m_pOGC->edge_copies.contains(e, this)) {
		edge ce = m_pOGC->edge_copies(e, this);
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
	OGDF_ASSERT(m_pOGC->node_copies(vOrig, this) == nullptr);
	m_pOGC->node_copies(vOrig, this) = v;
	return v;
}

edge OverlappingGraphCopy::newEdge(edge eOrig) {
	OGDF_ASSERT(eOrig != nullptr);
	OGDF_ASSERT(eOrig->graphOf() == &original());
	node s = m_pOGC->node_copies(eOrig->source(), this);
	node t = m_pOGC->node_copies(eOrig->target(), this);
	edge e = Graph::newEdge(s, t);
	m_eOrig[e] = eOrig;
	OGDF_ASSERT(m_pOGC->edge_copies(eOrig, this) == nullptr);
	m_pOGC->edge_copies(eOrig, this) = e;
	return e;
}

void OverlappingGraphCopy::delNode(node v) {
	if (!isDummy(v)) {
		m_pOGC->node_copies.remove(original(v), this);
	}
	Graph::delNode(v);
}

void OverlappingGraphCopy::delEdge(edge e) {
	if (!isDummy(e)) {
		m_pOGC->edge_copies.remove(original(e), this);
	}
	Graph::delEdge(e);
}

void OverlappingGraphCopy::unmap() {
	if (m_pOGC == nullptr) {
		return;
	}
	for (edge e : edges) {
		if (!isDummy(e)) {
			m_pOGC->edge_copies.remove(original(e), this);
		}
	}
	for (node n : nodes) {
		if (!isDummy(n)) {
			m_pOGC->node_copies.remove(original(n), this);
		}
	}
}
