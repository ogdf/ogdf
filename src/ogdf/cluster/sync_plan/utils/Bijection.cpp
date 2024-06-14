#include "utils/Bijection.h"

#include <ogdf/basic/AdjEntryArray.h>

PipeBijIterator getPipeBijection(node u, node v) {
	OGDF_ASSERT(u->degree() == v->degree());
	return PipeBijIterator(u->adjEntries.begin(), u->adjEntries.end(), v->adjEntries.rbegin(),
			v->adjEntries.rend());
}

void getPipeBijection(node u, node v, PipeBij& out) {
	OGDF_ASSERT(u->graphOf() == v->graphOf());
	OGDF_ASSERT(u->degree() == v->degree());
	auto v_adj_it = v->adjEntries.rbegin();
	for (adjEntry u_adj : u->adjEntries) {
		OGDF_ASSERT(v_adj_it != v->adjEntries.rend());
		out.emplaceBack(u_adj, (*v_adj_it));
		v_adj_it++;
	}
	OGDF_ASSERT(v_adj_it == v->adjEntries.rend());
}

void getPipeBijection(node u, node v, AdjEntryArray<adjEntry>& out) {
	OGDF_ASSERT(u->graphOf() == v->graphOf());
	OGDF_ASSERT(u->graphOf() == out.graphOf());
	OGDF_ASSERT(u->degree() == v->degree());
	auto v_adj_it = v->adjEntries.rbegin();
	for (adjEntry u_adj : u->adjEntries) {
		OGDF_ASSERT(v_adj_it != v->adjEntries.rend());
		out[u_adj] = *v_adj_it;
		out[*v_adj_it] = u_adj;
		v_adj_it++;
	}
	OGDF_ASSERT(v_adj_it == v->adjEntries.rend());
}

void getPipeBijection(node u, node v, EdgeArray<edge>& out) {
	OGDF_ASSERT(u->graphOf() == v->graphOf());
	OGDF_ASSERT(u->graphOf() == out.graphOf());
	OGDF_ASSERT(u->degree() == v->degree());
	auto v_adj_it = v->adjEntries.rbegin();
	for (adjEntry u_adj : u->adjEntries) {
		OGDF_ASSERT(v_adj_it != v->adjEntries.rend());
		out[u_adj->theEdge()] = (*v_adj_it)->theEdge();
		v_adj_it++;
	}
	OGDF_ASSERT(v_adj_it == v->adjEntries.rend());
}

void getFrozenPipeBijection(node u, node v, FrozenPipeBij& bij) {
	OGDF_ASSERT(u->degree() == v->degree());
	auto v_adj_it = v->adjEntries.rbegin();
	for (adjEntry u_adj : u->adjEntries) {
		OGDF_ASSERT(v_adj_it != v->adjEntries.rend());
		bij.emplaceBack(u_adj->theEdge()->index(), (*v_adj_it)->theEdge()->index());
		++v_adj_it;
	}
	OGDF_ASSERT(v_adj_it == v->adjEntries.rend());
}

void freezePipeBijection(const PipeBij& in, FrozenPipeBij& out) {
	for (const PipeBijPair& pair : in) {
		out.emplaceBack(pair.first->theEdge()->index(), pair.second->theEdge()->index());
	}
}
