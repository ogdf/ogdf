#include "utils/GraphUtils.h"

#include "utils/GraphIterators.h"
#include "utils/Logging.h"

FilteringBFSIterator begin(FilteringBFS& bfs) { return FilteringBFSIterator(&bfs); }

FilteringBFSIterator end(FilteringBFS& bfs) { return FilteringBFSIterator(nullptr); }

void moveEnd(Graph& G, edge e, node keep_end, node new_end) {
	if (e->source() == keep_end) {
		G.moveTarget(e, new_end);
	} else {
		OGDF_ASSERT(e->target() == keep_end);
		G.moveSource(e, new_end);
	}
}

void moveEnd(Graph& G, adjEntry keep_adj, adjEntry new_adj, Direction dir) {
	if (keep_adj->isSource()) {
		G.moveTarget(keep_adj->theEdge(), new_adj, dir);
	} else {
		G.moveSource(keep_adj->theEdge(), new_adj, dir);
	}
}

edge splitEdge(Graph& G, edge old_edge, node new_adj_to_source, node new_adj_to_target,
		int new_edge_idx) {
	auto old_source = old_edge->source();
	auto old_target = old_edge->target();

	edge new_edge = new_edge_idx < 0 ? G.newEdge(new_adj_to_target, old_target)
									 : G.newEdge(new_adj_to_target, old_target, new_edge_idx);
	G.moveAdjAfter(new_edge->adjTarget(), old_edge->adjTarget());
	G.moveTarget(old_edge, new_adj_to_source);

	OGDF_ASSERT(old_edge->source() == old_source);
	OGDF_ASSERT(old_edge->target() == new_adj_to_source);
	OGDF_ASSERT(new_edge->source() == new_adj_to_target);
	OGDF_ASSERT(new_edge->target() == old_target);
	return new_edge;
}

adjEntry splitEdge(Graph& G, adjEntry adj, node new_adj_to_node, node new_adj_to_twin,
		int new_edge_idx) {
	bool reverse = !adj->isSource();
	edge e = adj->theEdge();
	node n = adj->theNode();
	node t = adj->twinNode();
	if (reverse) {
		G.reverseEdge(e);
	}
	// A ----------e---------> D
	// A ---e--> B   C ---c--> D
	edge c = splitEdge(G, e, new_adj_to_node, new_adj_to_twin, new_edge_idx);
	OGDF_ASSERT(e->source() == n);
	OGDF_ASSERT(e->target() == new_adj_to_node);
	OGDF_ASSERT(c->source() == new_adj_to_twin);
	OGDF_ASSERT(c->target() == t);
	adjEntry adjt = c->adjTarget();
	if (reverse) {
		G.reverseEdge(e);
		G.reverseEdge(c);
	}
	return adjt;
}

bool joinEdge(Graph& G, edge u_e, edge v_e, node u, node v) {
	return joinEdge(G, u_e->getAdj(u), v_e->getAdj(v), u, v);
}

bool joinEdge(Graph& G, adjEntry u_adj, adjEntry v_adj, node u, node v) {
	OGDF_ASSERT(u_adj->theNode() == u);
	OGDF_ASSERT(v_adj->theNode() == v);
	bool opposing = (u_adj->isSource() == v_adj->isSource());
	moveEnd(G, u_adj->twin(), v_adj->twin());
	G.delEdge(v_adj->theEdge());
	return opposing;
}

void assertStarCentreAndRay(node centre, node ray) {
	bool ray_found = false;
	for (auto adj : centre->adjEntries) {
		OGDF_ASSERT(centre->degree() > adj->twinNode()->degree());
		if (adj->twinNode() == ray) {
			ray_found = true;
		}
		for (auto twin_adj : adj->twinNode()->adjEntries) {
			OGDF_ASSERT(twin_adj->twinNode() == centre);
		}
	}
	OGDF_ASSERT(!ray || ray_found);
}

node getCentreOfStar(node g_n) {
	OGDF_ASSERT(g_n->degree() >= 1);
	node g_adj = g_n->adjEntries.head()->twinNode();
	if (g_n->degree() > g_adj->degree()) {
		assertStarCentreAndRay(g_n, g_adj);
		return g_n;
	} else {
		OGDF_ASSERT(g_n->degree() < g_adj->degree()); // a parallel path is deemed biconnected, not a star
		assertStarCentreAndRay(g_adj, g_n);
		return g_adj;
	}
}

OrderComp compareCyclicOrder(node n, List<adjEntry>& o, bool full_check) {
	OGDF_ASSERT(n->degree() == o.size());
	adjEntry n_it = n->firstAdj();
	ListIterator<adjEntry> o_it = o.search(n_it);
	OGDF_ASSERT(o_it.valid());
	OGDF_ASSERT(*o_it == n_it);

	ListIterator<adjEntry> o_succ = o.cyclicSucc(o_it);
	bool reverse;
	if (*o_succ == n_it->cyclicSucc()) {
		reverse = false;
	} else if (*o_succ == n_it->cyclicPred()) {
		reverse = true;
	} else {
		return DIFFERENT;
	}

#ifndef OGDF_DEBUG
	if (full_check)
#endif
	{
		for (adjEntry n_it : n->adjEntries) {
			if (*o_it != n_it) {
#ifdef OGDF_DEBUG
				std::cout << "node: " << printIncidentEdges(n->adjEntries) << std::endl;
				std::cout << "list: " << printIncidentEdges(o) << std::endl;
				if (!full_check) {
					std::cout << "full_check == false, so broken/differing order wouldn't have been found in release mode"
							  << std::endl;
					OGDF_ASSERT(full_check);
				}
#endif
				return DIFFERENT;
			}
			if (reverse) {
				o_it = o.cyclicPred(o_it);
			} else {
				o_it = o.cyclicSucc(o_it);
			}
		}
	}

	return reverse ? REVERSED : SAME;
}

std::pair<node, node> split(Graph& G, PipeBij& bij, const EdgeArray<int>* split_idcs,
		const EdgeArray<bool>* split_reverse, int src_idx, int tgt_idx) {
	node src = src_idx < 0 ? G.newNode() : G.newNode(src_idx);
	node tgt = tgt_idx < 0 ? G.newNode() : G.newNode(tgt_idx);
	for (auto& pair : bij) {
		OGDF_ASSERT(pair.second == nullptr);
		int split_idx = split_idcs == nullptr ? -1 : (*split_idcs)[pair.first];
		pair.second = splitEdge(G, pair.first, src, tgt, split_idx);
		if (split_reverse != nullptr && (*split_reverse)[pair.first]) {
			G.reverseEdge(pair.second->theEdge());
		}
	}
	G.reverseAdjEdges(tgt);
	return std::pair<node, node>(src, tgt);
}

void join(Graph& G, node u, node v, PipeBij& bij, List<bool>* reverse_v) {
	OGDF_ASSERT(u->degree() == bij.size());
	OGDF_ASSERT(v->degree() == bij.size());
	for (auto& pair : bij) {
		adjEntry f = pair.first->twin();
		adjEntry s = pair.second->twin();
		bool rev = joinEdge(G, pair.first, pair.second, u, v);
		if (reverse_v) {
			reverse_v->pushBack(rev);
		}
		pair.first = f;
		pair.second = nullptr;
	}
	OGDF_ASSERT(u->degree() == 0);
	OGDF_ASSERT(v->degree() == 0);
	G.delNode(u);
	G.delNode(v);
}

void moveAdjToFront(Graph& G, adjEntry f) {
	List<adjEntry> adjs;
	f->theNode()->allAdjEntries(adjs);
	while (adjs.front() != f) {
		adjs.pushBack(adjs.popFrontRet());
	}
	G.sort(f->theNode(), adjs);
}

void moveAdjToBack(Graph& G, adjEntry b) {
	List<adjEntry> adjs;
	b->theNode()->allAdjEntries(adjs);
	while (adjs.back() != b) {
		adjs.pushBack(adjs.popFrontRet());
	}
	G.sort(b->theNode(), adjs);
}
