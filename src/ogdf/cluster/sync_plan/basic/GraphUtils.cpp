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
#include <ogdf/basic/List.h>
#include <ogdf/basic/basic.h>
#include <ogdf/cluster/sync_plan/basic/GraphUtils.h>

#include <functional>

namespace ogdf {

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
		edge new_edge) {
	auto old_source = old_edge->source();
	auto old_target = old_edge->target();

	if (new_edge == nullptr) {
		new_edge = G.newEdge(new_adj_to_target, old_target);
	} else {
		G.moveSource(new_edge, new_adj_to_target);
		G.moveTarget(new_edge, old_target);
	}
	G.moveAdjAfter(new_edge->adjTarget(), old_edge->adjTarget());
	G.moveTarget(old_edge, new_adj_to_source);

	OGDF_ASSERT(old_edge->source() == old_source);
	OGDF_ASSERT(old_edge->target() == new_adj_to_source);
	OGDF_ASSERT(new_edge->source() == new_adj_to_target);
	OGDF_ASSERT(new_edge->target() == old_target);
	return new_edge;
}

adjEntry splitEdge(Graph& G, adjEntry adj, node new_adj_to_node, node new_adj_to_twin, edge new_edge) {
	bool reverse = !adj->isSource();
	edge e = adj->theEdge();
	node n = adj->theNode();
	node t = adj->twinNode();
	if (reverse) {
		G.reverseEdge(e);
	}
	// A ----------e---------> D
	// A ---e--> B   C ---c--> D
	edge c = splitEdge(G, e, new_adj_to_node, new_adj_to_twin, new_edge);
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
	return joinEdge(G, u_adj, v_adj, u, v, [&G](edge e) { G.delEdge(e); });
}

bool joinEdge(Graph& G, edge u_e, edge v_e, node u, node v,
		const std::function<void(edge)>& deleteEdge) {
	return joinEdge(G, u_e->getAdj(u), v_e->getAdj(v), u, v, deleteEdge);
}

bool joinEdge(Graph& G, adjEntry u_adj, adjEntry v_adj, node u, node v,
		const std::function<void(edge)>& deleteEdge) {
	OGDF_ASSERT(u_adj->theNode() == u);
	OGDF_ASSERT(v_adj->theNode() == v);
	bool opposing = (u_adj->isSource() == v_adj->isSource());
	moveEnd(G, u_adj->twin(), v_adj->twin());
	deleteEdge(v_adj->theEdge());
	return opposing;
}

void assertStarCentreAndRay(node centre, node ray) {
#ifdef OGDF_DEBUG
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
#endif
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
		return OrderComp::DIFFERENT;
	}

#ifndef OGDF_DEBUG
	if (full_check)
#endif
	{
		for (adjEntry a : n->adjEntries) {
			if (*o_it != a) {
				// #ifdef OGDF_DEBUG
				// std::cout << "node: " << ogdf::sync_plan::printIncidentEdges(n->adjEntries) << std::endl;
				// std::cout << "list: " << ogdf::sync_plan::printIncidentEdges(o) << std::endl;
				// if (!full_check) {
				// 	std::cerr << "Order differs in the middle, but full_check == false, so broken/differing order wouldn't have been found in release mode!"
				// 			  << std::endl;
				// 	OGDF_ASSERT(false);
				// }
				// #endif
				OGDF_ASSERT(full_check);
				return OrderComp::DIFFERENT;
			}
			if (reverse) {
				o_it = o.cyclicPred(o_it);
			} else {
				o_it = o.cyclicSucc(o_it);
			}
		}
	}

	return reverse ? OrderComp::REVERSED : OrderComp::SAME;
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

}
