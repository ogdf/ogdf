/** \file
 * \brief TODO Document
 *
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
#pragma once

#include <ogdf/basic/CombinatorialEmbedding.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/cluster/sync_plan/PMatching.h>
#include <ogdf/cluster/sync_plan/SyncPlan.h>

#include <cstdint>
#include <functional>
#include <vector>

namespace ogdf::sync_plan {

void randomSyncPlanInstance(SyncPlan& pq, int pipe_count, int min_deg = 3) {
	for (int i = 0; i < pipe_count; ++i) {
		node u = pq.G->chooseNode(
				[&](node n) { return !pq.matchings.isMatchedPVertex(n) && n->degree() >= min_deg; });
		if (u == nullptr) {
			return;
		}
		node v = pq.G->chooseNode([&](node n) {
			return !pq.matchings.isMatchedPVertex(n) && n->degree() == u->degree() && n != u;
		});
		if (v == nullptr) {
			return;
		}
		pq.matchings.matchNodes(u, v);
	}
}

void addEdges(Graph* g, std::vector<edge>& added, int cnt) {
	CombinatorialEmbedding E(*g);
	for (int i = 0; i < cnt; ++i) {
		face f = E.chooseFace([](face f) { return f->size() > 3; });
		if (f == nullptr) {
			return;
		}
		adjEntry a = f->firstAdj();
		for (int i = 0; i < randomNumber(0, f->size() - 1); ++i) {
			a = a->faceCycleSucc();
		}
		adjEntry b = a;
		for (int i = 0; i < randomNumber(2, f->size() - 2); ++i) {
			b = b->faceCycleSucc();
		}
		OGDF_ASSERT(b != a);
		OGDF_ASSERT(b != a->faceCycleSucc());
		OGDF_ASSERT(b != a->faceCyclePred());
		if (a->theNode() == b->theNode()) {
			continue;
		}
		if (g->searchEdge(a->theNode(), b->theNode())) {
			continue;
		}
		added.push_back(E.splitFace(a, b));
	}
}

void randomSEFEInstanceBySharedGraph(Graph* sefe, EdgeArray<uint8_t>& edge_types, int edges1,
		int edges2) {
	OGDF_ASSERT(!sefe->empty());
	OGDF_ASSERT(isConnected(*sefe));
	OGDF_ASSERT(sefe->representsCombEmbedding());
	OGDF_ASSERT(edge_types.graphOf() == sefe);
	for (edge e : sefe->edges) {
		edge_types[e] = 3;
	}

	std::vector<edge> added1;
	addEdges(sefe, added1, edges1);
	Graph::HiddenEdgeSet h1(*sefe);
	for (edge e : added1) {
		edge_types[e] = 1;
		h1.hide(e);
	}

	std::vector<edge> added2;
	addEdges(sefe, added2, edges2);
	for (edge e : added2) {
		edge_types[e] = 2;
	}

	h1.restore();
}

void randomSEFEInstanceByUnionGraph(const Graph* sefe, EdgeArray<uint8_t>& edge_types,
		double frac_shared = 0.34, double frac_g1 = 0.33) {
	OGDF_ASSERT(edge_types.graphOf() == sefe);
	for (edge e : sefe->edges) {
		double r = randomDouble(0, 1);
		if (r < frac_shared) {
			edge_types[e] = 3;
		} else if (r < frac_shared + frac_g1) {
			edge_types[e] = 1;
		} else {
			edge_types[e] = 2;
		}
	}
}

}
