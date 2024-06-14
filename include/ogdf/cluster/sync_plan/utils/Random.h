#pragma once

#include <ogdf/basic/CombinatorialEmbedding.h>
#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/Graph_d.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/simple_graph_alg.h>

#include "../PQPlanarity.h"
#include "GraphIterators.h"

using namespace ogdf;
using namespace std;

void randomPQPlanInstance(PQPlanarity& pq, int pipe_count, int min_deg = 3) {
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

void addEdges(Graph* g, vector<edge>& added, int cnt) {
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

	vector<edge> added1;
	addEdges(sefe, added1, edges1);
	Graph::HiddenEdgeSet h1(*sefe);
	for (edge e : added1) {
		edge_types[e] = 1;
		h1.hide(e);
	}

	vector<edge> added2;
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
