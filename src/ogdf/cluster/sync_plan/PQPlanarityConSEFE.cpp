#include <ogdf/basic/Graph_d.h>

#include <stdexcept>

#include "PQPlanarity.h"
#include "utils/GraphUtils.h"

PQPlanarity::PQPlanarity(const Graph* sefe, Graph* work, EdgeArray<uint8_t>& edge_types)
	: G(work)
	, GA(nullptr)
	, matchings(G)
	, partitions(G)
	, components(G)
	, is_wheel(*G, false)
#ifdef OGDF_DEBUG
	, consistency(*this)
#endif
{
	OGDF_ASSERT(work->empty());
	OGDF_ASSERT(edge_types.graphOf() == sefe);
	// TODO check that the shared graph is connected?
	// auto *op = new UndoInitConSEFE(); // TODO implement undo op to mirror embedding back to sefe Graph

	NodeArray<node> G1excl(*sefe, nullptr);
	NodeArray<node> G2excl(*sefe, nullptr);
	NodeArray<node> G1shared(*sefe, nullptr);
	NodeArray<node> G2shared(*sefe, nullptr);
	for (node n : sefe->nodes) {
		G1excl[n] = G->newNode();
		G2excl[n] = G->newNode();
		G1shared[n] = G->newNode();
		G2shared[n] = G->newNode();
	}

	for (edge e : sefe->edges) {
		OGDF_ASSERT(!e->isSelfLoop());
		switch (edge_types[e]) {
		case 1:
			G->newEdge(G1excl[e->source()], G1excl[e->target()]);
			G->newEdge(G1shared[e->source()], G->newNode());
			G->newEdge(G1shared[e->target()], G->newNode());
			break;
		case 2:
			G->newEdge(G2excl[e->source()], G2excl[e->target()]);
			G->newEdge(G2shared[e->source()], G->newNode());
			G->newEdge(G2shared[e->target()], G->newNode());
			break;
		case 3:
			G->newEdge(G1excl[e->source()], G1excl[e->target()]);
			G->newEdge(G2excl[e->source()], G2excl[e->target()]);
			G->newEdge(G1shared[e->source()], G2shared[e->source()]);
			G->newEdge(G1shared[e->target()], G2shared[e->target()]);
			break;
		default:
			throw runtime_error("illegal edge_type");
		}
	}

	for (node n : sefe->nodes) {
		G->reverseAdjEdges(G1shared[n]);
		G->reverseAdjEdges(G2shared[n]);
		matchings.matchNodes(G1excl[n], G1shared[n]);
		matchings.matchNodes(G2excl[n], G2shared[n]);
	}

	initComponents();
	matchings.rebuildHeap();
	// pushUndoOperationAndCheck(op);
}