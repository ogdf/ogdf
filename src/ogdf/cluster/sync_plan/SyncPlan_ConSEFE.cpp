/** \file
 * \brief Implementation of the Connected 2-SEFE-related functionality of SyncPlan.
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
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/cluster/sync_plan/PMatching.h>
#include <ogdf/cluster/sync_plan/SyncPlan.h>
#include <ogdf/cluster/sync_plan/SyncPlanComponents.h>

#include <cstdint>
#include <functional>
#include <ostream>
#include <stdexcept>
#include <vector>

using namespace ogdf::sync_plan::internal;

namespace ogdf::sync_plan {

struct UndoInitConSEFE : public SyncPlan::UndoOperation {
	Graph& sefe;
	Graph& work;
	EdgeArray<uint8_t>& edge_types;

	NodeArray<node> G1shared;
	NodeArray<node> G2shared;
	EdgeArray<adjEntry> Gedge;

	UndoInitConSEFE(Graph& _sefe, Graph& _work, EdgeArray<uint8_t>& _edge_types)
		: sefe(_sefe)
		, work(_work)
		, edge_types(_edge_types)
		, G1shared(sefe)
		, G2shared(sefe)
		, Gedge(work) { }

	void undo(SyncPlan& pq) override {
		std::vector<adjEntry> rot;
		for (node n : sefe.nodes) {
			rot.clear();
			rot.reserve(n->degree());
			adjEntry adj1_start = nullptr;
			for (adjEntry a : G1shared[n]->adjEntries) {
				if (edge_types[Gedge[a]] == 3) {
					adj1_start = a;
					break;
				}
			}
			if (adj1_start != nullptr) {
				adjEntry adj1 = adj1_start;
				adjEntry adj2 = adj1->twin();
				do {
					OGDF_ASSERT(adj1->theNode() == G1shared[n]);
					OGDF_ASSERT(adj2->theNode() == G2shared[n]);
					OGDF_ASSERT(Gedge[adj1] == Gedge[adj2]);
					rot.push_back(Gedge[adj1]);

					adj1 = adj1->cyclicSucc();
					while (edge_types[Gedge[adj1]] != 3 && adj1 != adj1_start) {
						OGDF_ASSERT(adj1->theNode() == G1shared[n]);
						OGDF_ASSERT(edge_types[Gedge[adj1]] == 1);
						rot.push_back(Gedge[adj1]);
						adj1 = adj1->cyclicSucc();
					}

					adj2 = adj2->cyclicPred();
					while (edge_types[Gedge[adj2]] != 3) {
						OGDF_ASSERT(adj2 != adj1_start->twin());
						OGDF_ASSERT(adj2->theNode() == G2shared[n]);
						OGDF_ASSERT(edge_types[Gedge[adj2]] == 2);
						rot.push_back(Gedge[adj2]);
						adj2 = adj2->cyclicPred();
					}
				} while (adj1 != adj1_start);
			} else {
				for (adjEntry adj1 : G1shared[n]->adjEntries) {
					OGDF_ASSERT(adj1->theNode() == G1shared[n]);
					OGDF_ASSERT(edge_types[Gedge[adj1]] == 1);
					rot.push_back(Gedge[adj1]);
				}
				for (adjEntry adj2 : G2shared[n]->adjEntries) {
					OGDF_ASSERT(adj2->theNode() == G2shared[n]);
					OGDF_ASSERT(edge_types[Gedge[adj2]] == 2);
					rot.push_back(Gedge[adj2]);
				}
			}
			sefe.sort(n, rot.begin(), rot.end());
		}
	}

	std::ostream& print(std::ostream& os) const override { return os << "UndoInitConSEFE"; }
};

SyncPlan::SyncPlan(Graph* sefe, Graph* work, EdgeArray<uint8_t>& edge_types)
	: G(work)
	, matchings(G)
	, partitions(G)
	, components(G)
	, deletedEdges(*G)
#ifdef OGDF_DEBUG
	, deletedNodes(*G)
#endif
	, GA(nullptr)
	, is_wheel(*G, false)
#ifdef OGDF_DEBUG
	, consistency(*this)
#endif
{
	OGDF_ASSERT(sefe != nullptr);
	OGDF_ASSERT(work != nullptr);
	OGDF_ASSERT(work->empty());
	OGDF_ASSERT(edge_types.graphOf() == sefe);
#ifdef OGDF_DEBUG
	{
		GraphCopySimple GC;
		GC.setOriginalGraph(*sefe);
		NodeArray<node> nodeMap(sefe);
		EdgeArray<edge> edgeMap(sefe);
		GC.insert(
				sefe->nodes.begin(), sefe->nodes.end(),
				[&edge_types](edge e) { return edge_types[e] == 3; }, nodeMap, edgeMap);
		safeForEach(GC.nodes, [&GC](node n) {
			if (n->degree() == 0) {
				GC.delNode(n);
			}
		});
		OGDF_ASSERT(isConnected(GC));
	}
#endif
	auto* op = new UndoInitConSEFE(*sefe, *work, edge_types);

	NodeArray<node> G1excl(*sefe, nullptr);
	NodeArray<node> G2excl(*sefe, nullptr);
	NodeArray<node>& G1shared = op->G1shared;
	NodeArray<node>& G2shared = op->G2shared;
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
			op->Gedge[G->newEdge(G1shared[e->source()], G->newNode())] = e->adjSource();
			op->Gedge[G->newEdge(G1shared[e->target()], G->newNode())] = e->adjTarget();
			break;
		case 2:
			G->newEdge(G2excl[e->source()], G2excl[e->target()]);
			op->Gedge[G->newEdge(G2shared[e->source()], G->newNode())] = e->adjSource();
			op->Gedge[G->newEdge(G2shared[e->target()], G->newNode())] = e->adjTarget();
			break;
		case 3:
			G->newEdge(G1excl[e->source()], G1excl[e->target()]);
			G->newEdge(G2excl[e->source()], G2excl[e->target()]);
			op->Gedge[G->newEdge(G1shared[e->source()], G2shared[e->source()])] = e->adjSource();
			op->Gedge[G->newEdge(G1shared[e->target()], G2shared[e->target()])] = e->adjTarget();
			break;
		default:
			throw std::runtime_error("illegal edge_type");
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
	pushUndoOperationAndCheck(op);
}

}
