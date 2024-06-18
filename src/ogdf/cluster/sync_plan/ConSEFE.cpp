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
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/basic.h>
#include <ogdf/cluster/sync_plan/PMatching.h>
#include <ogdf/cluster/sync_plan/SyncPlan.h>
#include <ogdf/cluster/sync_plan/SyncPlanComponents.h>

#include <cstdint>
#include <stdexcept>

using namespace ogdf::sync_plan::internal;

namespace ogdf::sync_plan {

SyncPlan::SyncPlan(const Graph* sefe, Graph* work, EdgeArray<uint8_t>& edge_types)
	: G(work)
	, matchings(G)
	, partitions(G)
	, components(G)
	, GA(nullptr)
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
	// pushUndoOperationAndCheck(op);
}

}
