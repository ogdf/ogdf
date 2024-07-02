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
#include <ogdf/basic/GraphSets.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/pctree/util/FilteringBFS.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/cluster/sync_plan/SyncPlanComponents.h>
#include <ogdf/cluster/sync_plan/SyncPlan_operation/Encapsulate.h>
#include <ogdf/decomposition/BCTree.h>

#include <functional>
#include <ostream>
#include <utility>

using namespace ogdf::sync_plan::internal;

namespace ogdf::sync_plan {

std::function<std::ostream&(std::ostream&)> SyncPlanComponents::fmtBCNode(node bc) const {
	OGDF_ASSERT(bc == nullptr || bc->graphOf() == &bcTree());
	return [bc, this](std::ostream& ss) -> std::ostream& {
		if (bc != nullptr) {
			ss << "{" << bcNodeType(bc) << " #" << bc->index() << " °" << bc->degree() << " @"
			   << bcConnectedId(bc);
			node repr = bcRepr(bc);
			if (repr != nullptr) {
				ss << " r" << repr;
			}
			ss << " s" << bcSize(bc);
			ss << "}";
		} else {
			ss << "{NULL}";
		}
		return ss;
	};
}

BCTree::GNodeType SyncPlanComponents::nodeType(node n) const {
	switch (bc_type[biconnectedComponent(n)]) {
	default:
		OGDF_ASSERT(false);
	case BCTree::BNodeType::BComp:
		return BCTree::GNodeType::Normal;
	case BCTree::BNodeType::CComp:
		return BCTree::GNodeType::CutVertex;
	}
}

node SyncPlanComponents::biconnectedComponent(node n) const {
	OGDF_ASSERT(n->graphOf() != &BC);
	OGDF_ASSERT(n->graphOf() == G);
	node bc = g_bc[n];
	OGDF_ASSERT(bc->graphOf() == &BC);
	return bc;
}

BCTree::BNodeType SyncPlanComponents::bcNodeType(node n) const {
	OGDF_ASSERT(n->graphOf() != G);
	OGDF_ASSERT(n->graphOf() == &BC);
	return bc_type[n];
}

int SyncPlanComponents::bcSize(node n) const {
	OGDF_ASSERT(n->graphOf() != G);
	OGDF_ASSERT(n->graphOf() == &BC);
	return bc_size[n];
}

int SyncPlanComponents::bcConnectedId(node n) const {
	OGDF_ASSERT(n->graphOf() != G);
	OGDF_ASSERT(n->graphOf() == &BC);
	return bc_conn_id[n];
}

node SyncPlanComponents::bcRepr(node bc) const {
	OGDF_ASSERT(bc->graphOf() != G);
	OGDF_ASSERT(bc->graphOf() == &BC);
	node n = bc_g[bc];
	if (isCutComponent(bc) || bc->degree() == 0) {
		OGDF_ASSERT(n->graphOf() == G);
		OGDF_ASSERT(g_bc[n] == bc);
	} else {
		OGDF_ASSERT(n == nullptr);
	}
	return n;
}

node SyncPlanComponents::findOtherRepr(node bc) const {
	for (adjEntry adj : bcRepr(bc)->adjEntries) {
		if (biconnectedComponent(adj->twinNode()) == bc) {
			return adj->twinNode();
		}
	}
	return nullptr;
}

std::pair<edge, edge> SyncPlanComponents::graphEdgeToBCEdge(node bc_src, node bc_tgt) const {
	OGDF_ASSERT(bc_src->graphOf() == &BC);
	OGDF_ASSERT(bc_tgt->graphOf() == &BC);
	if (bc_src->degree() > bc_tgt->degree()) {
		const std::pair<edge, edge>& res = graphEdgeToBCEdge(bc_tgt, bc_src);
		if (res.second == nullptr) {
			return res;
		} else {
			return std::make_pair(res.second, res.first);
		}
	}
	int c = counter++;
	for (adjEntry adj : bc_src->adjEntries) {
		if (adj->twinNode() == bc_tgt) {
			return std::make_pair(adj->theEdge(), nullptr);
		}
		int deg = adj->twinNode()->degree();
		if (deg == 2) {
			adjEntry next_adj = adj->twin()->cyclicSucc();
			if (next_adj->twinNode() == bc_tgt) {
				return std::make_pair(adj->theEdge(), next_adj->theEdge());
			}
		} else if (deg > 2) {
			marker[adj->twinNode()].first = c;
			marker[adj->twinNode()].second = adj->theEdge();
		}
	}
	for (adjEntry adj : bc_tgt->adjEntries) {
		if (marker[adj->twinNode()].first == c) {
			OGDF_ASSERT(marker[adj->twinNode()].second->isIncident(bc_src));
			return std::make_pair(marker[adj->twinNode()].second, adj->theEdge());
		}
	}
	OGDF_ASSERT(false);
	return std::make_pair(nullptr, nullptr);
}

node SyncPlanComponents::findCommonBiconComp(node bc_cut1, node bc_cut2) const {
	OGDF_ASSERT(isCutComponent(bc_cut1));
	OGDF_ASSERT(isCutComponent(bc_cut2));
	const std::pair<edge, edge>& bcEdge = graphEdgeToBCEdge(bc_cut1, bc_cut2);
	OGDF_ASSERT(bcEdge.second != nullptr);
	node common = bcEdge.first->commonNode(bcEdge.second);
	OGDF_ASSERT(common != nullptr);
	OGDF_ASSERT(!isCutComponent(common));
	return common;
}

void SyncPlanComponents::biconReprNodes(node bicon, SList<node>& nodes) const {
	node repr = bcRepr(bicon);
	if (repr != nullptr) {
		nodes.pushBack(repr);
	} else {
		for (adjEntry adj : bicon->adjEntries) {
			node r = bcRepr(adj->twinNode());
			OGDF_ASSERT(r != nullptr);
			nodes.pushBack(r);
		}
	}
	OGDF_ASSERT(!nodes.empty());
}

FilteringBFS SyncPlanComponents::nodesInBiconnectedComponent(node bicon) const {
	OGDF_ASSERT(bicon->graphOf() == &BC);
	SList<node> nodes;
	biconReprNodes(bicon, nodes);
	return FilteringBFS(*G, nodes,
			[this, bicon](adjEntry adj) -> bool { return g_bc[adj->twinNode()] == bicon; });
}

void SyncPlanComponents::makeRepr(node bc, node n) {
	OGDF_ASSERT(bc->graphOf() != G);
	OGDF_ASSERT(bc->graphOf() == &BC);
	bc_g[bc] = n;
	OGDF_ASSERT(bcRepr(bc) == n);
}

void SyncPlanComponents::nodeInserted(node g_n, node bc_n) {
	OGDF_ASSERT(g_n->graphOf() != &BC);
	OGDF_ASSERT(g_n->graphOf() == G);
	OGDF_ASSERT(bc_n->graphOf() != G);
	OGDF_ASSERT(bc_n->graphOf() == &BC);
	g_bc[g_n] = bc_n;
	bc_size[bc_n]++;
}

void SyncPlanComponents::insert(BCTree& tmp_bc) {
	OGDF_ASSERT(&tmp_bc.originalGraph() == G);
	NodeArray<node> bc_map(tmp_bc.bcTree(), nullptr);
	EdgeArray<edge> dummy(tmp_bc.bcTree(), nullptr);
	NodeArray<int> bc_conn(tmp_bc.bcTree(), -1);
	int new_conn = connectedComponents(tmp_bc.bcTree(), bc_conn);
	BC.insert(tmp_bc.bcTree(), bc_map, dummy);
	for (node h_node : tmp_bc.auxiliaryGraph().nodes) {
		node g_node = tmp_bc.original(h_node);
		node tmp_bc_node = tmp_bc.bcproper(g_node);
		node bc_node = bc_map[tmp_bc_node];
		g_bc[g_node] = bc_node;
		bc_g[bc_node] = g_node;
		OGDF_ASSERT(biconnectedComponent(g_node) == bc_node);
	}
	for (node tmp_bc_node : tmp_bc.bcTree().nodes) {
		BCTree::BNodeType type = tmp_bc.typeOfBNode(tmp_bc_node);
		node bc_node = bc_map[tmp_bc_node];
		bc_conn_id[bc_node] = conn_next_id + bc_conn[tmp_bc_node];
		bc_type[bc_node] = type;
		bc_size[bc_node] = tmp_bc.numberOfNodes(tmp_bc_node);
		if (type == BCTree::BNodeType::BComp && bc_node->degree() > 0) {
			bc_g[bc_node] = nullptr;
			OGDF_ASSERT(bcRepr(bc_node) == nullptr);
		} else {
			OGDF_ASSERT(biconnectedComponent(bcRepr(bc_node)) == bc_node);
		}
	}
	conn_count += new_conn;
	conn_next_id += new_conn;
}

void SyncPlanComponents::cutReplacedByWheel(node centre,
		const NodeArray<SList<adjEntry>>& block_neigh) {
	// we need to merge all blocks that have at least 2 edges connecting them to the center vertex
	node centre_bc = biconnectedComponent(centre);
	List<adjEntry> bc_adjs;
	centre_bc->allAdjEntries(bc_adjs);
	for (adjEntry bc_adj : bc_adjs) {
		node bicon_neigh = bc_adj->twinNode();
		OGDF_ASSERT(!isCutComponent(bicon_neigh));
		edge bc_edge = bc_adj->theEdge();

		if (block_neigh[bicon_neigh].size() == 1) {
			// a block with only one edge will lead to the respective rim vertex of the wheel again being a CV
			if (bc_adj->isSource()) {
				BC.reverseEdge(bc_edge);
			}
			edge new_edge = BC.split(bc_edge);
			OGDF_ASSERT(bicon_neigh == bc_edge->source());
			OGDF_ASSERT(bc_edge->target() == new_edge->source());
			OGDF_ASSERT(new_edge->target() == centre_bc);
			node new_cut = new_edge->source();

			bc_type[new_cut] = BCTree::BNodeType::CComp;
			node spike = block_neigh[bicon_neigh].front()->twinNode();
			bc_g[new_cut] = spike;
			g_bc[spike] = new_cut;
			bc_conn_id[new_cut] = bcConnectedId(centre_bc);
			bc_size[new_cut] = 1;
		} else {
			OGDF_ASSERT(block_neigh[bicon_neigh].size() >= 2);
			FilteringBFS bfs = nodesInBiconnectedComponent(bicon_neigh);
			// centre is now separated from the bicon nodes by the wheel nodes,
			// so we need to add those to init the bfs
			for (adjEntry adj : block_neigh[bicon_neigh]) {
				bfs.append(adj->twinNode());
			}
			for (node n : bfs) {
				if (biconnectedComponent(n) == bicon_neigh) {
					g_bc[n] = centre_bc;
				}
			}

			if (bc_edge->source() != centre_bc) {
				BC.reverseEdge(bc_edge);
			}
			OGDF_ASSERT(bc_edge->source() == centre_bc);
			bc_size[bc_edge->source()] += bc_size[bc_edge->target()] - 1;
			OGDF_IF_DBG(node contraced_node =) BC.contract(bc_edge);
			OGDF_ASSERT(contraced_node == centre_bc);
		}
	}
	bc_type[centre_bc] = BCTree::BNodeType::BComp;
	bc_g[centre_bc] = centre_bc->degree() > 0 ? nullptr : centre;
}

void SyncPlanComponents::relabelExplodedStar(node center1, node center2, List<node> remnants) {
	List<adjEntry> todelete;

	OGDF_ASSERT(center1->graphOf() == &BC);
	center1->allAdjEntries(todelete);
	for (adjEntry adj : todelete) {
		OGDF_ASSERT(adj->twinNode()->degree() == 1);
		BC.delNode(adj->twinNode());
	}
	OGDF_ASSERT(center1->degree() == 0);
	BC.delNode(center1);
	conn_count -= 1;

	if (center2 != nullptr) {
		OGDF_ASSERT(center2->graphOf() == &BC);
		center2->allAdjEntries(todelete);
		for (adjEntry adj : todelete) {
			OGDF_ASSERT(adj->twinNode()->degree() == 1);
			BC.delNode(adj->twinNode());
		}
		OGDF_ASSERT(center2->degree() == 0);
		BC.delNode(center2);
		conn_count -= 1;
	}

	for (node n : remnants) {
		OGDF_ASSERT(n->graphOf() == G);
	}
	BCTree bc(*G, remnants);
	insert(bc);
}

void SyncPlanComponents::preJoin(node keep, node merge) {
	// room for improvement: use union-find for relabeling
	int cc = bcConnectedId(keep);
	OGDF_ASSERT(cc != bcConnectedId(merge));
	for (node n : nodesInBiconnectedComponent(merge)) {
		if (biconnectedComponent(n) == merge) {
			g_bc[n] = keep;
		}
	}
	for (node bicon : FilteringBFS(BC, {merge})) {
		bc_conn_id[bicon] = cc;
	}

	bc_size[keep] += bc_size[merge] - 2;
	node n_bc = BC.contract(BC.newEdge(keep, merge));
	OGDF_ASSERT(n_bc == keep);
	conn_count--;
}

void SyncPlanComponents::postSplitOffEncapsulatedBlock(node cut, EncapsulatedBlock& block) {
	node bc_cut = biconnectedComponent(cut);
	int conn_id = bcConnectedId(bc_cut);
	OGDF_ASSERT(isCutComponent(bc_cut));
	OGDF_ASSERT(!isCutComponent(block.bicon));

	g_bc[block.bicon_rep] = block.bicon;
	BC.delEdge(BC.searchEdge(bc_cut, block.bicon));
	if (block.bicon->degree() == 0) {
		bc_g[block.bicon] = block.bicon_rep;
	}

	for (node bicon : FilteringBFS(BC, {(block.bicon)})) {
		bc_conn_id[bicon] = conn_next_id;
	}
	OGDF_ASSERT(bcConnectedId(block.bicon) == conn_next_id);
	OGDF_ASSERT(bcConnectedId(bc_cut) == conn_id);
	conn_next_id++;
	conn_count++;

	node star_bc = BC.newNode();
	g_bc[block.star_rep] = star_bc;
	bc_type[star_bc] = BCTree::BNodeType::BComp;
	bc_conn_id[star_bc] = bc_conn_id[bc_cut];
	bc_size[star_bc] = 2;
	BC.newEdge(bc_cut, star_bc);
}

void SyncPlanComponents::labelIsolatedNodes() {
	for (node n : G->nodes) {
		if (n->degree() == 0) {
			if (g_bc[n] == nullptr) {
				node new_bc = BC.newNode();
				g_bc[n] = new_bc;
				bc_type[new_bc] = BCTree::BNodeType::BComp;
				bc_g[new_bc] = n;
				bc_conn_id[new_bc] = conn_next_id;
				bc_size[new_bc] = 1;
				conn_count++;
				conn_next_id++;
			} else {
				OGDF_ASSERT(biconnectedComponent(n)->degree() == 0);
			}
		}
	}
}

BiconnectedIsolation::BiconnectedIsolation(SyncPlanComponents& comps, node bicon)
	: m_comps(comps)
	, m_bicon(bicon)
	, m_to_restore(comps.graph())
	, m_adjEntries(comps.graph())
	, m_hiddenEdges(comps.graph()) {
	OGDF_ASSERT(!comps.isCutComponent(bicon));
	if (m_bicon->degree() == 0) {
		return;
	}

	NodeArray<bool> include_bc(comps.bcTree(), false);
	include_bc[bicon] = true;
	for (adjEntry adj : bicon->adjEntries) {
		include_bc[adj->twinNode()] = true;
	}

	SList<node> pending;
	comps.biconReprNodes(bicon, pending);

	while (!pending.empty()) {
		node n = pending.popFrontRet();
		if (!m_adjEntries[n].empty()) {
			continue;
		}
		n->allAdjEntries(m_adjEntries[n]);

		for (adjEntry adj : m_adjEntries[n]) {
			node twin = adj->twinNode();
			if (!include_bc[comps.biconnectedComponent(twin)]) {
				if (m_adjEntries[twin].empty()) {
					twin->allAdjEntries(m_adjEntries[twin]);
				}
				m_hiddenEdges.hide(adj->theEdge());
				m_to_restore.insert(adj->theNode());
				m_to_restore.insert(adj->twinNode());
			} else if (m_adjEntries[twin].empty()) {
				pending.pushBack(twin);
			}
		}
	}

#ifdef OGDF_DEBUG
	EdgeArray<int> bc_id(comps.graph());
	biconnectedComponents(comps.graph(), bc_id);

	node repr = comps.bcRepr(bicon);
	if (repr == nullptr) {
		repr = comps.bcRepr(bicon->adjEntries.head()->twinNode());
	}
	int the_id = bc_id[repr->adjEntries.head()];
	for (node n : FilteringBFS(comps.graph(), {repr})) {
		for (adjEntry adj : n->adjEntries) {
			OGDF_ASSERT(bc_id[adj] == the_id);
		}
	}
	for (node n : m_to_restore.nodes()) {
		if (include_bc[comps.biconnectedComponent(n)]) {
			for (adjEntry adj : n->adjEntries) {
				OGDF_ASSERT(bc_id[adj] == the_id);
			}
		}
	}
#endif
}

void BiconnectedIsolation::restore(bool restore_embedding) {
	if (m_hiddenEdges.size() == 0) {
		return;
	}
	m_hiddenEdges.restore();
	if (restore_embedding) {
		for (node n : m_to_restore.nodes()) {
			m_comps.graph().sort(n, m_adjEntries[n]);
		}
	}
}

}
