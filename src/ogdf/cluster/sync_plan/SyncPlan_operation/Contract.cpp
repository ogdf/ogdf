/** \file
 * \brief Implementation of the SyncPlan::contract operation and its UndoOperation.
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
#include <ogdf/basic/GraphSets.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/Logger.h>
#include <ogdf/basic/basic.h>
#include <ogdf/cluster/sync_plan/PMatching.h>
#include <ogdf/cluster/sync_plan/SyncPlan.h>
#include <ogdf/cluster/sync_plan/SyncPlanComponents.h>
#include <ogdf/cluster/sync_plan/utils/Bijection.h>
#include <ogdf/cluster/sync_plan/utils/Logging.h>

#include <functional>
#include <ostream>
#include <string>
#include <utility>

using namespace ogdf::sync_plan::internal;

namespace ogdf::sync_plan {
using internal::operator<<;

struct FindType {
	bool selected : 1;
	bool discovered : 1;
	bool visited_src : 1;
	bool visited_tgt : 1;

	FindType() : selected(false), discovered(false), visited_src(false), visited_tgt(false) { }

	friend std::ostream& operator<<(std::ostream& os, const FindType& type) {
		os << (type.selected ? "S" : "-") << (type.discovered ? "D" : "-")
		   << (type.visited_src ? "V" : "-") << (type.visited_tgt ? "T" : "-");
		return os;
	}

	inline void discover() {
		OGDF_ASSERT(selected);
		OGDF_ASSERT(!discovered);
		discovered = true;
	}

	inline bool discoverable() const { return selected && !discovered; }

	inline bool isVisited(adjEntry adj) const {
		if (adj->isSource()) {
			return visited_src;
		} else {
			return visited_tgt;
		}
	}

	inline void visit(adjEntry adj) {
		OGDF_ASSERT(!isVisited(adj));
		if (adj->isSource()) {
			visited_src = true;
		} else {
			visited_tgt = true;
		}
	}
};

adjEntry getSelectedAdj(const NodeArray<FindType>& node_types, edge e) {
	if (node_types[e->source()].selected) {
		OGDF_ASSERT(!node_types[e->target()].selected);
		return e->adjSource();
	} else {
		OGDF_ASSERT(node_types[e->target()].selected);
		return e->adjTarget();
	}
}

#ifdef OGDF_DEBUG
#	define cean_log                           \
		if (log.is_lout(Logger::Level::Minor)) \
		log.lout(Logger::Level::Minor)
#else
#	define cean_log \
		if (false)   \
		log.lout(Logger::Level::Minor)
#endif

//! Find a cut through all edges in an embedded bipartite graph by finding a spanning tree in the dual trough all SELECTED nodes.
//! Assumes that the SELECTED nodes form one half of the node bipartition and all edges are SELECTED.
void findBipartiteEdgeCut(Logger& log, NodeArray<FindType>& node_types,
		EdgeArray<FindType>& edge_types, adjEntry start_adj, PipeBij& out_edges) {
	node_types[start_adj->theNode()].discover();
	cean_log << "findBipartiteEdgeCut(node " << start_adj->theNode()->index()
			 << " starting with edge " << start_adj->theEdge()->index() << " leading to node "
			 << start_adj->twinNode()->index() << ", " << out_edges.size() << " out edges)"
			 << std::endl;
#ifdef OGDF_DEBUG
	Logger::Indent _(&log);
#endif
	adjEntry node_adj = start_adj;
	do {
		cean_log << "adding edge " << node_adj->theEdge()->index() << " (" << edge_types[node_adj]
				 << ") leading to node " << node_adj->twinNode()->index() << " ("
				 << node_types[node_adj->twinNode()] << ")"
				 << " to out list at position " << out_edges.size() << std::endl;
		edge_types[node_adj].discover();
		out_edges.emplaceFront(node_adj, nullptr);

		if (!edge_types[node_adj].isVisited(node_adj)) {
			edge_types[node_adj].visit(node_adj);
			cean_log << "scanning face" << std::endl;
#ifdef OGDF_DEBUG
			Logger::Indent __(&log);
#endif
			// following the face of start_adj->cyclicPred() clockwise would lead us back to start_adj
			for (adjEntry face_adj = node_adj->clockwiseFaceSucc(); face_adj != node_adj;
					face_adj = face_adj->clockwiseFaceSucc()) {
				cean_log << "edge " << face_adj->theEdge()->index() << " (" << edge_types[face_adj]
						 << ") from node " << face_adj->theNode()->index() << " ("
						 << node_types[face_adj->theNode()] << ") to node "
						 << face_adj->twinNode()->index() << " ("
						 << node_types[face_adj->twinNode()] << ")"
						 << (edge_types[face_adj].isVisited(face_adj) ? " already visited" : "")
						 << std::endl;
				if (edge_types[face_adj].isVisited(face_adj)) {
					break;
				} else {
					edge_types[face_adj].visit(face_adj);
				}
				if (node_types[face_adj->theNode()].discoverable()) {
					findBipartiteEdgeCut(log, node_types, edge_types, face_adj->cyclicSucc(),
							out_edges);
				}
			}
		}

		node_adj = node_adj->cyclicSucc();
	} while (node_adj != start_adj);
	cean_log << "done" << std::endl;
}

//! Find a cut through all SELECTED edges in an embedded biconnected graph by simply walking along faces.
//! Assumes that the two respective subgraphs separated by the SELECTED edges are connected.
void findBiconnectedEdgeCut(Logger& log, NodeArray<FindType>& node_types,
		EdgeArray<FindType>& edge_types, adjEntry start_adj, PipeBij& out_edges) {
	cean_log << "findBiconnectedEdgeCut(" << start_adj << ")" << std::endl;
#ifdef OGDF_DEBUG
	Logger::Indent _(&log);
	int i = 0;
#endif
	adjEntry node_adj = start_adj;
	while (true) {
		if (edge_types[node_adj].discovered) {
			cean_log << "adj " << node_adj << ": reached discovered edge, done" << std::endl;
			OGDF_ASSERT(node_adj->theEdge() == start_adj->theEdge());
			OGDF_ASSERT(node_adj == start_adj);
			return;
		} else if (edge_types[node_adj].selected) {
			cean_log << "adj " << node_adj << ": edge " << node_adj->theEdge() << " and node "
					 << node_adj->theNode() << " are selected, "
					 << "adding to out list at position " << out_edges.size()
					 << " and switching face" << std::endl;
			edge_types[node_adj].discover();
			out_edges.emplaceFront(getSelectedAdj(node_types, node_adj->theEdge()), nullptr);
			node_adj = node_adj->twin();
		} else {
			cean_log << "adj " << node_adj << std::endl;
		}
		node_adj = node_adj->faceCycleSucc();
		OGDF_IF_DBG(i++);
		OGDF_ASSERT(i < 2 * start_adj->graphOf()->numberOfEdges());
	}
}

class UndoContract : public SyncPlan::UndoOperation {
public:
	int u_idx, v_idx;
	FrozenPipeBij bij;
	List<int> u_neigh_idcs;
	List<bool> reverse_v_edges;
	bool biconnected;

	UndoContract(node uIdx, node vIdx, const PipeBij& tBij, const List<node>& uNeighs,
			bool _biconnected)
		: u_idx(uIdx->index()), v_idx(vIdx->index()), biconnected(_biconnected) {
		freezePipeBijection(tBij, bij);
		for (node u_neigh : uNeighs) {
			u_neigh_idcs.pushBack(u_neigh->index());
		}
	}

	void undo(SyncPlan& pq) override {
		// SYNCPLAN_PROFILE_START("undo-contract")
		pq.log.lout(Logger::Level::High)
				<< "UNDO "
				<< (biconnected ? "CONTRACT/JOIN BICONNECTED" : "ENCAPSULATE AND CONTRACT/JOIN STARS")
				<< " with " << bij.size() << " cut edges from " << u_neigh_idcs.size()
				<< " neighbours." << std::endl;
		Logger::Indent _(&pq.log);

		EdgeArray<FindType> edge_types(*pq.G);
		EdgeArray<edge> partner_edge_idx(*pq.G, nullptr);
		EdgeArray<bool> reverse_partner_edge(*pq.G, false);
		auto reverse_v_it = reverse_v_edges.begin();
		for (const FrozenPipeBijPair& pair : bij) {
			OGDF_ASSERT(reverse_v_it != reverse_v_edges.end());
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
			OGDF_ASSERT(!edge_types[pair.first].selected);
			edge_types[pair.first].selected = true;
			partner_edge_idx[pair.first] = pq.edgeFromIndex(pair.second);
			reverse_partner_edge[pair.first] = *reverse_v_it;
#pragma GCC diagnostic pop
			++reverse_v_it;
		}
		OGDF_ASSERT(reverse_v_it == reverse_v_edges.end());

		NodeArray<FindType> node_types(*pq.G);
		for (int u_neigh_idx : u_neigh_idcs) {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
			node_types[u_neigh_idx].selected = true;
#pragma GCC diagnostic pop
		}

		PipeBij cut_edges;
		if (bij.size() <= 2) {
			for (const FrozenPipeBijPair& pair : bij) {
				cut_edges.emplaceBack(getSelectedAdj(node_types, pq.edgeFromIndex(pair.first)),
						nullptr);
			}
		} else if (biconnected) {
			edge e = pq.edgeFromIndex(bij.front().first);
			OGDF_ASSERT(edge_types[e].selected);
			adjEntry start = getSelectedAdj(node_types, e);

			pq.log.lout(Logger::Level::Medium)
					<< "Collecting biconnected cut edges starting at edge #" << e->index() << " " << e
					<< " adjacent to " << pq.fmtPQNode(start->theNode(), false) << "." << std::endl;
			// SYNCPLAN_PROFILE_START("undo-contract-collect-bicon")
			findBiconnectedEdgeCut(pq.log, node_types, edge_types, start->twin(), cut_edges);
			// SYNCPLAN_PROFILE_STOP("undo-contract-collect-bicon")
			pq.log.lout(Logger::Level::Medium)
					<< "Collected " << cut_edges.size() << " cut edges." << std::endl;
		} else {
			for (const FrozenPipeBijPair& pair : bij) {
				edge e = pq.edgeFromIndex(pair.first);
				if (edge_types[e].discovered) {
					continue;
				}
				OGDF_ASSERT(edge_types[e].selected);
				adjEntry start = getSelectedAdj(node_types, e);

				pq.log.lout(Logger::Level::Medium)
						<< "Collecting bipartite cut edges starting at edge #" << e->index() << " "
						<< e << " adjacent to " << pq.fmtPQNode(start->theNode(), false) << "."
						<< std::endl;
				// SYNCPLAN_PROFILE_START("undo-contract-collect-stars")
				findBipartiteEdgeCut(pq.log, node_types, edge_types, start, cut_edges);
				// SYNCPLAN_PROFILE_STOP("undo-contract-collect-stars")
				pq.log.lout(Logger::Level::Medium)
						<< "Collected " << cut_edges.size() << " cut edges." << std::endl;
#ifndef OGDF_DEBUG
				if (cut_edges.size() == bij.size()) {
					break;
				}
#endif
			}
		}
		pq.log.lout(Logger::Level::Medium) << printEdges(cut_edges) << std::endl;
		OGDF_ASSERT(cut_edges.size() == bij.size());

		for (const FrozenPipeBijPair& pair : bij) {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
			pq.deletedEdges.restore(partner_edge_idx[pair.first]);
#pragma GCC diagnostic pop
		}
		std::pair<node, node> nodes = split(*pq.G, cut_edges, &partner_edge_idx,
				&reverse_partner_edge, pq.nodeFromIndex(u_idx), pq.nodeFromIndex(v_idx));
		OGDF_ASSERT(pq.deletedNodes.remove(nodes.first));
		OGDF_ASSERT(pq.deletedNodes.remove(nodes.second));
		pq.matchings.matchNodes(nodes.first, nodes.second);
		pq.log.lout(Logger::Level::Medium) << pq.matchings.printBijection(nodes.first) << std::endl;
		pq.log.lout(Logger::Level::Medium)
				<< "Revived matched nodes " << pq.fmtPQNode(nodes.first, false) << " and "
				<< pq.fmtPQNode(nodes.second, false) << "." << std::endl;
#ifdef OGDF_DEBUG
		PipeBij actual_edges;
		pq.matchings.getIncidentEdgeBijection(nodes.first, actual_edges);
		for (PipeBijPair& pair : actual_edges) {
			pair.first = pair.first->twin();
			pair.second = pair.second->twin();
		}
		OGDF_ASSERT(cut_edges == actual_edges);

		pq.verifyPipeBijection(nodes.first, nodes.second, bij);
#endif
		// SYNCPLAN_PROFILE_STOP("undo-contract")
	}

	std::ostream& print(std::ostream& os) const override {
		return os << "UndoContract(u=" << u_idx << ", v=" << v_idx
				  << ", bij=" << printFrozenBijection(bij) << ", u_neighs=" << u_neigh_idcs
				  << ", reverse_v_edges=" << reverse_v_edges << ")";
	}

	string name() const override {
		if (bij.size() <= 2) {
			return "UndoContractSmall";
		} else if (biconnected) {
			return "UndoContractBiconnected";
		} else {
			return "UndoContractBipartite";
		}
	}
};

SyncPlan::Result SyncPlan::contract(node u) {
	if (!matchings.isMatchedPVertex(u)) {
		return SyncPlan::Result::NOT_APPLICABLE;
	}
	node v = matchings.getTwin(u);
	if (!canContract(matchings.getPipe(u))) {
		return SyncPlan::Result::NOT_APPLICABLE;
	}

	// SYNCPLAN_PROFILE_START("contract")
	node u_bc = components.biconnectedComponent(u), v_bc = components.biconnectedComponent(v);
	bool biconnected = !components.isCutComponent(u_bc);
	log.lout(Logger::Level::High)
			<< (biconnected ? "CONTRACT/JOIN BICONNECTED" : "ENCAPSULATE AND CONTRACT/JOIN STARS")
			<< " along pipe matching " << fmtPQNode(u) << " with " << fmtPQNode(v) << std::endl;
	Logger::Indent _(&log);

	if (biconnected) {
		// SYNCPLAN_PROFILE_START("contract-bicon")
		node new_repr = nullptr;
		if (u_bc->degree() == 0 && v_bc->degree() == 0) {
			new_repr = components.bcRepr(u_bc);
			if (new_repr == u) {
				log.lout() << "Need to find a new representative instead of u for u_bc." << std::endl;
				new_repr = components.bcRepr(v_bc);
				if (new_repr == v) {
					new_repr = components.findOtherRepr(u_bc);
					log.lout() << "Using neighbour " << fmtPQNode(new_repr) << " of u." << std::endl;
				} else {
					log.lout() << "Using representative of v_bc, which is " << new_repr << "."
							   << std::endl;
				}
			}
			OGDF_ASSERT(new_repr != nullptr);
		}
		log.lout() << "Relabeling all nodes of the second biconnected component (" << v_bc << ")."
				   << std::endl;
		components.preJoin(u_bc, v_bc);
		components.makeRepr(u_bc, new_repr);
		// SYNCPLAN_PROFILE_STOP("contract-bicon")
	} else {
		// SYNCPLAN_PROFILE_START("contract-encapsulate")
		log.lout() << "Encapsulating u and v." << std::endl;
		OGDF_IF_DBG(Result result =) encapsulate(u);
		OGDF_ASSERT(result == SyncPlan::Result::SUCCESS);
		OGDF_IF_DBG(result =) encapsulate(v);
		OGDF_ASSERT(result == SyncPlan::Result::SUCCESS);
		log.lout() << "Encapsulation complete, continuing with contraction." << std::endl;
		// SYNCPLAN_PROFILE_STOP("contract-encapsulate")
	}

	log.lout(Logger::Level::Medium) << "Current incident edge bijection:" << std::endl;
	log.lout(Logger::Level::Medium) << matchings.printBijection(u) << std::endl;
	PipeBij bij;
	matchings.getIncidentEdgeBijection(u, bij);
	List<node> u_neighs;
	for (adjEntry adj : u->adjEntries) {
		u_neighs.pushBack(adj->twinNode());
	}
	auto* op = new UndoContract(u, v, bij, u_neighs, biconnected);

	log.lout() << "Joining edges." << std::endl;
	matchings.removeMatching(u);
	join(
			*G, u, v, bij,
			[&](node n) {
#ifdef OGDF_DEBUG
				deletedNodes.insert(n);
#endif
			},
			[&](edge e) { deletedEdges.hide(e); }, &op->reverse_v_edges);
	log.lout(Logger::Level::Medium) << printEdges(bij) << std::endl;

	if (!biconnected) {
		log.lout() << "Relabeling exploded stars with 2*" << u_neighs.size()
				   << " rays: " << u_neighs << std::endl;
		components.relabelExplodedStar(u_bc, v_bc, u_neighs);
	}

	pushUndoOperationAndCheck(op);
	// SYNCPLAN_PROFILE_STOP("contract")
	return SyncPlan::Result::SUCCESS;
}

}
