/** \file
 * \brief Tests for ogdf::DualGraph.
 *
 * \author Mirko Wagner, Tilo Wiedera, Max Ilsen
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

#include <ogdf/basic/DualGraph.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/graph_generators.h>

#include <graphs.h>
#include <testing.h>

//! Tests consistency of a DynamicDualGraph \p dual and its corresponding primal
//! embedding \p emb as well as the primal graph &graph.
void describeDualGraph(const DynamicDualGraph& dual, const ConstCombinatorialEmbedding& emb,
		const Graph& graph) {
	it("returns its primal embedding", [&] { AssertThat(&dual.getPrimalEmbedding(), Equals(&emb)); });

	it("returns its primal graph", [&] { AssertThat(&dual.getPrimalGraph(), Equals(&graph)); });

	it("has a matching number of nodes, faces, and edges", [&] {
		AssertThat(dual.numberOfFaces(), Equals(graph.numberOfNodes()));
		AssertThat(dual.getGraph().numberOfNodes(), Equals(emb.numberOfFaces()));
		AssertThat(dual.getGraph().numberOfEdges(), Equals(graph.numberOfEdges()));
	});

	it("maps primal faces to dual nodes", [&] {
		for (face f : emb.faces) {
			node v = dual.dualNode(f);
			AssertThat(dual.primalFace(v), Equals(f));
			AssertThat(v->degree(), Equals(f->size()));
		}
	});

	it("maps primal nodes to dual faces", [&] {
		for (node v : graph.nodes) {
			face f = dual.dualFace(v);
			AssertThat(dual.primalNode(f), Equals(v));
			AssertThat(f->size(), Equals(v->degree()));
		}
	});

	it("maps edges and faces", [&] {
		for (edge e : graph.edges) {
			edge g = dual.dualEdge(e);
			AssertThat(g, !Equals(e));
			AssertThat(dual.primalEdge(g), Equals(e));

			face f = dual.primalFace(g->source());
			AssertThat(f, Equals(emb.rightFace(e->adjSource())));
			AssertThat(f, Equals(emb.leftFace(e->adjTarget())));

			f = dual.primalFace(g->target());
			AssertThat(f, Equals(emb.leftFace(e->adjSource())));
			AssertThat(f, Equals(emb.rightFace(e->adjTarget())));
		}
	});
}

//! Creates a DualGraph of \p graph and runs several tests on it.
void describeDualGraph(const Graph& graph) {
	if (graph.numberOfEdges() < 1) {
		return;
	}

	GraphCopy copy;
	CombinatorialEmbedding emb;
	std::unique_ptr<DynamicDualGraph> dual;

	describe("initialization", [&] {
		copy.init(graph);
		planarEmbed(copy);
		emb.init(copy);
		dual.reset(new DynamicDualGraph {emb});
		describeDualGraph(*dual, emb, copy);
	});

	if (graph.numberOfEdges() > 2) {
		describe("split and unsplit edges", [&] {
			copy.init(graph);
			planarEmbed(copy);
			emb.init(copy);
			dual.reset(new DynamicDualGraph {emb});
			internal::GraphObjectContainer<ogdf::EdgeElement>::iterator it {copy.edges.begin()};
			for (int i {0}; i < min(10, copy.numberOfEdges()); ++i) {
				edge e {*it};
				it++;
				edge eOut {dual->splitPrimal(e)};
				dual->unsplitPrimal(e, eOut);
			}
			describeDualGraph(*dual, emb, copy);
		});

		describe("split and unsplit faces", [&] {
			copy.init(graph);
			planarEmbed(copy);
			emb.init(copy);
			dual.reset(new DynamicDualGraph {emb});
			internal::GraphObjectContainer<ogdf::FaceElement>::iterator it {emb.faces.begin()};
			for (int i {0}; i < min(10, emb.numberOfFaces()); ++i) {
				face f {*it};
				it++;
				adjEntry firstAdj {f->firstAdj()};
				edge e {dual->splitFacePrimal(firstAdj, firstAdj->clockwiseFacePred())};
				dual->joinFacesPrimal(e);
			}
			describeDualGraph(*dual, emb, copy);
		});

		describe("split and join nodes", [&] {
			copy.init(graph);
			planarEmbed(copy);
			emb.init(copy);
			dual.reset(new DynamicDualGraph {emb});
			internal::GraphObjectContainer<ogdf::NodeElement>::iterator it {copy.nodes.begin()};
			for (int i {0}; i < min(10, copy.numberOfNodes()); ++i) {
				node v {*it};
				it++;
				adjEntry firstAdj {v->firstAdj()};
				dual->splitNodePrimal(firstAdj, firstAdj->cyclicPred());
				dual->contractPrimal(firstAdj->cyclicPred()->theEdge(), true);
				// TODO If contractPrimal-param keepSelfLoops is set to false,
				// this fails.
			}
			describeDualGraph(*dual, emb, copy);
		});

		describe("add and remove degree-1-nodes", [&] {
			copy.init(graph);
			planarEmbed(copy);
			emb.init(copy);
			dual.reset(new DynamicDualGraph {emb});
			internal::GraphObjectContainer<ogdf::NodeElement>::iterator it {copy.nodes.begin()};
			for (int i {0}; i < min(10, copy.numberOfNodes()); ++i) {
				node v {*it};
				it++;
				adjEntry firstAdj {v->firstAdj()};
				node newNode {copy.newNode()};
				if (randomNumber(0, 1) == 0) {
					dual->addEdgeToIsolatedNodePrimal(firstAdj, newNode);
				} else {
					dual->addEdgeToIsolatedNodePrimal(newNode, firstAdj);
				}
				dual->removeDeg1Primal(newNode);
			}
			describeDualGraph(*dual, emb, copy);
		});
	}
}

go_bandit([]() {
	describe("DualGraph", [] {
		forEachGraphDescribe({GraphProperty::planar, GraphProperty::connected},
				[&](const Graph& graph) { describeDualGraph(graph); });
	});
});
