/** \file
 * \brief Tests for functionality from ogdf/graphalg/Triconnectivity.h
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

#include <ogdf/basic/Array.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/graph_generators.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/graphalg/Triconnectivity.h>

#include <algorithm>
#include <functional>
#include <iterator>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <graphs.h>

#include <testing.h>

#define MAX_SIZE 1000
#define SKIP_PARALLEL true
#define SKIP_COMP false
#define SKIP_SPLIT false

struct RandomSkeleton : public Graph {
	RandomSkeleton* m_parent = nullptr;
	edge m_parent_edge = nullptr;
	EdgeArrayP<RandomSkeleton> m_virtual_edges;

	RandomSkeleton() : m_virtual_edges(*this) { }

	RandomSkeleton(Triconnectivity::CompType type, int edge_cnt, RandomSkeleton* parent,
			double virt_prob, double size_fact)
		: m_parent(parent), m_virtual_edges(*this) {
		OGDF_ASSERT(edge_cnt >= 3);
		if (type == ogdf::Triconnectivity::CompType::triconnected) {
			randomPlanarTriconnectedGraph(*this, edge_cnt / 3, edge_cnt);
		} else if (type == ogdf::Triconnectivity::CompType::bond) {
			node n1 = newNode();
			node n2 = newNode();
			for (int i = 0; i < edge_cnt; ++i) {
				newEdge(n1, n2);
			}
		} else {
			OGDF_ASSERT(type == ogdf::Triconnectivity::CompType::polygon);
			node n1 = newNode();
			node n = n1;
			for (int i = 0; i < edge_cnt - 1; ++i) {
				n = newEdge(n, newNode())->target();
			}
			newEdge(n, n1);
		}
		if (m_parent != nullptr) {
			m_parent_edge = chooseEdge();
		}
		for (edge e : edges) {
			if (e == m_parent_edge || randomDouble(0, 1) >= virt_prob) {
				continue;
			}
			int cnt = (int)(edge_cnt * size_fact * randomDouble(0.8, 1.2));
			if (cnt < 3) {
				continue;
			}
			m_virtual_edges[e] = std::make_unique<RandomSkeleton>(
					static_cast<Triconnectivity::CompType>(randomNumber(0, 2)), cnt, this,
					virt_prob * 0.8, size_fact);
		}
	}

	void createGraph(Graph& G, edge parent = nullptr) {
		NodeArray<node> copyN(*this);
		EdgeArray<edge> copyE(*this);
		G.insert(*this, copyN, copyE);
		for (edge e : this->edges) {
			if (m_virtual_edges[e]) {
				m_virtual_edges[e]->createGraph(G, copyE[e]);
			}
		}
		if (parent) {
			G.contract(G.newEdge(parent->source(), copyN[m_parent_edge->source()]));
			G.contract(G.newEdge(parent->target(), copyN[m_parent_edge->target()]));
			G.delEdge(copyE[m_parent_edge]);
			G.delEdge(parent);
		}
	}

	RandomSkeleton& makeVirtual(edge e) {
		auto& ptr = m_virtual_edges[e];
		OGDF_ASSERT(!ptr);
		OGDF_ASSERT(e != m_parent_edge);
		ptr = std::make_unique<RandomSkeleton>();
		ptr->m_parent = this;
		return *ptr;
	}
};

void assertSameUndirected(const GraphCopySimple& GC, edge e, node orig_source, node orig_target) {
	if (GC.original(e->source()) == orig_source) {
		AssertThat(GC.original(e->target()), Equals(orig_target));
	} else {
		AssertThat(GC.original(e->source()), Equals(orig_target));
		AssertThat(GC.original(e->target()), Equals(orig_source));
	}
}

void assertSameUndirected(const GraphCopySimple& GC, edge e, edge orig) {
	assertSameUndirected(GC, e, orig->source(), orig->target());
}

int checkEdgeCounts(const Graph& G, const Triconnectivity& T) {
	AssertThat(T.checkComp(), IsTrue());
	const GraphCopySimple& GC = *dynamic_cast<const GraphCopySimple*>(T.m_pG);
	AssertThat(GC.numberOfNodes(), Equals(G.numberOfNodes()));
	for (edge e : G.edges) {
		assertSameUndirected(GC, GC.copy(e), e);
	}
	int virt_edges = 0, real_edges = 0, comps = 0;
	for (int i = 0; i < T.m_numComp; ++i) {
		using namespace std;
		auto& comp = T.m_component[i];
		if (comp.m_edges.empty()) {
			continue;
		}
		comps++;
		// cout << endl << "Comp " << i << " " << comp.m_type << endl;
		for (edge e : comp.m_edges) {
			if (GC.isDummy(e)) {
				// cout << "\tvirt (";
				virt_edges++;
			} else {
				// cout << "\treal (";
				real_edges++;
				assertSameUndirected(GC, e, GC.original(e));
			}
			// cout << GC.original(e->source()) << " " << GC.original(e->target()) << ")";
			// cout << endl;
		}
		// cout << endl;
		{
			GraphCopySimple GCcomp;
			GCcomp.setOriginalGraph(GC);
			NodeArray<node> nodeMap(GC);
			EdgeArray<edge> edgeMap(GC);
			GCcomp.insert(GC.nodes, comp.m_edges, nodeMap, edgeMap);
			for (node n : std::vector<node> {GCcomp.nodes.begin(), GCcomp.nodes.end()}) {
				if (n->degree() < 1) {
					GCcomp.delNode(n);
				}
			}
			// for (edge e : GCcomp.edges) {
			// 	cout << "(" << GC.original(GCcomp.original(e->source())) << " "
			// 		 << GC.original(GCcomp.original(e->target())) << ") ";
			// }
			// cout << endl;
			AssertThat(GCcomp.numberOfEdges(), Equals(comp.m_edges.size()));
			if (GCcomp.numberOfNodes() > 3) {
				AssertThat(isBiconnected(GCcomp), IsTrue());
			} else {
				AssertThat(isConnected(GCcomp), IsTrue());
			}
			if (comp.m_type == ogdf::Triconnectivity::CompType::bond) {
				AssertThat(GCcomp.numberOfNodes(), Equals(2));
				AssertThat(isRegular(GCcomp, comp.m_edges.size()), IsTrue());

				node s = GC.original(comp.m_edges.front()->source());
				node t = GC.original(comp.m_edges.front()->target());
				for (edge e : comp.m_edges) {
					assertSameUndirected(GC, e, s, t);
				}
			} else if (comp.m_type == ogdf::Triconnectivity::CompType::triconnected) {
				AssertThat(isTriconnectedPrimitive(GCcomp), IsTrue());
			} else {
				AssertThat(comp.m_type, Equals(ogdf::Triconnectivity::CompType::polygon));
				AssertThat(isRegular(GCcomp, 2), IsTrue());
			}
		}
	}
	AssertThat(real_edges, Equals(G.numberOfEdges()));
	AssertThat(virt_edges, Equals(comps * 2 - 2));
	AssertThat(T.m_pG->numberOfEdges(), Equals(G.numberOfEdges() + comps - 1));
	return comps;
}

void checkBicon(const Graph& G, const string& graphName, bool checkIsTric) {
	it(
			"computes components",
			[&G, checkIsTric, graphName]() {
				Triconnectivity T(G);
				int comps = checkEdgeCounts(G, T);
				if (checkIsTric) {
					AssertThat(comps, Equals(1));
				} else {
					AssertThat(comps, IsGreaterThan(1));
				}
			},
			SKIP_COMP);

	it(
			"computes split pairs",
			[&G, checkIsTric]() {
				bool isTric = true;
				node s1 = nullptr;
				node s2 = nullptr;
				Triconnectivity _(G, isTric, s1, s2);
				AssertThat(isTric, Equals(checkIsTric));
				if (checkIsTric) {
					AssertThat(s1, IsNull());
					AssertThat(s2, IsNull());
				} else {
					AssertThat(s1, !IsNull());
					AssertThat(s2, !IsNull());
					AssertThat(s1, !Equals(s2));
					OGDF_ASSERT(s1->graphOf() == &G);
					OGDF_ASSERT(s2->graphOf() == &G);
				}
			},
			SKIP_SPLIT);
}

go_bandit([]() {
	describe("Triconnectivity", []() {
		describe("for triconnected graphs", []() {
			forEachGraphDescribe(
					{GraphProperty::triconnected, GraphProperty::simple},
					[](const Graph& G) {
						it(
								"computes components",
								[&G]() {
									Triconnectivity T(G);
									int comps = checkEdgeCounts(G, T);
									AssertThat(comps, Equals(1));
									AssertThat(T.m_component[0].m_type,
											Equals(Triconnectivity::CompType::triconnected));
									AssertThat(T.m_component[0].m_edges.size(),
											Equals(G.numberOfEdges()));
									AssertThat(T.m_pG->numberOfEdges(), Equals(G.numberOfEdges()));
								},
								SKIP_COMP);

						it(
								"computes split pairs",
								[&G]() {
									bool isTric = false;
									node s1 = nullptr;
									node s2 = nullptr;
									Triconnectivity _(G, isTric, s1, s2);
									AssertThat(isTric, IsTrue());
									AssertThat(s1, IsNull());
									AssertThat(s2, IsNull());
								},
								SKIP_SPLIT);
					},
					GraphSizes(), 3, MAX_SIZE);
		});

		describe(
				"for triconnected graphs with parallel edges",
				[]() {
					forEachGraphDescribe(
							{GraphProperty::triconnected, GraphProperty::loopFree},
							[](const Graph& G) {
								std::vector<edge> edges {G.edges.begin(), G.edges.end()};
								int parallels = 0;
								for (int i = 0; i < G.edges.size(); ++i) {
									if (randomDouble(0, 1) < 0.1) {
										parallels++;
										edge e = edges[i];
										for (int j = 0; j < randomNumber(1, 10); j++) {
											edges.push_back(e);
										}
									}
								}
								if (parallels < 1) {
									parallels++;
									edge e = edges[randomNumber(0, G.numberOfEdges() - 1)];
									for (int j = 0; j < randomNumber(1, 10); j++) {
										edges.push_back(e);
									}
								}
								auto rng = std::default_random_engine {(unsigned int)randomSeed()};
								std::shuffle(edges.begin(), edges.end(), rng);

								Graph G2;
								NodeArray<node> origN(G, nullptr);
								EdgeArray<edge> origE(G, nullptr);
								G2.insert<internal::GraphObjectContainer<NodeElement>::iterator, //
										std::vector<edge>::iterator, //
										false, false, true>( //
										G.nodes.begin(), G.nodes.end(), //
										edges.begin(), edges.end(), //
										origN, origE);
								makeIndicesNonContinuous(G2, 0.2);
								for (edge e : G2.edges) {
									if (randomDouble(0, 1) < 0.3) {
										G2.reverseEdge(e);
									}
								}

								EdgeArray<SListPure<edge>> para_edges(G2);
								getParallelFreeUndirected(G2, para_edges);

								it(
										"computes components",
										[&G2, &G, &para_edges, parallels]() {
											Triconnectivity T(G2);
											const GraphCopySimple& GC =
													*dynamic_cast<const GraphCopySimple*>(T.m_pG);
											int comps = checkEdgeCounts(G, T);
											AssertThat(comps, Equals(1 + parallels));
											int found_ps = 0, found_rs = 0;
											for (int i = 0; i < T.m_numComp; ++i) {
												auto& comp = T.m_component[i];
												if (comp.m_type
														== Triconnectivity::CompType::triconnected) {
													found_rs++;
													AssertThat(comp.m_edges.size(),
															Equals(G.numberOfEdges()));
												} else {
													AssertThat(comp.m_type,
															Equals(Triconnectivity::CompType::bond));
													found_ps++;
													edge orig = nullptr;
													edge virt = nullptr;
													for (edge e : comp.m_edges) {
														if (GC.isDummy(e)) {
															AssertThat(virt, IsNull());
															virt = e;
														} else if (para_edges[GC.original(e)].size()
																> 0) {
															AssertThat(orig, IsNull());
															AssertThat(comp.m_edges.size(),
																	Equals(para_edges[GC.original(e)]
																					.size()
																			+ 2));
															orig = e;
														}
													}
													AssertThat(orig, !IsNull());
													AssertThat(virt, !IsNull());
												}
											}
											AssertThat(found_ps, Equals(parallels));
											AssertThat(found_rs, Equals(1));
										},
										SKIP_COMP);

								it(
										"computes split pairs",
#ifdef OGDF_DEBUG
										[&G2, &G]() {
#else
										[&G2]() {
#endif
											bool isTric = true;
											node s1 = nullptr;
											node s2 = nullptr;
											Triconnectivity _(G2, isTric, s1, s2);
											AssertThat(isTric, IsFalse());
											AssertThat(s1, !IsNull());
											AssertThat(s2, !IsNull());
											AssertThat(s1, !Equals(s2));
											OGDF_ASSERT(s1->graphOf() == &G);
											OGDF_ASSERT(s2->graphOf() == &G);
										},
										SKIP_SPLIT);
							},
							GraphSizes(), 3, MAX_SIZE);

					it("correctly handles the instance from GH issue #207", []() {
						RandomSkeleton Po;
						customGraph(Po, 2, {{0, 1}, {0, 1}, {0, 1}});

						auto& Sl = Po.makeVirtual(Po.edges.head());
						auto& Sr = Po.makeVirtual(Po.edges.tail());
						customGraph(Sl, 4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}});
						Sl.m_parent_edge = Sl.edges.head();
						customGraph(Sr, 4, {{0, 1}, {1, 2}, {2, 3}, {0, 3}});
						Sr.m_parent_edge = Sr.edges.head();

						auto& Pr = Sr.makeVirtual(*std::next(Sr.edges.begin(), 2));
						customGraph(Pr, 2, {{0, 1}, {0, 1}, {1, 0}});
						Pr.m_parent_edge = Pr.edges.head();

						auto& Pl = Sl.makeVirtual(*std::next(Sr.edges.begin(), 2));
						customGraph(Pl, 2, {{0, 1}, {0, 1}, {0, 1}});
						Pl.m_parent_edge = Pl.edges.head();

						Graph G;
						Po.createGraph(G);

						Triconnectivity T(G);
						int comps = checkEdgeCounts(G, T);
						AssertThat(comps, Equals(4));
					});
				},
				SKIP_PARALLEL);

		describe("for biconnected graphs", []() {
			forEachGraphDescribe(
					{GraphProperty::biconnected, GraphProperty::loopFree},
					[](const Graph& G, const std::string& graphName,
							const std::set<GraphProperty>& props) {
						// // it seems that all our biconnected graphs are already parallel free
						// GraphCopySimple GC(G);
						// makeParallelFree(GC);
						// if (GC.numberOfEdges() != G.numberOfEdges()) {
						// 	describe("without parallels", [&G, &GC, &graphName]() {
						// 		AssertThat(isSimpleUndirected(GC), IsTrue());
						// 		bool checkIsTric = isTriconnectedPrimitive(G);
						// 		checkBicon(GC, graphName, checkIsTric);
						// 	});
						// }
						AssertThat(isSimpleUndirected(G), IsTrue());
						bool checkIsTric = isTriconnectedPrimitive(G);
						checkBicon(G, graphName, checkIsTric);
					},
					GraphSizes(), 3, MAX_SIZE);
		});

		describe("for two joined rigids", []() {
			forEachGraphDescribe(
					{GraphProperty::triconnected, GraphProperty::simple},
					[](const Graph& G1) {
						Graph G;
						NodeArray<node> copyN(G1);
						EdgeArray<edge> copyE(G1);
						G.insert(G1, copyN, copyE);

						node s1 = G.edges.head()->source();
						node t1 = G.edges.head()->target();

						RandomSkeleton rigid2;
						randomPlanarTriconnectedGraph(rigid2, 20, 60);
						rigid2.m_parent_edge = rigid2.edges.head();
						rigid2.createGraph(G, G.edges.head());

						it(
								"computes components",
								[&G, &G1, &rigid2]() {
									Triconnectivity T(G);
									int comps = checkEdgeCounts(G, T);
									AssertThat(comps, Equals(2));

									Triconnectivity::CompStruct *comp1 = nullptr, *comp2 = nullptr;
									for (int i = 0; i < T.m_numComp; ++i) {
										if (T.m_component[i].m_edges.empty()) {
											continue;
										}
										if (!comp1) {
											comp1 = &T.m_component[i];
										} else {
											AssertThat(comp2, IsNull());
											comp2 = &T.m_component[i];
										}
									}
									AssertThat(comp2, !IsNull());

									AssertThat(comp1->m_type,
											Equals(Triconnectivity::CompType::triconnected));
									AssertThat(comp2->m_type,
											Equals(Triconnectivity::CompType::triconnected));
									if (comp1->m_edges.size() == G1.numberOfEdges()) {
										AssertThat(comp2->m_edges.size(),
												Equals(rigid2.numberOfEdges()));
									} else {
										AssertThat(comp1->m_edges.size(),
												Equals(rigid2.numberOfEdges()));
										AssertThat(comp2->m_edges.size(), Equals(G1.numberOfEdges()));
									}
								},
								SKIP_COMP);

						it(
								"computes split pairs",
								[&G, s1, t1]() {
									bool isTric = false;
									node s2 = nullptr;
									node t2 = nullptr;
									Triconnectivity _(G, isTric, s2, t2);
									AssertThat(isTric, IsFalse());
									if (s1 == s2) {
										AssertThat(t2, Equals(t1));
									} else {
										AssertThat(s2, Equals(t1));
										AssertThat(t2, Equals(s1));
									}
								},
								SKIP_SPLIT);
					},
					GraphSizes(), 3, MAX_SIZE);
		});

		// TODO check that random SPQR tree structure matches
		// TODO check that separation pair occurs in random SPQR tree
		// TODO fix and reenable parallel edges case
		// {
		// 	RandomSkeleton skel(static_cast<Triconnectivity::CompType>(randomNumber(0, 2)), 20,
		// 			nullptr, 0.3, 0.8);
		// 	Graph G;
		// 	skel.createGraph(G, nullptr);
		// 	std::cout << G.numberOfNodes() << " " << G.numberOfEdges() << std::endl;
		// }
	});
});
