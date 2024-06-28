/** \file
 * \brief Regression test for TODO
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

#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/graph_generators/clustering.h>
#include <ogdf/cluster/CconnectClusterPlanar.h>
#include <ogdf/cluster/sync_plan/ClusterPlanarity.h>
#include <ogdf/cluster/sync_plan/SyncPlan.h>

#include <random>

#include <graphs.h>

using namespace ogdf::sync_plan;

void validateGraphCopy(const GraphCopySimple& Gcopy) {
	const Graph& G = Gcopy.original();
	AssertThat(Gcopy.numberOfNodes(), Equals(G.numberOfNodes()));
	AssertThat(Gcopy.numberOfEdges(), Equals(G.numberOfEdges()));
	for (node n : G.nodes) {
		AssertThat(Gcopy.original(Gcopy.copy(n)), Equals(n));
		AssertThat(Gcopy.copy(n)->outdeg(), Equals(n->outdeg()));
		AssertThat(Gcopy.copy(n)->indeg(), Equals(n->indeg()));
	}
	for (edge e : G.edges) {
		AssertThat(Gcopy.original(Gcopy.copy(e)), Equals(e));
		AssertThat(Gcopy.copy(e)->source(), Equals(Gcopy.copy(e->source())));
		AssertThat(Gcopy.copy(e)->target(), Equals(Gcopy.copy(e->target())));
	}
}

void hideType(const Graph& G, const EdgeArray<uint8_t>& types, Graph::HiddenEdgeSet& hes, int type) {
	safeForEach(G.edges, [&hes, &types](edge e) {
		if (types[e] != 3) {
			hes.hide(e);
		}
	});
}

using CP = SyncPlanClusterPlanarityModule;
using CconCP = CconnectClusterPlanarityModule;

static int N = 5;

go_bandit([]() {
	// ogdf::Logger::globalLogLevel(ogdf::Logger::Level::Minor);
	describe("Synchronized Planarity", []() {
		forEachGraphDescribe(
				{GraphProperty::planar, GraphProperty::loopFree},
				[&](const Graph& G) {
					it("correctly tests random SyncPlan instances", [&G]() {
						GraphCopySimple Gcopy(G);
						planarEmbed(Gcopy);
						SyncPlan SP(&Gcopy);
						ogdf::setSeed(1234);
						randomSyncPlanInstance(SP, G.numberOfNodes() / N);
						AssertThat(SP.makeReduced(), IsTrue());
						AssertThat(SP.solveReduced(), IsTrue());
						SP.embed();
						AssertThat(Gcopy.representsCombEmbedding(), IsTrue());
						validateGraphCopy(Gcopy);
					});
					if (!isConnected(G) || !isSimpleUndirected(G)) {
						return;
					}
					it("correctly tests random cluster-planar graphs", [&G]() {
						GraphCopySimple Gcopy(G);
						planarEmbed(Gcopy);
						ClusterGraph CG(Gcopy);
						RandomClusterConfig conf;
						conf.expected_nodes(N);
						ogdf::setSeed(1234);
						randomPlanarClustering(CG, conf);

						CP cp;
						std::vector<std::pair<adjEntry, adjEntry>> augmentation;
						cp.setStoredAugmentation(&augmentation);
						AssertThat(cp.clusterPlanarEmbed(CG, Gcopy), IsTrue());
						AssertThat(CG.representsCombEmbedding(), IsTrue());
						validateGraphCopy(Gcopy);

						EdgeSet<> added(Gcopy);
						insertAugmentationEdges(CG, Gcopy, augmentation, &added, true, true);

						Graph::HiddenEdgeSet hes(Gcopy);
						for (edge e : added) {
							hes.hide(e);
						}
						validateGraphCopy(Gcopy);
					});
					it("correctly tests random cluster-connected graphs", [&G]() {
						ClusterGraph CG(G);
						ogdf::setSeed(1234);
						randomCConnectedClustering(CG, G.numberOfNodes() / N);
						AssertThat(CP().isClusterPlanar(CG), Equals(CconCP().isClusterPlanar(CG)));
					});
					it("correctly tests random cluster-planar cluster-connected graphs", [&G]() {
						GraphCopySimple Gcopy(G);
						planarEmbed(Gcopy);
						ClusterGraph CG(Gcopy);
						RandomClusterConfig conf;
						conf.expected_nodes(N);
						conf.cconnected = true;
						ogdf::setSeed(1234);
						randomPlanarClustering(CG, conf);
						validateGraphCopy(Gcopy);
						AssertThat(CconCP().isClusterPlanar(CG), IsTrue());
						AssertThat(CP().isClusterPlanarDestructive(CG, Gcopy), IsTrue());
						validateGraphCopy(Gcopy);
						AssertThat(CP().clusterPlanarEmbedClusterPlanarGraph(CG, Gcopy), IsTrue());
						validateGraphCopy(Gcopy);
						AssertThat(CG.representsConnectedCombEmbedding(), IsTrue());
					});
					it(
							"correctly tests random union graph SEFE instances",
							[&G]() {
								GraphCopySimple Gcopy(G);
								shuffleEmbedding(Gcopy);
								ogdf::setSeed(1234);
								EdgeArray<uint8_t> types(Gcopy);
								randomSEFEInstanceByUnionGraph(&Gcopy, types);

								Graph work;
								SyncPlan SP(&Gcopy, &work, types);
								AssertThat(SP.makeReduced(), IsTrue());
								AssertThat(SP.solveReduced(), IsTrue());
								SP.embed();
								AssertThat(Gcopy.representsCombEmbedding(), IsTrue());
								validateGraphCopy(Gcopy);
							},
							true); // we need to ensure that the shared graph is connected
					it("correctly tests random shared graph SEFE instances", [&G]() {
						ogdf::setSeed(1234);
						GraphCopySimple Gcopy(G);
						planarEmbed(Gcopy);

						int free = (G.numberOfNodes() * 3 - 6) - G.numberOfEdges();
						if (free <= 0) {
							return;
						}
						free = (int)(free * max(1.0, randomDouble(0.5, 1.5)));
						int graph1 = randomNumber((int)(free * 0.1), (int)(free * 0.9));

						EdgeArray<uint8_t> types(Gcopy);
						randomSEFEInstanceBySharedGraph(&Gcopy, types, graph1, free - graph1);

						Graph work;
						SyncPlan SP(&Gcopy, &work, types);
						AssertThat(SP.makeReduced(), IsTrue());
						AssertThat(SP.solveReduced(), IsTrue());
						SP.embed();
						Graph::HiddenEdgeSet hes(Gcopy);
						hideType(Gcopy, types, hes, 2);
						AssertThat(Gcopy.representsCombEmbedding(), IsTrue());
						hideType(Gcopy, types, hes, 1);
						AssertThat(Gcopy.representsCombEmbedding(), IsTrue());
						validateGraphCopy(Gcopy);
					});
				},
				GraphSizes(), 3);

		describe("on random cluster graphs", []() {
			for (int c : {5, 10, 20, 30}) {
				for (int n : {10, 20, 50}) {
					for (double e : {2.0, 2.5, 3.0}) {
						it("handles instances with " + to_string(n) + " nodes and "
										+ to_string((int)(e * n)) + " edges in " + to_string(c)
										+ " clusters",
								[&]() {
									Graph G;
									ClusterGraph CG(G);
									randomClusterPlanarGraph(G, CG, c, n, (int)(e * n));
									AssertThat(SyncPlanClusterPlanarityModule()
													   .clusterPlanarEmbedClusterPlanarGraph(CG, G),
											IsTrue());
								});
					}
				}
			}
		});

		describe("on grid graphs", []() {
			for (int w : {20, 21}) {
				for (int h : {19, 20, 21}) {
					for (bool loopW : {true, false}) {
						for (bool loopH : {true, false}) {
							it("handles width " + to_string(w) + " and height " + to_string(h) + " "
											+ (loopW ? "with" : "without") + " horizontal and "
											+ (loopH ? "with" : "without") + " vertical loops",
									[&]() {
										Graph G;
										gridGraph(G, w, h, loopW, loopH);
										ClusterGraph CG(G);
										cluster c1 = CG.createEmptyCluster();
										cluster c2 = CG.createEmptyCluster();
										for (node n : G.nodes) {
											if (n->index() % 2 == 0) {
												CG.reassignNode(n, c1);
											} else {
												CG.reassignNode(n, c2);
											}
										}
										AssertThat(
												SyncPlanClusterPlanarityModule()
														.clusterPlanarEmbedClusterPlanarGraph(CG, G),
												IsTrue());
										if (loopH && loopW) {
											AssertThat(isRegular(G, 4), IsTrue());
										}
									});
						}
					}
				}
			}
		});

		describe("on level graphs", []() {
			for (bool radial : {true, false}) {
				for (int n : {50, 100, 500}) {
					for (int l : {10, 25, 50}) {
						for (double e : {1.0, 0.8, 0.6}) {
							int L = n / l;
							int M = (int)(n * e);
							it(string("handles ") + (radial ? "radial " : "") + " instances with "
											+ to_string(n) + " nodes and " + to_string(M)
											+ " edges on " + to_string(L) + " levels",
									[&]() {
										Graph LG;
										std::vector<std::vector<node>> emb;
										randomProperMaximalLevelPlaneGraph(LG, emb, n, L, radial);
										pruneEdges(LG, M, 1);
										if (!radial) {
											node u = nullptr;
											for (auto& lvl : emb) {
												node v = LG.newNode();
												lvl.push_back(v);
												if (u != nullptr) {
													LG.newEdge(u, v);
												}
												u = v;
											}
										}

										Graph G;
										ClusterGraph CG(G);
										EdgeArray<node> map(G);
										reduceLevelPlanarityToClusterPlanarity(LG, emb, G, CG, map);

										CP cp;
										AssertThat(cp.clusterPlanarEmbed(CG, G), IsTrue());
									});
						}
					}
				}
			}
		});
	});
});
