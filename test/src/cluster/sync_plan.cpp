/** \file
 * \brief Regression test for SyncPlan
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
#include <ogdf/basic/GraphSets.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/graph_generators.h>
#include <ogdf/basic/graph_generators/clustering.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/cluster/CconnectClusterPlanar.h>
#include <ogdf/cluster/ClusterGraph.h>
#include <ogdf/cluster/ClusterGraphAttributes.h>
#include <ogdf/cluster/HananiTutteCPlanarity.h>
#include <ogdf/cluster/sync_plan/ClusterPlanarity.h>
#include <ogdf/cluster/sync_plan/PMatching.h>
#include <ogdf/cluster/sync_plan/PipeOrder.h>
#include <ogdf/cluster/sync_plan/SyncPlan.h>
#include <ogdf/fileformats/GraphIO.h>

#include <algorithm>
#include <cstdint>
#include <exception>
#include <functional>
#include <initializer_list>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <graphs.h>
#include <resources.h>

#include <testing.h>

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
	safeForEach(G.edges, [&hes, &types, type](edge e) {
		if (types[e] == type) {
			hes.hide(e);
		}
	});
}

static std::mt19937 configRandom;

void configSyncPlan(SyncPlan& SP) {
	auto s = configRandom();
	SP.setAllowContractBBPipe(s & 1);
	s = s >> 1;
	SP.setIntersectTrees(s & 1);
	s = s >> 1;
	SP.setBatchSpqr(s & 1);
	s = s >> 1;
	bool next = (s >> 2) & 1;
	switch (s & 3) {
	case 0: {
		SP.matchings.setPipeQueue(new PipeQueueByDegree(next));
	} break;
	case 1:
		SP.matchings.setPipeQueue(new PipeQueueRandom);
		break;
	case 2:
		SP.matchings.setPipeQueue(new PipeQueueByDegreePreferContract(&SP, true, next));
		break;
	case 3:
		SP.matchings.setPipeQueue(new PipeQueueByDegreePreferContract(&SP, false, next));
		break;
	default:
		OGDF_ASSERT(false);
	}
}

using CP = SyncPlanClusterPlanarityModule;
using CconCP = CconnectClusterPlanarityModule;

static int N = 5;

go_bandit([]() {
	// ogdf::Logger::globalLogLevel(ogdf::Logger::Level::Minor);
	describe("Synchronized Planarity", []() {
		before_each([]() { ogdf::setSeed(1234567890); });
		forEachGraphDescribe(
				{GraphProperty::planar, GraphProperty::loopFree},
				[&](const Graph& G) {
					it("correctly tests random SyncPlan instances", [&G]() {
						GraphCopySimple Gcopy(G);
						planarEmbed(Gcopy);
						SyncPlan SP(&Gcopy);
						configSyncPlan(SP);
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
						randomPlanarClustering(CG, conf);

						CP cp;
						std::vector<std::pair<adjEntry, adjEntry>> augmentation;
						cp.setStoreAugmentation(&augmentation);
						AssertThat(cp.clusterPlanarEmbed(CG, Gcopy), IsTrue());
						AssertThat(CG.representsCombEmbedding(), IsTrue());
						validateGraphCopy(Gcopy);

						EdgeSet added(Gcopy);
						insertAugmentationEdges(CG, Gcopy, augmentation, &added, true, true);

						Graph::HiddenEdgeSet hes(Gcopy);
						for (edge e : added) {
							hes.hide(e);
						}
						validateGraphCopy(Gcopy);
					});
					it("correctly tests random cluster-connected graphs", [&G]() {
						ClusterGraph CG(G);
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
						randomPlanarClustering(CG, conf);
						validateGraphCopy(Gcopy);
						AssertThat(CconCP().isClusterPlanar(CG), IsTrue());
						AssertThat(CP().isClusterPlanarDestructive(CG, Gcopy), IsTrue());
						validateGraphCopy(Gcopy);
						AssertThat(CP().clusterPlanarEmbedClusterPlanarGraph(CG, Gcopy), IsTrue());
						validateGraphCopy(Gcopy);
						AssertThat(CG.representsConnectedCombEmbedding(), IsTrue());
					});

					// we need to ensure that the shared graph is connected
					// it(
					// 		"correctly tests random union graph SEFE instances",
					// 		[&G]() {
					// 			GraphCopySimple Gcopy(G);
					// 			shuffleEmbedding(Gcopy);
					// 			EdgeArray<uint8_t> types(Gcopy);
					// 			randomSEFEInstanceByUnionGraph(&Gcopy, types);
					//
					// 			Graph work;
					// 			SyncPlan SP(&Gcopy, &work, types);
					// 			AssertThat(SP.makeReduced(), IsTrue());
					// 			AssertThat(SP.solveReduced(), IsTrue());
					// 			SP.embed();
					// 			AssertThat(Gcopy.representsCombEmbedding(), IsTrue());
					// 			validateGraphCopy(Gcopy);
					// 		});

					if ((G.numberOfNodes() * 3 - 6) - G.numberOfEdges() > 0) {
						it("correctly tests random shared graph SEFE instances", [&G]() {
							GraphCopySimple Gcopy(G);
							planarEmbed(Gcopy);

							int free = (G.numberOfNodes() * 3 - 6) - G.numberOfEdges();
							free = (int)(free * max(1.0, randomDouble(0.5, 1.5)));
							int graph1 = randomNumber((int)(free * 0.1), (int)(free * 0.9));

							EdgeArray<uint8_t> types(Gcopy);
							randomSEFEInstanceBySharedGraph(&Gcopy, types, graph1, free - graph1);

							Graph work;
							SyncPlan SP(&Gcopy, &work, types);
							configSyncPlan(SP);
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
					}
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
									AssertThat(CP().clusterPlanarEmbedClusterPlanarGraph(CG, G),
											IsTrue());
								});
					}
				}
			}
		});

		describe("on grid graphs", []() {
			for (int w : {10, 11}) {
				for (int h : {9, 10, 11}) {
					for (bool loopW : {true, false}) {
						for (bool loopH : {true, false}) {
							if (loopW && loopH) {
								continue; // non-planar
							}
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
										AssertThat(CP().clusterPlanarEmbedClusterPlanarGraph(CG, G),
												Equals(!(loopH && w % 2 == 0)));
										// even width with horizontal loops is non-cplanar
									});
						}
					}
				}
			}
		});

		describe("on level graphs", []() {
			for (bool radial : {true, false}) {
				for (int n : {50, 100}) {
					for (int l : {5, 10, 25}) {
						for (double e : {1.0, 0.8, 0.6}) {
							int L = n / l;
							int M = (int)(n * e);
							if (L < 3) {
								continue;
							}
							it(string("handles ") + (radial ? "radial " : "") + "instances with "
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

										AssertThat(CP().clusterPlanarEmbed(CG, G), IsTrue());
									});
						}
					}
				}
			}
		});

		describe("interesting special cases", []() {
			it("handles a vertex enclosed by a c-connected cluster", []() {
				// create the cluster graph from Figure 9.2a, doi.org/10.15475/cpatp.2024
				Graph G;
				ClusterGraph CG(G);

				//   0   1
				// 4   5
				//   2   3
				Array<node> V;
				customGraph(G, 6,
						{
								{0, 1},
								{1, 3},
								{3, 2},
								{2, 0},
								{4, 0},
								{4, 1},
								{4, 2},
								{4, 3},
								{5, 0},
								{5, 1},
								{5, 2},
								{5, 3},
						},
						V);

				CG.createCluster({V[0], V[1], V[2], V[3]});

				AssertThat(CP().clusterPlanarEmbed(CG, G), IsFalse());
			});

			it("handles a CD-tree that becomes cyclic through reduction", []() {
				// create the cluster graph from Figure 6.7, doi.org/10.15475/cpatp.2024
				Graph G;
				ClusterGraph CG(G);

				//     3 4 5    upper
				// 0 1       2  outer
				//     6 7 8    lower
				Array<node> V;
				customGraph(G, 9,
						{{0, 4}, {0, 7}, {1, 3}, {1, 6}, {2, 5}, {2, 8}, {3, 6}, {4, 7}, {5, 8}}, V);

				cluster wrapper = CG.createEmptyCluster();
				// cluster upper =
				CG.createCluster({V[3], V[4], V[5]}, wrapper);
				// cluster lower =
				CG.createCluster({V[6], V[7], V[8]}, wrapper);

				AssertThat(CP().clusterPlanarEmbed(CG, G), IsTrue());
			});

			it("handles the 3-cluster instance where the Hanani-Tutte approach gives a false-positive",
					[]() {
						// load the cluster graph from Figure 16, doi.org/10.37236/5002
						Graph G;
						ClusterGraph CG(G);

						circulantGraph(G, 9, {1});
						std::vector<node> V {G.nodes.begin(), G.nodes.end()};

						CG.createCluster({V[0], V[3], V[6]});
						CG.createCluster({V[1], V[4], V[7]});
						CG.createCluster({V[2], V[5], V[8]});

						AssertThat(CP().isClusterPlanar(CG), IsFalse());
						AssertThat(HananiTutteCPlanarity().isCPlanar(CG, true, false,
										   HananiTutteCPlanarity::Solver::HananiTutte),
								Equals(HananiTutteCPlanarity::Verification::cPlanar));
						AssertThat(HananiTutteCPlanarity().isCPlanar(CG, true, false,
										   HananiTutteCPlanarity::Solver::HananiTutteVerify),
								Equals(HananiTutteCPlanarity::Verification::verificationFailed));
						AssertThrows(std::runtime_error, HananiTutteCPlanarity().isClusterPlanar(CG));
					});

			it("handles the 4-cluster instance where the Hanani-Tutte approach gives a false-positive",
					[]() {
						// load the cluster graph from Figure 18, doi.org/10.37236/5002
						Graph G;
						ClusterGraph CG(G);

						circulantGraph(G, 6 * 3, {1});
						std::vector<node> V {G.nodes.begin(), G.nodes.end()};

						cluster c0 = CG.createEmptyCluster();
						cluster c1 = CG.createEmptyCluster();
						cluster c2 = CG.createEmptyCluster();
						cluster c3 = CG.createEmptyCluster();
						for (int i = 0; i < 3; ++i) {
							CG.reassignNode(V[0 + i * 6], c1);
							CG.reassignNode(V[1 + i * 6], c0);
							CG.reassignNode(V[2 + i * 6], c2);
							CG.reassignNode(V[3 + i * 6], c0);
							CG.reassignNode(V[4 + i * 6], c3);
							CG.reassignNode(V[5 + i * 6], c0);
						}

						AssertThat(CP().isClusterPlanar(CG), IsFalse());
						AssertThat(HananiTutteCPlanarity().isCPlanar(CG, true, false,
										   HananiTutteCPlanarity::Solver::HananiTutte),
								Equals(HananiTutteCPlanarity::Verification::cPlanar));
						AssertThat(HananiTutteCPlanarity().isCPlanar(CG, true, false,
										   HananiTutteCPlanarity::Solver::HananiTutteVerify),
								Equals(HananiTutteCPlanarity::Verification::verificationFailed));
					});

			it("handles the instance where the Hanani-Tutte implementation gives a false-negative",
					[]() {
						// load the cluster graph from Figure 9.2b, doi.org/10.15475/cpatp.2024
						Graph G;
						ClusterGraph CG(G);
						ClusterGraphAttributes CGA(CG, ClusterGraphAttributes::all);

						std::stringstream ss {ResourceFile::get("misc/ht-err-minor.gml")->data()};
						GraphIO::readGML(CGA, CG, G, ss);

						AssertThat(CP().isClusterPlanar(CG), IsTrue());
						AssertThat(HananiTutteCPlanarity().isCPlanar(CG, true, false,
										   HananiTutteCPlanarity::Solver::HananiTutte),
								Equals(HananiTutteCPlanarity::Verification::nonCPlanarVerified));
					});
		});
	});
});
