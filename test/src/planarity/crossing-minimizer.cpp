/** \file
 * \brief Tests for various crossing minimization modules.
 *
 * \author Tilo Wiedera
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
#include <ogdf/basic/graph_generators.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/energybased/FMMMLayout.h>
#include <ogdf/planarity/FixedEmbeddingInserter.h>
#include <ogdf/planarity/MultiEdgeApproxInserter.h>
#include <ogdf/planarity/PlanarizerChordlessCycle.h>
#include <ogdf/planarity/PlanarizerMixedInsertion.h>
#include <ogdf/planarity/PlanarizerStarReinsertion.h>
#include <ogdf/planarity/SubgraphPlanarizer.h>
#include <ogdf/planarity/VariableEmbeddingInserter.h>
#include <ogdf/planarity/VariableEmbeddingInserterDyn.h>

#include <set>

#include <resources.h>

/**
 * Verifies that \p graph resembles a planarization of the original graph.
 *
 * \param graph a supposed planarization to be verified
 * \param cost a pointer to the cost of each edge in the original graph
 * \return the weighted crossing number of the given planarization
 */
int verifyCrossings(const GraphCopy& graph, const EdgeArray<int>* cost) {
	int result = 0;

	const Graph& original = graph.original();
	int numberOfDummies = graph.numberOfNodes() - original.numberOfNodes();
	AssertThat(graph.numberOfEdges() - original.numberOfEdges(), Equals(2 * numberOfDummies));

	int dummyCounter = 0;
	for (node v : graph.nodes) {
		if (graph.isDummy(v)) {
			dummyCounter++;

			AssertThat(v->degree(), Equals(4));
			AssertThat(v->indeg(), Equals(2));

			edge e = graph.original(v->firstAdj()->theEdge());
			edge f = graph.original(v->lastAdj()->theEdge());
			edge e2 = graph.original(v->lastAdj()->cyclicPred()->theEdge());
			edge f2 = graph.original(v->firstAdj()->cyclicSucc()->theEdge());
			AssertThat(e, Equals(e2));
			AssertThat(f, Equals(f2));

			List<edge> inEdges;
			v->inEdges(inEdges);
			AssertThat(graph.original(inEdges.front()), !Equals(graph.original(inEdges.back())));

			AssertThat(e, !IsNull());
			AssertThat(f, !IsNull());
			result += cost ? (*cost)[e] * (*cost)[f] : 1;
		}
	}

	AssertThat(dummyCounter, Equals(numberOfDummies));

	for (edge e : graph.edges) {
		node s = e->source();
		node t = e->target();

		AssertThat(graph.isDummy(e), IsFalse());

		if (!graph.isDummy(s)) {
			AssertThat(s, Equals(graph.copy(graph.original(e)->source())));
		}

		if (!graph.isDummy(t)) {
			AssertThat(t, Equals(graph.copy(graph.original(e)->target())));
		}
	}

	return result;
}

/**
 * Tests a planarization algorithm on a single instance.
 *
 * \param cmm an algorithm to be tested
 * \param graph a graph that should be planarized
 * \param expected the crossing number of the input graph
 * \param isOptimal whether the algorithm is supposed to yield an optimal solution for this instance
 * \param cost costs of all edges. If \c nullptr is given each edge is assumed to have cost 1
 */
void testComputation(CrossingMinimizationModule& cmm, const Graph& graph, int expected,
		bool isOptimal, const EdgeArray<int>* cost = nullptr) {
	using ReturnType = CrossingMinimizationModule::ReturnType;

	PlanRep planRep(graph);
	planRep.initCC(0);
	int actual(17); // an arbitrary nonzero number
	ReturnType result = cmm.call(planRep, 0, actual, cost);
	if (isOptimal) {
		AssertThat(result, Equals(ReturnType::Optimal));
	} else {
		AssertThat(result,
				Equals(ReturnType::Optimal) || Equals(ReturnType::Feasible)
						|| Equals(ReturnType::TimeoutFeasible));
	}

	if (isOptimal) {
		AssertThat(actual, Equals(expected));
	} else {
		AssertThat(actual, !IsLessThan(expected));
	}

	bool planar = isPlanar(planRep);

	// optimal algorithms don't need to return planarizations
	if (!isOptimal) {
		AssertThat(planar, IsTrue());
	}

	if (planar) {
		AssertThat(verifyCrossings(planRep, cost), Equals(actual));
	}
	if (planar && isLoopFree(graph)) {
		AssertThat(isLoopFree(planRep), IsTrue());
	}
}

/**
 * Tests an ogdf::CrossingMinimizationModule \p cmm for correctness.
 *
 * \param cmm an algorithm to be tested
 * \param title a human-readable title of the algorithm
 * \param isOptimal whether the algorithm is optimal
 */
void testModule(CrossingMinimizationModule& cmm, const std::string title, bool isOptimal) {
	describe(title, [&]() {
		Graph graph;

		it("works on a K4", [&]() {
			completeGraph(graph, 4);
			testComputation(cmm, graph, 0, isOptimal);
		});

		it("works on a K5", [&]() {
			completeGraph(graph, 5);
			testComputation(cmm, graph, 1, isOptimal);
		});

		it("works on a K6", [&]() {
			completeGraph(graph, 6);
			testComputation(cmm, graph, 3, isOptimal);
		});

		it("works on a K3,3", [&]() {
			completeBipartiteGraph(graph, 3, 3);
			testComputation(cmm, graph, 1, isOptimal);
		});

		it("works on a K4,3", [&]() {
			completeBipartiteGraph(graph, 4, 3);
			testComputation(cmm, graph, 2, isOptimal);
		});

		it("works on a K4,4", [&]() {
			completeBipartiteGraph(graph, 4, 4);
			testComputation(cmm, graph, 4, isOptimal);
		});

		it("works on a petersen graph", [&]() {
			petersenGraph(graph, 5, 2);
			testComputation(cmm, graph, 2, isOptimal);
		});

		it("works on a generalized petersen graph (9,2)", [&]() {
			petersenGraph(graph, 9, 2);
			testComputation(cmm, graph, 3, isOptimal);
		});

		it("works on a weighted K3,3", [&]() {
			completeBipartiteGraph(graph, 3, 3);

			EdgeArray<int> cost(graph, 2);
			testComputation(cmm, graph, 4, isOptimal, &cost);

			cost[graph.chooseEdge()] = 1;
			testComputation(cmm, graph, 2, isOptimal, &cost);
		});

		// TODO test forbidden edges ?

		if (isOptimal) {
#ifdef OGDF_USE_ASSERT_EXCEPTIONS
			// optimal algorithms should throw exceptions on non-pre-processed instances

			it("aborts if the graph contains self-loops", [&]() {
				completeGraph(graph, 5);
				node v = graph.chooseNode();
				graph.newEdge(v, v);
				AssertThrows(AssertionFailed, testComputation(cmm, graph, 1, true));
			});

			it("aborts if the graph contains parallel edges", [&]() {
				completeGraph(graph, 5);
				graph.newEdge(graph.firstNode(), graph.lastNode());
				AssertThrows(AssertionFailed, testComputation(cmm, graph, 1, true));
			});

			it("aborts if the graph contains nodes with degree 2", [&]() {
				completeGraph(graph, 5);
				node v = graph.newNode();
				graph.newEdge(graph.chooseNode(), v);
				graph.newEdge(graph.chooseNode(), v);
				AssertThrows(AssertionFailed, testComputation(cmm, graph, 1, true));
			});

			it("aborts if the graph isn't biconnected", [&]() {
				completeGraph(graph, 5);
				List<node> nodes = {graph.chooseNode()};

				for (int i = 0; i < 4; i++) {
					nodes.pushBack(graph.newNode());
				}

				for (node v : nodes) {
					for (node w : nodes) {
						if (w->index() < v->index()) {
							graph.newEdge(v, w);
						}
					}
				}

				AssertThrows(AssertionFailed, testComputation(cmm, graph, 1, true));
			});
#endif
		} else {
			// we assume non-optimal algorithms to be faster
			// some algorithms require biconnectivity!

			it("works on a generalized petersen graph (15,3)", [&]() {
				petersenGraph(graph, 15, 3);
				testComputation(cmm, graph, 5, isOptimal);
			});

			it("works on a K10", [&]() {
				completeGraph(graph, 10);
				testComputation(cmm, graph, 60, false);
			});

			it("works on a random simple biconnected graph with 50 nodes", [&]() {
				randomBiconnectedGraph(graph, 50, 70);
				makeSimpleUndirected(graph);
				testComputation(cmm, graph, -1, false);
			});
		}
	});
}

template<typename EdgeInserter>
void setRemoveReinsert(EdgeInserter& edgeInserter, RemoveReinsertType rrType) {
	edgeInserter.removeReinsert(rrType);
}

template<>
void setRemoveReinsert(MultiEdgeApproxInserter& edgeInserter, RemoveReinsertType rrType) {
	edgeInserter.removeReinsertVar(rrType);
	edgeInserter.removeReinsertFix(rrType);
}

/**
 * Test the ogdf::SubgraphPlanarizer with a specific type of edge remove-reinsert post-processing.
 */
template<typename EdgeInserter>
void testSPRRType(SubgraphPlanarizer& heuristic, EdgeInserter* edgeInserter,
		RemoveReinsertType rrType, const std::string name) {
	auto performTest = [&]() {
		setRemoveReinsert(*edgeInserter, rrType);
		heuristic.permutations(1);
		testModule(heuristic, "single run", false);
		heuristic.permutations(4);
		testModule(heuristic, "4 permutations", false);
	};

	string title = "remove-reinsert: " + name;
	describe(title, performTest);
}

/**
 * Test the ogdf::SubgraphPlanarizer with a specific ogdf::EdgeInsertionModule.
 */
template<typename EdgeInserter>
void testSPEdgeInserter(EdgeInserter* edgeInserter, const std::string name) {
	describe("edge insertion: " + name, [&]() {
		SubgraphPlanarizer heuristic;
		heuristic.setInserter(edgeInserter);

		testSPRRType(heuristic, edgeInserter, RemoveReinsertType::None, "none");
		testSPRRType(heuristic, edgeInserter, RemoveReinsertType::Inserted, "inserted");
		testSPRRType(heuristic, edgeInserter, RemoveReinsertType::MostCrossed, "most-crossed");
		testSPRRType(heuristic, edgeInserter, RemoveReinsertType::All, "all");
		testSPRRType(heuristic, edgeInserter, RemoveReinsertType::Incremental, "incremental");
		testSPRRType(heuristic, edgeInserter, RemoveReinsertType::IncInserted, "inc-inserted");
	});
}

/**
 * Test variants of the ogdf::SubgraphPlanarizer.
 */
void testSubgraphPlanarizer() {
	describe("SubgraphPlanarizer", []() {
		testSPEdgeInserter(new FixedEmbeddingInserter, "FixedEmbedding");
		testSPEdgeInserter(new MultiEdgeApproxInserter, "MultiEdgeApprox");
		testSPEdgeInserter(new VariableEmbeddingInserter, "VariableEmbedding");
		testSPEdgeInserter(new VariableEmbeddingInserterDyn, "VariableEmbeddingDyn");
	});
}

/**
 * Test planarizers based on star insertion.
 */
void testStarInsertionPlanarizers() {
	describe("PlanarizerStarReinsertion", []() {
		for (int maxIter : {0, 1, 100, -1}) {
			PlanarizerStarReinsertion heuristic;
			heuristic.maxIterations(maxIter);
			testModule(heuristic, "maxIterations = " + std::to_string(maxIter), false);
		}
	});
	describe("PlanarizerStarReinsertion with PlanarizerMixedInsertion initialization", []() {
		for (int maxIter : {0, 1, 100, -1}) {
			PlanarizerStarReinsertion heuristic2;
			heuristic2.setPlanarization(new PlanarizerMixedInsertion);
			heuristic2.maxIterations(maxIter);
			testModule(heuristic2, "maxIterations = " + std::to_string(maxIter), false);
		}
	});
	describe("PlanarizerStarReinsertion with PlanarizerChordlessCycle initialization", []() {
		for (int maxIter : {0, 1, 100, -1}) {
			PlanarizerStarReinsertion heuristic2;
			heuristic2.setPlanarization(new PlanarizerChordlessCycle);
			heuristic2.maxIterations(maxIter);
			testModule(heuristic2, "maxIterations = " + std::to_string(maxIter), false);
		}
	});
	describe("PlanarizerMixedInsertion", []() {
		int i {0};
		for (PlanarizerMixedInsertion::NodeSelectionMethod method :
				{PlanarizerMixedInsertion::NodeSelectionMethod::Random,
						PlanarizerMixedInsertion::NodeSelectionMethod::HigherDegree,
						PlanarizerMixedInsertion::NodeSelectionMethod::LowerDegree,
						PlanarizerMixedInsertion::NodeSelectionMethod::HigherNonPlanarDegree,
						PlanarizerMixedInsertion::NodeSelectionMethod::LowerNonPlanarDegree,
						PlanarizerMixedInsertion::NodeSelectionMethod::BothEndpoints}) {
			PlanarizerMixedInsertion heuristic;
			heuristic.nodeSelectionMethod(method);
			testModule(heuristic, "nodeSelectionMethod = " + std::to_string(i), false);
			i++;
		}
	});
	describe("PlanarizerChordlessCycle", []() {
		PlanarizerChordlessCycle heuristic;
		testModule(heuristic, "default configuration", false);
	});
}

go_bandit([]() {
	testSubgraphPlanarizer();
	testStarInsertionPlanarizers();
});
