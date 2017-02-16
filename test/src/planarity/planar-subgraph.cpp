#include <random>

#include <bandit/bandit.h>
#include <resources.h>

#include <ogdf/planarity/PlanarSubgraphBoyerMyrvold.h>
#include <ogdf/planarity/PlanarSubgraphFast.h>
#include <ogdf/planarity/MaximalPlanarSubgraphSimple.h>
#include <ogdf/planarity/MaximumPlanarSubgraph.h>
#include <ogdf/planarity/PlanarSubgraphCactus.h>
#include <ogdf/planarity/PlanarSubgraphTree.h>
#include <ogdf/basic/graph_generators.h>

using namespace ogdf;
using namespace bandit;
using std::minstd_rand;

template <typename TCost>
void testSubgraphInstance(Graph &graph, PlanarSubgraphModule<TCost> &psm, PlanarityModule &tester, bool assertMaximality, bool weighEdges, bool connects) {
	makeSimpleUndirected(graph);
	makeConnected(graph);

	EdgeArray<TCost> costs;
	edge mustHaveEdge = nullptr;

	if(weighEdges) {
		costs.init(graph);

		for(edge e : graph.edges) {
			costs[e] = TCost(1);
		}

		mustHaveEdge = graph.chooseEdge();
		costs[mustHaveEdge] = TCost(1 + graph.numberOfEdges());
	}

	List<edge> removedEdges;
	if(weighEdges) {
		psm.call(graph, costs, removedEdges);
	} else {
		psm.call(graph, removedEdges);
	}

	cout << std::endl << "      removed " << removedEdges.size() << " edges" << std::endl;
	bool connected = isConnected(graph);

	Graph::HiddenEdgeSet set(graph);
	List<edge> hiddenEdges;
	for(edge e : removedEdges) {
		set.hide(e);
		hiddenEdges.pushBack(e);
		AssertThat(e, Is().Not().EqualTo(mustHaveEdge));
	}

	if(connects) {
		AssertThat(isConnected(graph), Equals(connected));
	}

	AssertThat(tester.isPlanar(graph), Equals(true));

	if(assertMaximality) {
		ListConstIterator<edge> it = removedEdges.begin();
		for(edge e : hiddenEdges) {
			set.restore(e);
			AssertThat(tester.isPlanar(graph), Equals(false));
			set.hide(e);
			it = it.succ();
		}
	}
}

void testSubgraphInstanceForIntAndDouble(Graph &graph, PlanarSubgraphModule<int> &psmi, PlanarSubgraphModule<double> &psmd) {
	makeSimpleUndirected(graph);
	makeConnected(graph);
	EdgeArray<int> costsInt(graph);
	EdgeArray<double> costsDouble(graph);
	int cnt = 0;
	for(edge e : graph.edges) {
		cnt++;
		costsInt(e) = cnt;
		costsDouble(e) = double(cnt);
	}
	List<edge> removedEdgesInt;
	List<edge> removedEdgesDouble;
	psmi.call(graph, costsInt, removedEdgesInt);
	psmd.call(graph, costsDouble, removedEdgesDouble);
	int sumOfRemovedInt = 0;
	double sumRemovedDouble = 0;
	for (edge e: removedEdgesInt) { sumOfRemovedInt += costsInt(e); }
	for (edge e: removedEdgesDouble) { sumRemovedDouble += costsDouble(e); }
	AssertThat( sumOfRemovedInt, Equals(int(sumRemovedDouble)) );
}

void performGenericTests(const string &name, bool optimal, bool respectsEdgeWeight, bool skip, std::function<void(Graph &, bool)> callFunc) {
	describe(name, [&]() {
		BoothLueker bl;
		minstd_rand rng(42);

		for(int n = 4; n <= (optimal ? 10 : 30); n+=2) {
			int m = n*(n-1)/3;
			it(string("works on a dense random graph with " + to_string(n) + " nodes and " + to_string(m) + " edges"), [&]() {
				Graph graph;
				randomGraph(graph, n, m);
				callFunc(graph, false);
			});
		}

		for(int n = 30; n <= (optimal ? 0 : 50); n+=5) {
			int m = 4*n;
			it(string("works on a sparse random graph with " + to_string(n) + " nodes and " + to_string(m) + " edges"), [&]() {
				Graph graph;
				randomGraph(graph, n, m);
				callFunc(graph, false);
			});
		}

		auto testOnComplete = [&callFunc](std::initializer_list<int> vertices) {
			for (int n : vertices) {
				it(string("works on a K" + to_string(n) + " with weighted edges"), [&]() {
					Graph graph;
					completeGraph(graph, n);
					callFunc(graph, true);
				});
			}
		};

		if(respectsEdgeWeight) {
			if (optimal) {
				testOnComplete({5, 6});
			} else {
				testOnComplete({5, 7, 9, 11, 13, 15, 17, 19});
			}

			for (int n = 10; n <= (optimal ? 12 : 50); n+= (optimal ? 1 : 5)) {
				int m = 3 * n;
				it(string(
					"works on a sparse weighted random graph with " + to_string(n) + " nodes and " + to_string(m) +
						" edges"), [&]() {
					Graph graph;
					randomGraph(graph, n, m);
					callFunc(graph, true);
				});
			}
		}

		std::vector<string> instances = {
			"north/g.61.11.gml",
			"rome/grafo3703.45.lgr.gml.pun",
			"rome/grafo5745.50.lgr.gml.pun"
		};

		for(int i = 0; i < (optimal ? 0 : 2); i++) {
			string tags = "";
			bool weighted = i;

			if(weighted) {
				tags += " weighted";
			}

			if(!weighted || respectsEdgeWeight) {
				for_each_graph_it(("works on" + tags), instances,
				                  [&](Graph &graph, const string &filename) {
					callFunc(graph, weighted);
				});
			}
		}
	}, skip);
}

template <typename TCost>
void testSubgraphAlgorithm(const string &name, PlanarSubgraphModule<TCost> &psm, bool optimal, bool maximal, bool respectsEdgeWeight, bool connects, bool skip = false) {
	maximal |= optimal;
	BoothLueker bl;
	performGenericTests(name, optimal, respectsEdgeWeight, skip, [&](Graph &graph, bool weighEdges) {
		testSubgraphInstance<TCost>(graph, psm, bl, maximal, weighEdges, connects);
	});
}

void testSubgraphAlgorithmForIntAndDouble(const string &name, PlanarSubgraphModule<int> &psmi, PlanarSubgraphModule<double> &psmd) {
	performGenericTests(name + " int VS double", false, true, false, [&](Graph &graph, bool weighEdges) {
		testSubgraphInstanceForIntAndDouble(graph, psmi, psmd);
	});
}

template<template<typename> class Algorithm>
void describeAlgorithm(const string &name, bool optimal, bool maximal, bool respectsEdgeWeight, bool connects, bool skip = false) {
	Algorithm<int> algoInt;
	testSubgraphAlgorithm(name + "<int>", algoInt, optimal, maximal, respectsEdgeWeight, connects, skip);

	Algorithm<double> algoDouble;
	testSubgraphAlgorithm(name + "<double>", algoDouble, optimal, maximal, respectsEdgeWeight, connects, skip);
}

go_bandit([]() {
	describe("Planar Subgraphs", []() {
		PlanarSubgraphBoyerMyrvold bms;
		MaximumPlanarSubgraph mps;

		testSubgraphAlgorithm("PlanarSubgraphBoyerMyrvold",              bms, false, false, true,  true, true);
		describeAlgorithm<PlanarSubgraphFast>("PlanarSubgraphFast",           false, false, false, true);
		describeAlgorithm<PlanarSubgraphCactus>("PlanarSubgraphCactus",       false, false, false, true);
		describeAlgorithm<PlanarSubgraphTree>("PlanarSubgraphTree",           false, false, false, true);
		testSubgraphAlgorithm("MaximumPlanarSubgraph",                   mps, true,  true,  true,  true);
		describeAlgorithm<PlanarSubgraphEmpty>("PlanarSubgraphEmpty",         false, false, false, false);

		MaximalPlanarSubgraphSimple<int> mpss;
		PlanarSubgraphCactus<int> psc;
		MaximalPlanarSubgraphSimple<int> mpssPsc(psc);
		PlanarSubgraphFast<int> fps;
		MaximalPlanarSubgraphSimple<int> mpssFps(fps);
		MaximalPlanarSubgraphSimple<int> mpssBms(bms);

		testSubgraphAlgorithm("MaximalPlanarSubgraphSimple",        mpss,    false, true, false, true);
		testSubgraphAlgorithm("Maximal PlanarSubgraphCactus",       mpssPsc, false, true, false, true);
		testSubgraphAlgorithm("Maximal PlanarSubgraphFast",         mpssFps, false, true, false, true);
		testSubgraphAlgorithm("Maximal PlanarSubgraphBoyerMyrvold", mpssBms, false, true, false, true, true);

		PlanarSubgraphCactus<double> pscd;
		MaximalPlanarSubgraphSimple<double> mpssPscd(pscd);
		testSubgraphAlgorithmForIntAndDouble("Maximal PlanarSubgraphCactus", mpssPsc, mpssPscd);
	});
});
