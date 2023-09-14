/** \file
 * \brief Tests for planar separator algorithms.
 *
 * \author Thomas Klein
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

#include <ogdf/basic/graph_generators/deterministic.h>
#include <ogdf/graphalg/SeparatorDual.h>
#include <ogdf/graphalg/SeparatorDualFC.h>
#include <ogdf/graphalg/SeparatorHarPeled.h>
#include <ogdf/graphalg/SeparatorLiptonTarjan.h>
#include <ogdf/graphalg/SeparatorLiptonTarjanFC.h>
#include <ogdf/graphalg/ShortestPathAlgorithms.h>

#include <iostream>

#include <graphs.h>
#include <testing.h>

/**
 * Creates a random planar graph of a given size.
 *
 * @param n number of nodes of the graph
 * @return the graph
 */
Graph getRandomPlanarGraph(int n) {
	Graph graph;
	int edges = randomNumber(n, 3 * n - 6);
	randomPlanarConnectedGraph(graph, n, edges);
	return graph;
}

/**
 * Checks if the size of a separator was smaller than 2 * diameter + 1, the size
 * guarantee of FC-algorithms.
 * Calculating the actual diameter would be prohibitively expensive, so instead,
 * this method walks over all the shortest paths and checks if that distance was
 * large enough, returning at the first shortest path that was long enough. This
 * way, the test does not take very long, usually.
 *
 * @param sepSize the size (=number of nodes) of a separation
 * @param G the graph
 * @return true if the separator was small enough
 */
bool checkSizeAgainstDiameter(int sepSize, const Graph& G) {
	NodeArray<int> distance;
	distance.init(G);
	int edgeCosts = 1;

	for (node v : G.nodes) {
		bfs_SPSS(v, G, distance, edgeCosts);
		for (node w : G.nodes) {
			if (2 * distance[w] + 1 >= sepSize) {
				return true;
			}
		}
	}
	return false;
}

/**
 * Tests whether a node appears twice in any of the lists.
 *
 * @param G the graph
 * @param sep list of nodes in the separator
 * @param first list of nodes in the first half
 * @param second list of nodes in the second half
 * @return true if every node appeared exactly once
 */
bool testListCompleteness(const Graph& G, const List<node>& sep, const List<node>& first,
		const List<node>& second) {
	NodeArray<bool> marked;
	marked.init(G, false);

	auto listCheck = [&marked](const List<node>& list) -> bool {
		for (node no : list) {
			if (marked[no]) {
				return false;
			}
			marked[no] = true;
		}
		return true;
	};

	// making sure that no node was mentioned twice
	if (!listCheck(sep) || !listCheck(first) || !listCheck(second)) {
		return false;
	}

	// making sure that no node was forgotten
	for (node no : G.nodes) {
		if (!marked[no]) {
			return false;
		}
	}

	return true;
}

/**
 * Tests whether the separator actually separates the graph.
 *
 * @param G the graph
 * @param sep the separator
 * @param first the first half
 * @param second the second half
 * @return true if the graph is separated
 */
bool testSeparatorCorrectness(const Graph& G, const List<node>& sep, const List<node>& first,
		const List<node>& second) {
	if (G.empty()) {
		return true;
	}
	GraphCopy graphCopy(G);

	for (node no : sep) {
		graphCopy.delNode(graphCopy.copy(no));
	}

	NodeArray<bool> marked;
	marked.init(graphCopy, false);

	ArrayBuffer<node> buffer;

	// start BFS at each node of first
	for (node v : first) {
		node v_copy = graphCopy.copy(v);
		if (marked[v_copy]) {
			continue;
		}

		buffer.push(v_copy);
		marked[v_copy] = true;

		while (!buffer.empty()) {
			node w = buffer.popRet();
			for (adjEntry adj : w->adjEntries) {
				node x = adj->twinNode();
				if (!marked[x]) {
					marked[x] = true;
					buffer.push(x);
				}
			}
		}
	}

	// now make sure that none of the nodes in second were visited
	for (node no : second) {
		if (marked[graphCopy.copy(no)]) {
			return false;
		}
	}

	return true;
}

/**
 * Tests whether the size constraints on the separation are fulfilled, ie
 * 	1. the separator is not larger than the value guaranteed by the algorithm
 * 	2. no list contains more than 2/3 of the nodes
 * 	3. sep.size + first.size + second.size = n
 *
 * @param G the graph
 * @param sep the separator
 * @param first the first half
 * @param second the second half
 * @param maxSize the maximal size of the separator as guaranteed by the algorithm
 */
void assertListSizes(const Graph& G, const List<node>& sep, const List<node>& first,
		const List<node>& second, double maxSize) {
	int n = G.numberOfNodes();
	AssertThat(n == sep.size() + first.size() + second.size(), IsTrue());
	AssertThat(first.size() <= 2.0 / 3.0 * n, IsTrue());
	AssertThat(second.size() <= 2.0 / 3.0 * n, IsTrue());

	if (!G.empty()) {
		if (maxSize > 0) {
			AssertThat(sep.size() < maxSize, IsTrue());
		} else {
			AssertThat(checkSizeAgainstDiameter(sep.size(), G), IsTrue());
		}
	}
}

/**
 * Runs all standard tests on a separation.
 *
 * @param G the graph
 * @param sep the separator
 * @param first the first half
 * @param second the second half
 * @param maxSize the maximal size of the separator as guaranteed by the algorithm
 */
void testSeparator(const Graph& G, const List<node>& sep, const List<node>& first,
		const List<node>& second, double maxSize) {
	AssertThat(testListCompleteness(G, sep, first, second), IsTrue());
	AssertThat(testSeparatorCorrectness(G, sep, first, second), IsTrue());
	assertListSizes(G, sep, first, second, maxSize);
}

/**
 * Tests a given separator algorithm on a given graph.
 *
 * @param G the graph to be tested
 * @param sep the planar separator algorithm
 */
void testGraph(Graph& G, PlanarSeparatorModule& sep) {
	List<node> separator;
	List<node> first;
	List<node> second;

	makeSimpleUndirected(G);
	planarEmbedPlanarGraph(G);

	sep.separate(G, separator, first, second);
	double maxSize = sep.getMaxSeparatorSize(G.numberOfNodes());
	testSeparator(G, separator, first, second, maxSize);
}

/**
 * Applies a planar separator algorithm to simple planar graphs.
 *
 * @param sep the separator algorithm to be tested
 */
static void describeSimplePlanarGraphs(PlanarSeparatorModule& sep) {
	forEachGraphItWorks(
			std::set<GraphProperty>({GraphProperty::planar, GraphProperty::simple}),
			[&sep](Graph& testG) { testGraph(testG, sep); }, GraphSizes(), 0);
}

/**
 * Applies a planar separator algorithm to a few random instances.
 *
 * @param sep the separator algorithm to be tested
 */
static void describeRandomInstances(PlanarSeparatorModule& sep) {
	int numInstances = 5;
	int size = 500;

	for (int i = 0; i < numInstances; i++) {
		it("works on a random planar graph number " + std::to_string(i), [&] {
			setSeed(i);
			Graph G = getRandomPlanarGraph(size);
			testGraph(G, sep);
		});
	}
}

/**
 * Applies a planar separator algorithm to a globe instance.
 * This is mainly interesting for Har-Peled, which goes to exit point Region_R
 * for these instances.
 *
 * @param sep the separator algorithm to be tested
 */
static void describeGlobeInstances(PlanarSeparatorModule& sep) {
	it("works on a globe graph", [&] {
		Graph G;
		globeGraph(G, 50, 50);
		testGraph(G, sep);
	});
}

go_bandit([] {
	describe("PlanarSeparators", [] {
		// collection of all planar separator algorithms
		std::vector<PlanarSeparatorModule*> separators;

		SeparatorLiptonTarjan sepLipTar1;
		SeparatorLiptonTarjan sepLipTar2(true, 2);
		SeparatorDual sepDual1;
		SeparatorDual sepDual2(true, 2);
		SeparatorLiptonTarjanFC sepLTFC(true);
		SeparatorDualFC sepDFC(true);
		SeparatorHarPeled sepHarPel;

		separators.push_back(&sepLipTar1);
		separators.push_back(&sepLipTar2);
		separators.push_back(&sepDual1);
		separators.push_back(&sepDual2);
		separators.push_back(&sepLTFC);
		separators.push_back(&sepDFC);
		separators.push_back(&sepHarPel);

		// for all separators
		for (auto sep : separators) {
			setSeed(42);

			describe(sep->getName(), [&sep] {
				// test basic planar instances
				describeSimplePlanarGraphs(*sep);

				// test on a few randomly generated instances
				describeRandomInstances(*sep);

				// test on globe instances
				describeGlobeInstances(*sep);
			});
		}
	});
});
