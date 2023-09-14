/** \file
 * \brief Tests for minimum cut algorithms
 *
 * \author Stephan Beyer
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

#include <ogdf/basic/graph_generators.h>
#include <ogdf/graphalg/MinimumCutNagamochiIbaraki.h>
#include <ogdf/graphalg/MinimumCutStoerWagner.h>

#include <testing.h>

template<typename T>
static void constructWeights(EdgeArray<T>& weights, std::initializer_list<T>&& consWeights) {
	typename std::initializer_list<T>::const_iterator cw {consWeights.begin()};
	for (auto& w : weights) {
		w = *cw;
		++cw;
	}
}

template<typename T>
static void minCutTest(std::string&& graphDesc, int numNodes,
		std::initializer_list<std::pair<int, int>>&& consGraph,
		std::initializer_list<T>&& consWeights, T expectedValue,
		std::initializer_list<int>&& expectedEdges, std::initializer_list<int>&& expectedNodes,
		MinimumCutModule<T>& minCut, bool providesPartition) {
	describe("on " + graphDesc, [&] {
		Graph graph;
		customGraph(graph, numNodes, consGraph);
		OGDF_ASSERT(int(consWeights.size()) == graph.numberOfEdges());
		OGDF_ASSERT(int(expectedEdges.size()) <= graph.numberOfEdges());
		OGDF_ASSERT(int(expectedNodes.size()) <= graph.numberOfNodes());

		EdgeArray<T> weights {graph};
		constructWeights(weights, std::move(consWeights));
		T value {minCut.call(graph, weights)};

		it("provides the correct value", [&] {
			AssertThat(value, Equals(expectedValue));
			AssertThat(minCut.value(), Equals(expectedValue));
		});

		if (providesPartition) {
			it("provides the correct cut edges", [&] {
				const ArrayBuffer<edge>& edges {minCut.edges()};
				AssertThat(edges.size(), Equals(int(expectedEdges.size())));

				ArrayBuffer<int> indices {edges.size()};
				for (edge e : edges) {
					indices.push(e->index());
				}
				indices.quicksort();

				int i {0};
				for (int e : expectedEdges) {
					AssertThat(indices[i], Equals(e));
					++i;
				}
			});

			it("provides the correct cut nodes", [&] {
				const ArrayBuffer<node>& nodes {minCut.nodes()};
				AssertThat(nodes.size(),
						Equals(int(expectedNodes.size()))
								|| Equals(graph.numberOfNodes() - int(expectedNodes.size())));

				ArrayBuffer<int> indices1(nodes.size());
				ArrayBuffer<int> indices2(graph.numberOfNodes() - nodes.size());
				for (node v : nodes) {
					indices1.push(v->index());
				}
				indices1.quicksort();
				for (int index {}, j {}; index < graph.numberOfNodes(); ++index) {
					if (j >= indices1.size() || index == indices1[j]) {
						++j;
					} else {
						indices2.push(index);
					}
				}

				int i {0};
				for (int v : expectedNodes) {
					if (i >= indices1.size()) {
						AssertThat(v, Equals(indices2[i]));
					} else if (i >= indices2.size()) {
						AssertThat(v, Equals(indices1[i]));
					} else {
						AssertThat(v, Equals(indices1[i]) || Equals(indices2[i]));
					}
					++i;
				}
			});
		}
	});
}

template<typename T>
static void minCutTests(std::string&& type, MinimumCutModule<T>& minCut,
		bool providesPartition = false) {
	describe("with weights of type " + type, [&] {
		minCutTest<T>("an empty graph", 0, {}, {}, std::numeric_limits<T>::max(), {}, {}, minCut,
				providesPartition);
		minCutTest<T>("one isolated node", 1, {}, {}, std::numeric_limits<T>::max(), {}, {}, minCut,
				providesPartition);
		minCutTest<T>("two isolated nodes", 2, {}, {}, 0, {}, {0}, minCut, providesPartition);
		minCutTest<T>("two isolated nodes with self-loops", 2, {{0, 0}, {1, 1}}, {3, 3}, 0, {}, {0},
				minCut, providesPartition);
		minCutTest<T>("a single-edge graph with zero weight", 2, {{0, 1}}, {0}, 0, {0}, {0}, minCut,
				providesPartition);
		minCutTest<T>("a single-edge graph with non-zero weight", 2, {{0, 1}}, {3}, 3, {0}, {0},
				minCut, providesPartition);
		minCutTest<T>("a graph with parallel edges", 2, {{0, 1}, {0, 1}, {0, 1}}, {1, 1, 1}, 3,
				{0, 1, 2}, {0}, minCut, providesPartition);
		minCutTest<T>("a 4-cycle", 4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}}, {5, 2, 4, 1}, 3, {1, 3},
				{0, 1}, minCut, providesPartition);
		minCutTest<T>("a graph with unit weights", 8,
				{{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}, {4, 0}, {4, 5}, {5, 6}, {6, 4},
						{4, 7}, {5, 7}, {6, 7}, {5, 2}},
				{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 2, {6, 13}, {0, 1, 2, 3}, minCut,
				providesPartition);

		it("on a big graph", [&] {
			Graph graph;
			completeGraph(graph, 500);

			EdgeArray<T> weights {graph};
			for (T& w : weights) {
				w = randomNumber(1, 100);
			}

			minCut.call(graph, weights);
			AssertThat(minCut.value(), IsGreaterThan(0));
			if (providesPartition) {
				AssertThat(minCut.nodes(), !IsEmpty());
			}
		});
	});
}

go_bandit([] {
	describe("StoerWagner", [] {
		MinimumCutStoerWagner<double> minCut;
		minCutTests<double>("double", minCut, true);

		MinimumCutStoerWagner<int> minCut2;
		minCutTests<int>("int", minCut2, true);
	});

	describe("NagamochiIbaraki", [] {
		describe("no preprocessing and no Padberg-Rinaldi heuristics", [] {
			MinimumCutNagamochiIbaraki minCut {false, false, Logger::Level::Force};
			minCutTests<int>("int", minCut);
		});
		describe("no preprocessing but with Padberg-Rinaldi heuristics", [] {
			MinimumCutNagamochiIbaraki minCut {true, false, Logger::Level::Force};
			minCutTests<int>("int", minCut);
		});
		describe("with preprocessing, but no Padberg-Rinaldi heuristics", [] {
			MinimumCutNagamochiIbaraki minCut {false, true, Logger::Level::Force};
			minCutTests<int>("int", minCut);
		});
		describe("with preprocessing and Padberg-Rinaldi heuristics", [] {
			MinimumCutNagamochiIbaraki minCut {true, true, Logger::Level::Force};
			minCutTests<int>("int", minCut);
		});
	});
});
