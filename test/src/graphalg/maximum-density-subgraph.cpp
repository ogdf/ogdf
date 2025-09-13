/** \file
 * \brief Tests for maximum density subgraph algorithm.
 *
 * \author Finn Stutzenstein
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
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/GraphSets.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/graph_generators/deterministic.h>
#include <ogdf/graphalg/MaximumDensitySubgraph.h>

#include <functional>
#include <set>
#include <string>
#include <utility>

#include <graphs.h>

#include <testing.h>

void test(const std::string& description, std::function<void(Graph& G, NodeSet& expected)> setup) {
	Graph G;
	NodeSet expected(G);
	setup(G, expected);

	it(description, [&] {
		GraphCopySimple GC(G);
		NodeSet subgraphNodes(G);
		bool success =
				maximumDensitySubgraph(GC, subgraphNodes, [&](node n) { return GC.original(n); });
		AssertThat(success, IsTrue());

		AssertThat(subgraphNodes.size(), Equals(expected.size()));
		for (node n : expected.nodes()) {
			AssertThat(subgraphNodes.isMember(n), IsTrue());
		}
	});
}

void testContainsAll(const std::string& description, std::function<void(Graph& G)> setup) {
	test(description, [&](Graph& G, NodeSet& expected) {
		setup(G);
		for (node n : G.nodes) {
			expected.insert(n);
		}
	});
}

go_bandit([] {
	describe("MaximumDensitySubgraph", [] {
		describe("No errors on general graphs", [] {
			forEachGraphItWorks(
					{GraphProperty::simple},
					[&](const Graph& G) {
						GraphCopySimple GC(G);
						NodeSet subgraphNodes(G);
						bool success = maximumDensitySubgraph(GC, subgraphNodes,
								[&](node n) { return GC.original(n); });
						AssertThat(success, IsTrue());

						if (G.numberOfEdges() > 0) {
							AssertThat(subgraphNodes.size(), Is().GreaterThan(1));
						}
					},
					GraphSizes(10, 50, 20), 0, 50);
		});

		describe("Small/special cases", [] {
			test("works for a single node", [](Graph& G, NodeSet& expected) { G.newNode(); });
			test("works for nodes without edges",
					[](Graph& G, NodeSet& expected) { emptyGraph(G, 3); });
			test("works for a single edge", [](Graph& G, NodeSet& expected) {
				customGraph(G, 2, {{0, 1}});
				expected.insert(G.firstNode());
				expected.insert(G.lastNode());
			});
			testContainsAll("contains all nodes of a K5", [](Graph& G) { completeGraph(G, 5); });
			testContainsAll("contains all nodes of a tree", [](Graph& G) { regularTree(G, 50, 4); });
			testContainsAll("contains all nodes of a wheel graph",
					[](Graph& G) { wheelGraph(G, 10); });
			testContainsAll("contains all nodes of a cube graph", [](Graph& G) { cubeGraph(G, 6); });
			test("finds the K5", [](Graph& G, NodeSet& expected) {
				completeGraph(G, 5);
				List<node> K5nodes;
				G.allNodes(K5nodes);

				for (node v : K5nodes) {
					expected.insert(v);
					G.newEdge(G.newNode(), v);
				}
			});
			it("timeouts", [] {
				Graph G;
				completeGraph(G, 5);
				NodeSet subgraphNodes(G);
				bool success = maximumDensitySubgraph(
						G, subgraphNodes, [&](node n) { return n; }, 0);
				AssertThat(success, IsFalse());
			});
		});
	});
});
