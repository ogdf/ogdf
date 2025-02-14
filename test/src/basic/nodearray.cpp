/** \file
 * \brief Tests for ogdf::NodeArray
 *
 * \author Mirko Wagner, Tilo Wiedera
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
#include <ogdf/basic/GraphSets.h>
#include <ogdf/basic/RegisteredSet.h>
#include <ogdf/basic/graph_generators.h>

#include "array_helper.h" // IWYU pragma: associated
#include <testing.h>

go_bandit([]() {
	auto chooseNode = [](const Graph& graph) { return graph.chooseNode(); };

	auto allNodes = [](const Graph& graph, List<node>& list) { graph.allNodes(list); };

	auto createNode = [](Graph& graph) { return graph.newNode(); };

	auto deleteNode = [](Graph& graph, node n) { return graph.delNode(n); };

	auto clearNodes = [](Graph& graph) { graph.clear(); };

	auto init = [](Graph& graph) { randomGraph(graph, 42, 168); };

	runBasicArrayTests<Graph, NodeArray, node>("NodeArray", init, chooseNode, allNodes, createNode);

	describe("NodeArray filled with pointers", [&]() {
		Graph G;
		init(G);

		it("initializes with nullptr values", [&]() {
			NodeArray<int*> arr1(G);
			NodeArray<int*, false> arr2(G);
			AssertThat(arr1[chooseNode(G)], IsNull());
			AssertThat(arr2[chooseNode(G)], IsNull());
		});

		it("initializes with a default value", [&]() {
			std::unique_ptr<int> p(new int(42));
			NodeArray<int*> arr(G, p.get());
			AssertThat(arr[chooseNode(G)], Equals(p.get()));
		});
	});

	runBasicSetTests<Graph, NodeSet, node>("NodeSet", init, chooseNode, allNodes, createNode,
			deleteNode, clearNodes);
});
