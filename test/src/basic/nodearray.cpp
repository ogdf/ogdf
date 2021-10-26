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
#include "array_helper.h"

using namespace ogdf;
using namespace snowhouse;
using namespace bandit;

go_bandit([]() {
	auto chooseNode = [](const Graph& graph) { return graph.chooseNode(); };

	auto allNodes = [](const Graph& graph, List<node>& list) { graph.allNodes(list); };

	auto createNode = [](Graph& graph) { return graph.newNode(); };

	auto init = [](Graph& graph) { randomGraph(graph, 42, 168); };

	describeArray<Graph, NodeArray, node, int>("NodeArray filled with ints", 42, 43, init,
			chooseNode, allNodes, createNode);
	describeArray<Graph, NodeArray, node, List<int>>("NodeArray filled with lists of ints",
			{1, 2, 3}, {42}, init, chooseNode, allNodes, createNode);

	describeArray<Graph, NodeArray, node, bool>("NodeArray filled with bools", false, true, init,
			chooseNode, allNodes, createNode);

	describeArrayWithoutDefault<Graph, NodeArrayWithoutDefault, node, std::unique_ptr<int>>(
			"NodeArray filled with unique pointers", init, chooseNode, allNodes, createNode);
	describeArrayWithoutDefault<Graph, NodeArrayWithoutDefault, node, std::vector<std::unique_ptr<int>>>(
			"NodeArray filled with vectors of unique pointers", init, chooseNode, allNodes,
			createNode);

	describe("NodeArray filled with pointers", [&]() {
		Graph G;
		init(G);

		it("initializes with nullptr values", [&]() {
			NodeArray<int*> arr1(G);
			NodeArrayWithoutDefault<int*> arr2(G);
			AssertThat(arr1[chooseNode(G)], IsNull());
			AssertThat(arr2[chooseNode(G)], IsNull());
		});

		it("initializes with a default value", [&]() {
			std::unique_ptr<int> p(new int(42));
			NodeArray<int*> arr(G, p.get());
			AssertThat(arr[chooseNode(G)], Equals(p.get()));
		});
	});
});
