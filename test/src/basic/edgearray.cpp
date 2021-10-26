/** \file
 * \brief Tests for ogdf::EdgeArray
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
using namespace bandit;

go_bandit([]() {
	auto chooseEdge = [](const Graph& graph) { return graph.chooseEdge(); };

	auto allEdges = [](const Graph& graph, List<edge>& list) { graph.allEdges(list); };

	auto createEdge = [](Graph& graph) {
		return graph.newEdge(graph.chooseNode(), graph.chooseNode());
	};

	auto init = [](Graph& graph) { randomGraph(graph, 42, 168); };

	describeArray<Graph, EdgeArray, edge, int>("EdgeArray filled with ints", 42, 43, init,
			chooseEdge, allEdges, createEdge);
	describeArray<Graph, EdgeArray, edge, List<int>>("EdgeArray filled with lists of ints",
			{1, 2, 3}, {42}, init, chooseEdge, allEdges, createEdge);

	describeArray<Graph, EdgeArray, edge, bool>("EdgeArray filled with bools", false, true, init,
			chooseEdge, allEdges, createEdge);

	describeArrayWithoutDefault<Graph, EdgeArrayWithoutDefault, edge, std::unique_ptr<int>>(
			"EdgeArray filled with unique pointers", init, chooseEdge, allEdges, createEdge);
	describeArrayWithoutDefault<Graph, EdgeArrayWithoutDefault, edge, std::vector<std::unique_ptr<int>>>(
			"EdgeArray filled with vectors of unique pointers", init, chooseEdge, allEdges,
			createEdge);
});
