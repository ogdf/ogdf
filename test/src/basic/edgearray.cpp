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
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphSets.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/RegisteredSet.h>
#include <ogdf/basic/graph_generators/randomized.h>

#include <functional>
#include <string>

#include "array_helper.h"
#include <testing.h>

go_bandit([]() {
	auto chooseEdge = [](const Graph& graph) { return graph.chooseEdge(); };

	auto allEdges = [](const Graph& graph, List<edge>& list) { graph.allEdges(list); };

	auto createEdge = [](Graph& graph) {
		return graph.newEdge(graph.chooseNode(), graph.chooseNode());
	};

	auto deleteEdge = [](Graph& graph, edge e) { return graph.delEdge(e); };

	auto clearEdges = [](Graph& graph) { return graph.clear(); };

	auto init = [](Graph& graph) { randomGraph(graph, 42, 168); };

	runBasicArrayTests<Graph, EdgeArray, edge>("EdgeArray", init, chooseEdge, allEdges, createEdge);
	runBasicSetTests<Graph, EdgeSet, edge>("EdgeSet", init, chooseEdge, allEdges, createEdge,
			deleteEdge, clearEdges);
});
