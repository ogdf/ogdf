/** \file
 * \brief Applies the node coloring heuristic "Recursive Largest First"
 * specified by Leighton (1979).
 *
 * \author Jan-Niklas Buckow
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
#include <ogdf/basic/List.h>
#include <ogdf/basic/basic.h>
#include <ogdf/graphalg/NodeColoringModule.h>
#include <ogdf/graphalg/NodeColoringRecursiveLargestFirst.h>

namespace ogdf {

using NColor = NodeColoringModule::NodeColor;

NColor NodeColoringRecursiveLargestFirst::call(const Graph& graph, NodeArray<NColor>& colors,
		NColor start) {
	auto numberColorsUsed = NColor(0);
	// Copy the input graph
	GraphCopy graphMain = GraphCopy(graph);
	preprocessGraph(graphMain);

	// Perform the coloring while the graph is not empty
	while (!graphMain.empty()) {
		// Search for big independent sets
		GraphCopy gSubgraph = graphMain;

		// Array to store the degrees into the not available nodes
		NodeArray<int> degreesUnavailable(gSubgraph, 0);

		// Color the next independent set
		while (!gSubgraph.empty()) {
			// Search a fitting max degree node either in the set of available nodes or in the set of unavailable
			// nodes
			node maxDegreeNode;
			getCandidate(maxDegreeNode, degreesUnavailable, gSubgraph);

			// Color the node
			colors[gSubgraph.original(maxDegreeNode)] = start + numberColorsUsed;

			// Delete the node and its neighbors from the subgraph
			List<node> nodes_to_delete;
			for (adjEntry adj : maxDegreeNode->adjEntries) {
				nodes_to_delete.emplaceBack(adj->twinNode());
			}
			for (node v : nodes_to_delete) {
				for (adjEntry adj : v->adjEntries) {
					degreesUnavailable[adj->twinNode()]++;
				}
				gSubgraph.delNode(v);
			}

			node maxDegOrigNode = gSubgraph.original(maxDegreeNode);
			gSubgraph.delNode(maxDegreeNode);

			// Delete the already colored node from the graph
			graphMain.delNode(graphMain.copy(maxDegOrigNode));
		}
		// Increment the color the next independent set
		numberColorsUsed++;
	}
	// Check the coloring
	OGDF_ASSERT(checkColoring(graph, colors));
	return numberColorsUsed;
}

void NodeColoringRecursiveLargestFirst::getCandidate(node& candidate,
		NodeArray<int>& degreesUnavailable, Graph& graph) {
	OGDF_ASSERT(degreesUnavailable.graphOf() == &graph);
	candidate = graph.firstNode();
	int maxDegreeUnavailable = 0;
	int maxDegree = 0;
	for (node v : graph.nodes) {
		if (degreesUnavailable[v] >= maxDegreeUnavailable) {
			if (degreesUnavailable[v] > maxDegreeUnavailable || v->degree() > maxDegree) {
				candidate = v;
				maxDegreeUnavailable = degreesUnavailable[v];
				maxDegree = v->degree();
			}
		}
	}
}

}
