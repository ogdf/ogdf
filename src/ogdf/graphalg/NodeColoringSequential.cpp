/** \file
 * \brief Applies the sequential coloring algorithm to a given graph.
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

#include <ogdf/basic/Array.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/basic.h>
#include <ogdf/graphalg/NodeColoringModule.h>
#include <ogdf/graphalg/NodeColoringSequential.h>

#include <algorithm>
#include <vector>

namespace ogdf {

using NColor = NodeColoringModule::NodeColor;

void NodeColoringSequential::sortByDegree(List<node>& nodes) {
	// Store all nodes in an array
	Array<node> nodesArray(nodes.size());
	int i = 0;
	for (node v : nodes) {
		nodesArray[i++] = v;
	}

	// Sort the nodes increasing with the degree
	NodeDegreeComparer comparer;
	nodesArray.quicksort(comparer);

	// Determine a list with decreasing degree nodes
	nodes.clear();
	for (i = 0; i < nodesArray.size(); i++) {
		nodes.emplaceFront(nodesArray[i]);
	}
}

NColor NodeColoringSequential::colorByIndex(const Graph& graph, NodeArray<NColor>& colors,
		NColor start) {
	// Store all nodes in an array
	List<node> nodesList;
	for (node v : graph.nodes) {
		nodesList.emplaceBack(v);
	}

	// Perform the sequential coloring with the given permutation
	auto numberColors = fromPermutation(graph, colors, nodesList, start);

	// Check the resulting coloring
	OGDF_ASSERT(checkColoring(graph, colors));
	return numberColors;
}

NColor NodeColoringSequential::colorByDegree(const Graph& graph, NodeArray<NColor>& colors,
		NColor start) {
	// Store all nodes in a list
	List<node> nodesList;
	for (node v : graph.nodes) {
		nodesList.emplaceBack(v);
	}

	// Sort the node list decreasingly by the degree
	sortByDegree(nodesList);

	// Perform the sequential coloring with the given permutation
	auto numberColors = fromPermutation(graph, colors, nodesList, start);

	// Check the resulting coloring
	OGDF_ASSERT(checkColoring(graph, colors));
	return numberColors;
}

NColor NodeColoringSequential::fromPermutation(const Graph& graph, NodeArray<NColor>& colors,
		List<node>& nodePermutation, NColor start, bool lookForNeighbors) {
	// Copy the graph
	GraphCopy graphMain = GraphCopy(graph);
	preprocessGraph(graphMain);

	// Create a node table
	NodeArray<node> nodeTableOrig2New(graph);
	for (node vNew : graphMain.nodes) {
		node vOrig = graphMain.original(vNew);
		nodeTableOrig2New[vOrig] = vNew;
	}

	// Vector to store the neighbor colors
	std::vector<NColor> neighborColors;
	neighborColors.reserve(graphMain.numberOfNodes());

	// Vector to store which of the neighbor colors are already used
	std::vector<bool> isUsed;
	isUsed.reserve(graphMain.numberOfNodes() + 1);

	// Store the biggest color used
	NColor biggestColor = NColor(0);

	// Mark every node as uncolored
	NodeArray<bool> isColored(graphMain, false);
	if (lookForNeighbors) {
		isColored = NodeArray<bool>(graphMain, true);
		for (node vOrig : nodePermutation) {
			isColored[nodeTableOrig2New[vOrig]] = false;
		}
	}

	// Process the nodes in the given order
	for (node vOrig : nodePermutation) {
		OGDF_ASSERT(vOrig->graphOf() == &graph);
		node vNew = nodeTableOrig2New[vOrig];

		// Check if the node is already colored
		if (isColored[vNew]) {
			continue;
		}

		// Get the colors of the neighbors
		neighborColors.clear();
		for (adjEntry adj : vNew->adjEntries) {
			node w = adj->twinNode();
			if (isColored[w]) {
				NColor neighborColor = colors[graphMain.original(w)] - start;
				neighborColors.push_back(neighborColor);
			}
		}

		// Find the smallest unused color index
		isUsed.clear();
		isUsed.resize(neighborColors.size() + 1, false);
		for (NColor color = NColor(0); color < neighborColors.size(); color++) {
			if (neighborColors[color] <= neighborColors.size()) {
				isUsed[neighborColors[color]] = true;
			}
		}

		unsigned int smallestColorIndex;
		for (smallestColorIndex = NColor(0); smallestColorIndex < isUsed.size();
				smallestColorIndex++) {
			if (!isUsed[smallestColorIndex]) {
				break;
			}
		}

		// Assign an not conflicted color to the node
		colors[graphMain.original(vNew)] = smallestColorIndex + start;
		isColored[vNew] = true;
		biggestColor = std::max(biggestColor, smallestColorIndex + 1);
	}
	OGDF_ASSERT(checkColoring(graph, colors));
	return biggestColor;
}

}
