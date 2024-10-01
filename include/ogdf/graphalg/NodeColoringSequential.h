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

#pragma once

#include <ogdf/basic/GraphCopy.h>
#include <ogdf/graphalg/NodeColoringModule.h>

namespace ogdf {

/**
 * Approximation algorithms for the node coloring problem in graphs.
 * This class applies the sequential coloring algorithm which greedily
 * assigns to each node the first available color.
 */
class NodeColoringSequential : public NodeColoringModule {
public:
	virtual NodeColor call(const Graph& graph, NodeArray<NodeColor>& colors,
			NodeColor start = 0) override {
		return colorByDegree(graph, colors, start);
	}

	/**
	 * Sorts a list of nodes decreasingly by the degree.
	 * @param nodes List of nodes to be sorted.
	 */
	virtual void sortByDegree(List<node>& nodes) {
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
		for (int i = 0; i < nodesArray.size(); i++) {
			nodes.emplaceFront(nodesArray[i]);
		}
	}

	/**
	 * Performs the sequential nodes coloring.
	 * The nodes are thereby sorted by the largest degree.
	 *
	 * @param graph The graph to be colored
	 * @param colors The resulting coloring
	 * @param start The starting color index
	 * @return The number of colors used
	 */
	virtual NodeColor colorByIndex(const Graph& graph, NodeArray<NodeColor>& colors,
			NodeColor start = 0) {
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

	/**
	 * Performs the sequential nodes coloring.
	 * The nodes are thereby sorted by the smallest index.
	 *
	 * @param graph The graph to be colored
	 * @param colors The resulting coloring
	 * @param start The starting color index
	 * @return The number of colors used
	 */
	virtual NodeColor colorByDegree(const Graph& graph, NodeArray<NodeColor>& colors,
			NodeColor start = 0) {
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

	/**
	 * Colors a graph with the sequential coloring algorithm from a given node permutation.
	 * @param graph The graph to be colored
	 * @param colors The resulting coloring
	 * @param nodePermutation The given node permutation
	 * @param start The starting color index
	 * @param lookForNeighbors Iff true, all nodes, including those not contained in the list will be considered
	 * @return The number of colors used
	 */
	virtual NodeColor fromPermutation(const Graph& graph, NodeArray<NodeColor>& colors,
			List<node>& nodePermutation, NodeColor start = 0, bool lookForNeighbors = false) {
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
		std::vector<NodeColor> neighborColors;
		neighborColors.reserve(graphMain.numberOfNodes());

		// Vector to store which of the neighbor colors are already used
		std::vector<bool> isUsed;
		isUsed.reserve(graphMain.numberOfNodes() + 1);

		// Store the biggest color used
		NodeColor biggestColor = NodeColor(0);

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
					NodeColor neighborColor = colors[graphMain.original(w)] - start;
					neighborColors.push_back(neighborColor);
				}
			}

			// Find the smallest unused color index
			isUsed.clear();
			isUsed.resize(neighborColors.size() + 1, false);
			for (NodeColor color = NodeColor(0); color < neighborColors.size(); color++) {
				if (neighborColors[color] <= neighborColors.size()) {
					isUsed[neighborColors[color]] = true;
				}
			}

			unsigned int smallestColorIndex;
			for (smallestColorIndex = NodeColor(0); smallestColorIndex < isUsed.size();
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

private:
	/**
	 * Class for comparing two nodes by the node degree.
	 */
	class NodeDegreeComparer {
	public:
		/**
		 * Compares two nodes by using the node degree.
		 * A higher nodes degree indicates a higher value.
		 *
		 * @param v1 First nodes
		 * @param v2 Second node
		 * @return The compare value
		 */
		static int compare(const node& v1, const node& v2) {
			OGDF_ASSERT(v1->graphOf() == v2->graphOf());
			return v1->degree() - v2->degree();
		}

		OGDF_AUGMENT_STATICCOMPARER(node);
	};
};
}
