/** \file
 * \brief Applies the node coloring approximation specified by Wigderson.
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
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/graphalg/NodeColoringModule.h>
#include <ogdf/graphalg/NodeColoringSequential.h>
#include <ogdf/graphalg/NodeColoringSimple.h>
#include <ogdf/graphalg/NodeColoringWigderson.h>

#include <cmath>

namespace ogdf {

using NColor = NodeColoringModule::NodeColor;

NColor NodeColoringWigderson::call(const Graph& graph, NodeArray<NColor>& colors, NColor start) {
	if (graph.empty()) {
		return 0;
	}

	// Copy the input graph
	GraphCopy graphMain = GraphCopy(graph);
	preprocessGraph(graphMain);
	NodeArray<NColor> colorsMain(graphMain);

	// Perform a search to find the smallest possible parameter k
	auto searchWrapper = SearchWrapperWigderson(*this, graphMain, colorsMain);
	int k;
	switch (m_searchProcedure) {
	case SearchProcedure::linearSearch:
		k = searchLinear(&searchWrapper, 2, graphMain.numberOfNodes());
		break;
	case SearchProcedure::binarySearch:
		k = searchBinary(&searchWrapper, 2, graphMain.numberOfNodes());
		break;

	case SearchProcedure::wigdersonSearch:
		k = searchWigderson(&searchWrapper);
		break;

	default:
		k = searchWigderson(&searchWrapper);
		break;
	}

	// Perform Wigderson with the smallest possible k
	NColor tmp = start;
	bool wigdersonResult = wigdersonCaller(graphMain, colorsMain, tmp, k);
	NColor numColors = getMaximumNodeColor(colorsMain) + 1 - start;
	start += numColors;
	OGDF_ASSERT(wigdersonResult);
	for (node v : graph.nodes) {
		colors[v] = colorsMain[graphMain.copy(v)];
	}
	return numColors;
}

bool NodeColoringWigderson::wigdersonCaller(const Graph& graph, NodeArray<NColor>& colors,
		NColor& color, int k) {
	// Store the starting color
	auto startColor = color;

	// Determine the degrees in the original graph
	NodeArray<int> degreesOriginal(graph);
	for (node v : graph.nodes) {
		degreesOriginal[v] = v->degree();
	}

	// List of nodes which are to be colored
	List<node> nodesToBeColored;

	// Apply the recursive Wigderson algorithm
	auto result = wigdersonRecursive(graph, colors, color, k, degreesOriginal, nodesToBeColored);

	// Check the result
	if (result) {
		m_sequentialColoring.fromPermutation(graph, colors, nodesToBeColored, startColor, true);
		auto maxColorIndex = getMaximumNodeColor(colors);
		color = maxColorIndex + 1;
		OGDF_ASSERT(checkColoring(graph, colors));
	}
	return result;
}

bool NodeColoringWigderson::wigdersonRecursive(const Graph& graph, NodeArray<NColor>& colors,
		NColor& color, int k, NodeArray<int>& degreesOriginal, List<node>& nodesToBeColored) {
	// 1. Preparation
	OGDF_ASSERT(k > 0);
	GraphCopy graphMain = GraphCopy(graph);
	int n = graphMain.numberOfNodes();

	// 2. Simple cases (recursion anchor)
	// Check if the graph is bipartite
	if (k <= 2) {
		NodeArray<bool> partitions(graph);
		bool bipartite = isBipartite(graph, partitions);
		if (!bipartite) {
			return false;
		}
		bool partitionUsedA = false;
		bool partitionUsedB = false;
		for (node v : graph.nodes) {
			if (partitions[v]) {
				colors[v] = color;
				partitionUsedA = true;
			} else {
				colors[v] = color + 1;
				partitionUsedB = true;
			}
		}
		color += partitionUsedA + partitionUsedB;
		return true;
	} else if (k >= std::log2(n)) {
		switch (m_recursionAnchorProcedure) {
		case RecursionAnchorProcedure::trivialColoring:
			color += m_coloringSimple.call(graph, colors, color);
			break;
		case RecursionAnchorProcedure::sequentialColoringSimple:
			color += m_sequentialColoring.colorByIndex(graph, colors, color);
			break;
		case RecursionAnchorProcedure::sequentialColoringDegree:
			color += m_sequentialColoring.colorByDegree(graph, colors, color);
			break;
		case RecursionAnchorProcedure::pooled:
			for (node v : graph.nodes) {
				nodesToBeColored.pushBack(v);
			}
			m_sequentialColoring.sortByDegree(nodesToBeColored);
			break;
		}
		return true;
	}

	// 3. Recursive coloring stage
	node maxDegreeNode = graphMain.firstNode();
	while (!graphMain.empty()
			&& (getMaximumDegreeNode(graphMain, maxDegreeNode) >= wigdersonFunction(n, k))) {
		// Get a good fitting maximum degree node
		if (m_maxDegreeProcedure != MaxDegreeProcedure::smallestIndex) {
			List<node> maxDegreeNodes;
			getMaximumDegreeNodes(graphMain, maxDegreeNodes);
			maxDegreeNode = maxDegreeNodes.front();

			int extremalDegree;
			if (m_maxDegreeProcedure == MaxDegreeProcedure::minDegreeRecursion
					|| m_maxDegreeProcedure == MaxDegreeProcedure::maxDegreeRecursion) {
				extremalDegree = graphMain.original(maxDegreeNode)->degree();
			} else if (m_maxDegreeProcedure == MaxDegreeProcedure::minDegreeOriginal
					|| m_maxDegreeProcedure == MaxDegreeProcedure::maxDegreeOriginal) {
				extremalDegree = degreesOriginal[graphMain.original(maxDegreeNode)];
			} else { // MaxDegreeProcedure::min/maxDegreeNeighbors
				extremalDegree = getNeighborDegrees(maxDegreeNode);
			}
			for (node v : maxDegreeNodes) {
				int degreeRecursiveGraph = graphMain.original(v)->degree();
				int degreeOriginalGraph = degreesOriginal[graphMain.original(v)];
				int degreeNeighbors = getNeighborDegrees(v);
				if (m_maxDegreeProcedure == MaxDegreeProcedure::minDegreeRecursion) {
					if (degreeRecursiveGraph < extremalDegree) {
						extremalDegree = degreeRecursiveGraph;
						maxDegreeNode = v;
					}
				} else if (m_maxDegreeProcedure == MaxDegreeProcedure::maxDegreeRecursion) {
					if (degreeRecursiveGraph > extremalDegree) {
						extremalDegree = degreeRecursiveGraph;
						maxDegreeNode = v;
					}
				} else if (m_maxDegreeProcedure == MaxDegreeProcedure::minDegreeOriginal) {
					if (degreeOriginalGraph < extremalDegree) {
						extremalDegree = degreeOriginalGraph;
						maxDegreeNode = v;
					}
				} else if (m_maxDegreeProcedure == MaxDegreeProcedure::maxDegreeOriginal) {
					if (degreeOriginalGraph > extremalDegree) {
						extremalDegree = degreeOriginalGraph;
						maxDegreeNode = v;
					}
				} else if (m_maxDegreeProcedure == MaxDegreeProcedure::minDegreeNeighbors) {
					if (degreeNeighbors < extremalDegree) {
						extremalDegree = degreeNeighbors;
						maxDegreeNode = v;
					}
				} else if (m_maxDegreeProcedure == MaxDegreeProcedure::maxDegreeNeighbors) {
					if (degreeNeighbors > extremalDegree) {
						extremalDegree = degreeNeighbors;
						maxDegreeNode = v;
					}
				}
			}
		}

		// Get the neighbor nodes
		List<node> neighbors;
		for (adjEntry adj : maxDegreeNode->adjEntries) {
			neighbors.emplaceBack(adj->twinNode());
		}
		// Create subgraph
		Graph subGraph;
		NodeArray<node> nodeTableOrig2New(graphMain, nullptr);
		EdgeArray<edge> edgeTableOrig2New(graphMain, nullptr);
		subGraph.insert(neighbors, graphMain.edges, nodeTableOrig2New, edgeTableOrig2New);
		NodeArray<node> nodeTableNew2Orig;
		reverseNodeTable(graphMain, subGraph, nodeTableOrig2New, nodeTableNew2Orig);
		// Determine the original node degrees of the subgraph
		NodeArray<int> degreesOriginalSubgraph(subGraph);
		for (node v : neighbors) {
			degreesOriginalSubgraph[nodeTableOrig2New[v]] = degreesOriginal[graphMain.original(v)];
		}
		// Variables for the subgraph
		NodeArray<NColor> colorsSubGraph(subGraph);
		List<node> subgraphNodesToBeColored;
		// Recursive algorithm call
		if (!wigdersonRecursive(subGraph, colorsSubGraph, color, k - 1, degreesOriginalSubgraph,
					subgraphNodesToBeColored)) {
			return false;
		}
		// Store the remaining nodes to be colored
		for (node toBeColored : subgraphNodesToBeColored) {
			nodesToBeColored.pushBack(graphMain.original(nodeTableNew2Orig[toBeColored]));
		}
		// Delete the already colored nodes
		colors[graphMain.original(maxDegreeNode)] = color++;
		for (node v : neighbors) {
			colors[graphMain.original(v)] = colorsSubGraph[nodeTableOrig2New[v]];
			graphMain.delNode(v);
		}
		graphMain.delNode(maxDegreeNode);
	}

	// 4. Brute force coloring stage
	NodeArray<NColor> colorsCopy(graphMain);
	switch (m_bruteForceProcedure) {
	case BruteForceProcedure::sequentialColoringSimple:
		color += m_sequentialColoring.colorByIndex(graphMain, colorsCopy, color);
		break;
	case BruteForceProcedure::sequentialColoringDegree:
		color += m_sequentialColoring.colorByDegree(graphMain, colorsCopy, color);
		break;
	case BruteForceProcedure::pooled:
		for (node v : graphMain.nodes) {
			nodesToBeColored.emplaceBack(graphMain.original(v));
			m_sequentialColoring.sortByDegree(nodesToBeColored);
		}
		break;
	}
	for (node v : graphMain.nodes) {
		colors[graphMain.original(v)] = colorsCopy[v];
	}

	// Return statement
	return true;
}

}
