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

#pragma once

#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/graphalg/NodeColoringModule.h>
#include <ogdf/graphalg/NodeColoringSequential.h>
#include <ogdf/graphalg/NodeColoringSimple.h>

namespace ogdf {

/**
 * Approximation algorithms for the node coloring problem in graphs.
 * This class implements the approximation given by Wigderson which
 * colors the graph recursively until a bipartite graph is reached or
 * the given approximation ratio can be reached by a primitive coloring
 * algorithm.
 */
class NodeColoringWigderson : public NodeColoringModule {
public:
	/**
	 * Declares procedure to find the maximum degree nodes.
	 */
	enum class MaxDegreeProcedure {
		smallestIndex, ///< Use the node with the smallest index
		minDegreeRecursion, ///< Use the node which has the smallest degree in the recursive graph
		maxDegreeRecursion, ///< Use the node which has the largest degree in the recursive graph
		minDegreeOriginal, ///< Use the node which has the smallest degree in the original graph
		maxDegreeOriginal, ///< Use the node which has the largest degree in the original graph
		minDegreeNeighbors, ///< Use the node which neighbors have in sum the smallest degrees
		maxDegreeNeighbors, ///< Use the node which neighbors have in sum the largest degrees
	};

	/**
	 * Declares the procedure dealing with the recursion anchor coloring.
	 */
	enum class RecursionAnchorProcedure {
		trivialColoring, ///< Uses trivial coloring, i.e. colors each node with a different color
		sequentialColoringSimple, ///< Uses sequential coloring with the trivial node permutation
		sequentialColoringDegree, ///< Uses sequential coloring with the nodes sorted decreasing by their degree
		pooled, ///< Pools all nodes from the recursion anchor and then applies sequential coloring
	};

	/**
	 * Declares the procedure dealing with the brute force coloring.
	 */
	enum class BruteForceProcedure {
		sequentialColoringSimple, ///< Uses sequential coloring with the trivial node permutation
		sequentialColoringDegree, ///< Uses sequential coloring with the nodes sorted decreasing by their degree
		pooled, ///< Pools all nodes from the recursion anchor and then applies sequential coloring
	};

	/**
	 * The constructor.
	 * Initializes the search procedure with the Wigderson-Search by default.
	 * Initializes the procedure to find maximum degree nodes with the smallest index procedure.
	 */
	NodeColoringWigderson()
		: m_coloringSimple()
		, m_sequentialColoring()
		, m_searchProcedure(SearchProcedure::wigdersonSearch)
		, m_maxDegreeProcedure(MaxDegreeProcedure::smallestIndex)
		, m_recursionAnchorProcedure(RecursionAnchorProcedure::pooled)
		, m_bruteForceProcedure(BruteForceProcedure::pooled) { }

	/**
	 * Sets the search procedure to find the smallest possible parameter k such as
	 * the graph is k-colorable.
	 * @param searchProcedure The desired search procedure
	 */
	void setSearchProcedure(SearchProcedure searchProcedure) {
		m_searchProcedure = searchProcedure;
	}

	/**
	 * Sets the procedure of finding maximum degree nodes.
	 * @param maxDegreeProcedure The desired maximum degree finding procedure
	 */
	void setMaxDegreeProcedure(MaxDegreeProcedure maxDegreeProcedure) {
		m_maxDegreeProcedure = maxDegreeProcedure;
	}

	/**
	 * Sets the recursion anchor procedure.
	 * @param recursionAnchorProcedure The desired recursion anchor coloring procedure
	 */
	void setRecursionAnchorProcedure(RecursionAnchorProcedure recursionAnchorProcedure) {
		m_recursionAnchorProcedure = recursionAnchorProcedure;
	}

	/**
	 * Sets the brute force procedure.
	 * @param bruteForceProcedure The desired brute force coloring procedure
	 */
	void setBruteForceProcedure(BruteForceProcedure bruteForceProcedure) {
		m_bruteForceProcedure = bruteForceProcedure;
	}

	virtual NodeColor call(const Graph& graph, NodeArray<NodeColor>& colors,
			NodeColor start = 0) override {
		// Copy the input graph
		GraphCopy graphMain = GraphCopy(graph);
		preprocessGraph(graphMain);
		NodeArray<NodeColor> colorsMain(graphMain);

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
		NodeColor tmp = start;
		bool wigdersonResult = wigdersonCaller(graphMain, colorsMain, tmp, k);
		NodeColor numColors = getMaximumNodeColor(colorsMain) + 1 - start;
		start += numColors;
		OGDF_ASSERT(wigdersonResult);
		for (node v : graph.nodes) {
			colors[v] = colorsMain[graphMain.copy(v)];
		}
		return numColors;
	}

private:
	NodeColoringSimple m_coloringSimple;
	NodeColoringSequential m_sequentialColoring;
	SearchProcedure m_searchProcedure;
	MaxDegreeProcedure m_maxDegreeProcedure;
	RecursionAnchorProcedure m_recursionAnchorProcedure;
	BruteForceProcedure m_bruteForceProcedure;

	/**
	 * Calculates the function for the graph maximum degree
	 * specified by Wigderson.
	 *
	 * @param n Number of nodes.
	 * @param k Parameter, such that the graph is k-colorable.
	 * @return The function value.
	 */
	int wigdersonFunction(int n, int k) const {
		OGDF_ASSERT(n >= 0);
		OGDF_ASSERT(k >= 1);
		return std::ceil(std::pow(n, 1.0 - (1.0 / (static_cast<double>(k) - 1.0))));
	}

	/**
	 * Calls initially the recursive algorithm of Wigderson.
	 * The result is a certain graph coloring.
	 * @param graph The input graph to be colored.
	 * @param colors The resulting coloring.
	 * @param color The start color.
	 * @param k The parameter k such that the graph is k-colorable.
	 * @return True, iff a coloring was found for the given parameter k.
	 */
	bool wigdersonCaller(const Graph& graph, NodeArray<NodeColor>& colors, NodeColor& color, int k) {
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

	/**
	 * Calls the recursive algorithm of Wigderson to find a feasible graph coloring.
	 * @param graph The input graph to be colored.
	 * @param colors The resulting coloring.
	 * @param color The start color.
	 * @param k The parameter k such that the graph is k-colorable.
	 * @param degreesOriginal The node degrees in the initial original graph.
	 * @param nodesToBeColored Resulting list of nodes which need to be colored.
	 * @return True, iff a coloring was found for the given parameter k.
	 */
	bool wigdersonRecursive(const Graph& graph, NodeArray<NodeColor>& colors, NodeColor& color,
			int k, NodeArray<int>& degreesOriginal, List<node>& nodesToBeColored) {
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
			NodeArray<node> nodeTableOrig2New;
			inducedSubGraph<ListIterator<node>>(graphMain, neighbors.begin(), subGraph,
					nodeTableOrig2New);
			NodeArray<node> nodeTableNew2Orig;
			reverseNodeTable(graphMain, subGraph, nodeTableOrig2New, nodeTableNew2Orig);
			// Determine the original node degrees of the subgraph
			NodeArray<int> degreesOriginalSubgraph(subGraph);
			for (node v : neighbors) {
				degreesOriginalSubgraph[nodeTableOrig2New[v]] =
						degreesOriginal[graphMain.original(v)];
			}
			// Variables for the subgraph
			NodeArray<NodeColor> colorsSubGraph(subGraph);
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
		NodeArray<NodeColor> colorsCopy(graphMain);
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

	/**
	 * Wraps the Wigderson recursive algorithm.
	 */
	struct SearchWrapperWigderson : public SearchWrapper {
		/**
		 * Creates the wrapper.
		 * @param coloringWigderson Reference to the NodeColoringWigderson
		 * @param graph The graph to color
		 * @param colors The array of colors to be assigned
		 */
		SearchWrapperWigderson(NodeColoringWigderson& coloringWigderson, const Graph& graph,
				NodeArray<NodeColor>& colors)
			: m_coloring(coloringWigderson), m_graph(graph), m_colors(colors) { }

		bool step(int k) override {
			NodeColor start = NodeColor(1);
			return m_coloring.wigdersonCaller(m_graph, m_colors, start, k);
		}

		NodeColoringWigderson& m_coloring;
		const Graph& m_graph;
		NodeArray<NodeColor>& m_colors;
	};
};
}
