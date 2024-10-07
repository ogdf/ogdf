/** \file
 * \brief Applies the node coloring approximation specified by Berger&Rompel.
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
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/graphalg/NodeColoringBergerRompel.h>
#include <ogdf/graphalg/NodeColoringJohnson.h>
#include <ogdf/graphalg/NodeColoringModule.h>
#include <ogdf/graphalg/NodeColoringSequential.h>
#include <ogdf/graphalg/NodeColoringSimple.h>
#include <ogdf/graphalg/NodeColoringWigderson.h>

#include <algorithm>
#include <cmath>

namespace ogdf {

using NColor = NodeColoringModule::NodeColor;

NColor NodeColoringBergerRompel::call(const Graph& graph, NodeArray<NColor>& colors, NColor start) {
	// Copy the input graph
	GraphCopy graphMain = GraphCopy(graph);
	preprocessGraph(graphMain);
	NodeArray<NColor> colorsMain(graphMain);

	// Perform a search to find the smallest possible parameter k
	auto searchWrapper = SearchWrapperBergerRompel(*this, graphMain, colorsMain, m_alpha);
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

	// Perform Berger&Rompel with the smallest possible k
	NColor numberColorsBegin = start;
	bool bergerRompelResult = bergerRompelParameterized(graphMain, colorsMain, start, k, m_alpha);
	OGDF_ASSERT(bergerRompelResult);
	for (node v : graph.nodes) {
		colors[v] = colorsMain[graphMain.copy(v)];
	}
	return start - numberColorsBegin;
}

bool NodeColoringBergerRompel::bergerRompelParameterized(const Graph& graph,
		NodeArray<NColor>& colors, NColor& color, int k, double alpha) {
	// 1. Check special easy cases
	OGDF_ASSERT(k >= 1);
	auto startColor = color;
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
		OGDF_ASSERT(checkColoring(graph, colors));
		return true;
	}

	// 2. Check if Johnson algorithm is sufficient
	const int n = graph.numberOfNodes();
	if (k > min(std::sqrt(n), std::pow(n, alpha))) {
		color += m_johnsonColoring.call(graph, colors, color);
		return true;
	}

	// 3. Copy the input graph
	GraphCopy graphMain = GraphCopy(graph);
	preprocessGraph(graphMain);

	// 4. Prepare the main loop
	const int m = std::floor(alpha * (std::log(n) / std::log(k)));

	// 5. Perform the coloring until the number of nodes is small enough
	while (graphMain.numberOfNodes() >= k * m) {
		// Search for big independent sets
		GraphCopy gSubgraph = graphMain;

		while (gSubgraph.numberOfNodes() >= k * m) {
			Array<Array<node>> buckets;
			createBuckets(gSubgraph, k * m, buckets);
			auto finished = false;
			for (Array<node>& bucket : buckets) {
				for (NodeColoringBergerRompel::SubsetIterator subsetIterator(bucket, m);
						subsetIterator.isOk(); subsetIterator.advance()) {
					auto subset = subsetIterator.currentSubset();
					List<node> neighbors;
					getNeighbors<ListIterator<node>>(gSubgraph, subset.begin(), neighbors);
					auto numberNodesSubgraph = gSubgraph.numberOfNodes();
					if (checkIndependentSet(gSubgraph, subset)
							&& neighbors.size() <= (numberNodesSubgraph - numberNodesSubgraph / k)) {
						// Color the nodes in the subset
						for (node v : subset) {
							colors[gSubgraph.original(v)] = color;
						}

						List<node> nodesToDelete;
						mergeNodeLists<ListIterator<node>>(gSubgraph, subset.begin(),
								neighbors.begin(), nodesToDelete);

						// Delete the already colored nodes from the original graph
						for (node v : subset) {
							graphMain.delNode(graphMain.copy(gSubgraph.original(v)));
						}

						// Delete nodes in the subgraph
						for (node v : nodesToDelete) {
							gSubgraph.delNode(v);
						}

						finished = true;
						break;
					}
				}
				if (finished) {
					break;
				}
			}
			if (!finished) {
				return false;
			}
		}
		// Increment the color the next independent set
		color++;
	}

	// 6. Brute force coloring stage
	NodeArray<NColor> colorsCopy(graphMain);
	// Case distinction
	switch (m_bruteForceProcedure) {
	case BruteForceProcedure::trivialColoring:
		color += m_simpleColoring.call(graphMain, colorsCopy, color);
		break;
	case BruteForceProcedure::sequentialColoringSimple:
		color += m_sequentialColoring.colorByIndex(graphMain, colorsCopy, color);
		break;
	case BruteForceProcedure::sequentialColoringDegree:
		color += m_sequentialColoring.colorByDegree(graphMain, colorsCopy, color);
		break;
	case BruteForceProcedure::johnsonColoring:
		color += m_johnsonColoring.call(graphMain, colorsCopy, color);
		break;
	case BruteForceProcedure::wigdersonColoring:
		color += m_wigdersonColoring.call(graphMain, colorsCopy, color);
		break;
	case BruteForceProcedure::pooled:
		List<node> nodesToBeColored;
		for (node v : graphMain.nodes) {
			nodesToBeColored.emplaceBack(graphMain.original(v));
		}
		m_sequentialColoring.sortByDegree(nodesToBeColored);
		m_sequentialColoring.fromPermutation(graph, colors, nodesToBeColored, startColor, true);
		auto maxColorIndex = getMaximumNodeColor(colors);
		color = maxColorIndex + 1 - startColor;
		break;
	}
	// Copy the coloring result to the original graph
	if (m_bruteForceProcedure != BruteForceProcedure::pooled) {
		for (node v : graphMain.nodes) {
			colors[graphMain.original(v)] = colorsCopy[v];
		}
	}

	// Check the coloring
	OGDF_ASSERT(checkColoring(graph, colors));
	return true;
}

}
