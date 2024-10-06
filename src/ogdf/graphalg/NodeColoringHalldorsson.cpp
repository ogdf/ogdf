/** \file
 * \brief Applies the node coloring approximation specified by Halldorsson.
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
#include <ogdf/graphalg/NodeColoringHalldorsson.h>
#include <ogdf/graphalg/NodeColoringModule.h>

#include <algorithm>
#include <cmath>

namespace ogdf {

using NColor = NodeColoringModule::NodeColor;

NColor NodeColoringHalldorsson::call(const Graph& graph, NodeArray<NColor>& colors, NColor start) {
	NColor numberOfColorsUsed = 0;
	// Copy the input graph
	GraphCopy graphMain = GraphCopy(graph);
	preprocessGraph(graphMain);
	NodeArray<NColor> colorsMain(graphMain);

	// Color each independent set until the graph is colored
	while (!graphMain.empty()) {
		List<node> halldorssonIndependentSet;
		performHalldorsson(graphMain, halldorssonIndependentSet);
		List<node> removelIndependentSet;
		cliqueRemoval(graphMain, removelIndependentSet);

		List<node>& largestIndependentSet = halldorssonIndependentSet;
		if (removelIndependentSet.size() > halldorssonIndependentSet.size()) {
			largestIndependentSet = removelIndependentSet;
		}

		for (node& v : largestIndependentSet) {
			colors[graphMain.original(v)] = start;
			graphMain.delNode(v);
		}

		start++;
		numberOfColorsUsed++;
	}

	OGDF_ASSERT(checkColoring(graph, colors));
	return numberOfColorsUsed;
}

void NodeColoringHalldorsson::performHalldorsson(const Graph& graph, List<node>& independentSet) {
	double alpha = 1.0;

	// Perform a search to find the smallest possible parameter k
	auto searchWrapper = SearchWrapperHalldorsson(*this, graph, independentSet, alpha);
	int k;
	switch (m_searchProcedure) {
	case SearchProcedure::linearSearch:
		k = searchLinear(&searchWrapper, 2, graph.numberOfNodes());
		break;
	case SearchProcedure::binarySearch:
		k = searchBinary(&searchWrapper, 2, graph.numberOfNodes());
		break;
	case SearchProcedure::wigdersonSearch:
		k = searchWigderson(&searchWrapper);
		break;
	default:
		k = searchWigderson(&searchWrapper);
		break;
	}

	// Perform the recursive algorithm with the smallest possible parameter k
	bool halldorssonResult = halldorssonRecursive(graph, independentSet, k, alpha);
	OGDF_ASSERT(halldorssonResult);
	OGDF_ASSERT(checkIndependentSet(graph, independentSet));
}

bool NodeColoringHalldorsson::halldorssonRecursive(const Graph& graph, List<node>& independentSet,
		int k, double alpha) {
	// 1. Check special easy cases
	OGDF_ASSERT(k >= 1);
	independentSet.clear();
	if (graph.numberOfNodes() <= 1) {
		for (node v : graph.nodes) {
			independentSet.emplaceBack(v);
		}
		return true;
	}

	// 2. Copy the input graph
	GraphCopy graphMain = GraphCopy(graph);
	preprocessGraph(graphMain);

	// 3. Prepare the main loop
	unsigned int n = graphMain.numberOfNodes();
	unsigned int m = std::max(
			std::min(std::floor(alpha * (std::log(n) / std::log(k))), static_cast<double>(n)), 1.0);
	Array<Array<node>> buckets;
	createBuckets(graphMain, std::min(k * m, n), buckets);
	for (Array<node>& bucket : buckets) {
		for (NodeColoringHalldorsson::SubsetIterator subsetIterator(bucket, m);
				subsetIterator.isOk(); subsetIterator.advance()) {
			independentSet.clear();
			auto subset = subsetIterator.currentSubset();

			// Check if the subset is independent
			if (checkIndependentSet(graphMain, subset)) {
				// Add the independent nodes to the result
				for (node& w : subset) {
					independentSet.emplaceBack(graphMain.original(w));
				}

				// Determine the complement neighbors of the independent set
				List<node> neighborsComplement;
				getNeighborsComplement<ListIterator<node>>(graphMain, subset.begin(),
						neighborsComplement);
				double size = static_cast<double>(n) / static_cast<double>(k) * std::log(n)
						/ (2.0 * std::log(std::log(n)));
				if (neighborsComplement.empty()) {
					return true;
				}

				// Determine the subgraph of the complement neighbors
				Graph subGraph;
				NodeArray<node> nodeTableOrig2New(graphMain, nullptr);
				EdgeArray<edge> edgeTableOrig2New(graphMain, nullptr);
				subGraph.insert(neighborsComplement, graphMain.edges, nodeTableOrig2New,
						edgeTableOrig2New);
				NodeArray<node> nodeTableNew2Orig;
				reverseNodeTable(graphMain, subGraph, nodeTableOrig2New, nodeTableNew2Orig);

				// Check the subgraph recursively for more independent sets
				if (neighborsComplement.size() >= size) {
					List<node> subIndependentSet;
					halldorssonRecursive(subGraph, subIndependentSet, k, alpha);
					for (node& w : subIndependentSet) {
						independentSet.emplaceBack(graphMain.original(nodeTableNew2Orig[w]));
					}
					return true;
				} else {
					List<node> subIndependentSet;
					cliqueRemoval(subGraph, subIndependentSet);
					int threshold =
							std::ceil(std::pow(std::log(n), 3.0) / (6.0 * std::log(std::log(n))));
					if (subIndependentSet.size() + independentSet.size() >= threshold) {
						for (node& w : subIndependentSet) {
							independentSet.emplaceBack(graphMain.original(nodeTableNew2Orig[w]));
						}
						return true;
					} else {
						continue;
					}
				}
			}
		}
	}
	// Failed to find such an independent set
	return false;
}
}
