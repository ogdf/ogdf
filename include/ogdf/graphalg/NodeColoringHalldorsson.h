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

#pragma once

#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/graphalg/NodeColoringJohnson.h>
#include <ogdf/graphalg/NodeColoringModule.h>
#include <ogdf/graphalg/NodeColoringSequential.h>

namespace ogdf {

/**
 * Approximation algorithms for the node coloring problem in graphs.
 * This class implements the approximation given by Halldorsson which
 * colors the graph by finding independent sets.
 */
class NodeColoringHalldorsson : public NodeColoringModule {
public:
	/**
	 * The constructor.
	 * Initializes the search procedure with the linear search by default.
	 */
	NodeColoringHalldorsson() { m_searchProcedure = SearchProcedure::wigdersonSearch; }

	/**
	 * Sets the search procedure to find the smallest possible parameter k such as
	 * the graph is k-colorable.
	 * @param searchProcedure The desired search procedure
	 */
	void setSearchProcedure(SearchProcedure searchProcedure) {
		m_searchProcedure = searchProcedure;
	}

	virtual NodeColor call(const Graph& graph, NodeArray<NodeColor>& colors,
			NodeColor start = 0) override {
		NodeColor numberOfColorsUsed = 0;
		// Copy the input graph
		GraphCopy graphMain = GraphCopy(graph);
		preprocessGraph(graphMain);
		NodeArray<NodeColor> colorsMain(graphMain);

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

private:
	SearchProcedure m_searchProcedure;

	/**
	 * Performs the Halldorsson algorithm for finding an independent set.
	 * This functions calls the recursive algorithm with varying parameters
	 * such that the biggest possible independent set will be found.
	 *
	 * @param graph The input graph
	 * @param independentSet The resulting independent set
	 */
	void performHalldorsson(const Graph& graph, List<node>& independentSet) {
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

	/**
	 * Performs the recursive Halldorsson algorithm to find large
	 * independent sets.
	 *
	 * @param graph The input graph
	 * @param independentSet The resulting independent set
	 * @param k Value, such that the graph is k-colorable
	 * @param alpha Control parameter for the proximity
	 * @return If the given proximity can be reached for the given k-value
	 */
	bool halldorssonRecursive(const Graph& graph, List<node>& independentSet, int k, double alpha) {
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
				std::min(std::floor(alpha * (std::log(n) / std::log(k))), static_cast<double>(n)),
				1.0);
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
						int threshold = std::ceil(
								std::pow(std::log(n), 3.0) / (6.0 * std::log(std::log(n))));
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

	/**
	 * Wraps the recursive Halldórsson algorithm
	 */
	struct SearchWrapperHalldorsson : public SearchWrapper {
		/**
		 * Creates the wrapper.
		 * @param coloringHalldorsson Reference to the NodeColoringHalldorsson
		 * @param graph The graph to search an independent set
		 * @param independentSet The resulting independent set
		 * @param alpha Control parameter alpha of the Halldorsson algorithm
		 */
		SearchWrapperHalldorsson(NodeColoringHalldorsson& coloringHalldorsson, const Graph& graph,
				List<node>& independentSet, double alpha)
			: m_coloring(coloringHalldorsson)
			, m_graph(graph)
			, m_independentSet(independentSet)
			, m_alpha(alpha) { }

		bool step(int k) override {
			return m_coloring.halldorssonRecursive(m_graph, m_independentSet, k, m_alpha);
		}

		NodeColoringHalldorsson& m_coloring;
		const Graph& m_graph;
		List<node>& m_independentSet;
		double m_alpha;
	};
};
}
