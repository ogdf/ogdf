/** \file
 * \brief Applies the node coloring approximation specified by Johnson.
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
 * This class implements the approximation given by Johnson which tries
 * to find heuristically large independent sets one other the other.
 */
class NodeColoringJohnson : public NodeColoringModule {
public:
	/**
	 * Declares procedure to find the minimum degree nodes.
	 */
	enum class MinDegreeProcedure {
		smallestIndex, ///< Use the node with the smallest index
		minDegreeOriginal, ///< Use the node which has the smallest degree in the original graph
		maxDegreeOriginal, ///< Use the node which has the largest degree in the original graph
		minDegreeSubgraph, ///< Use the node which has the smallest degree in the current subgraph
		maxDegreeSubgraph ///< Use the node which has the largest degree in the current subgraph
	};

	/**
	 * The constructor.
	 * Initializes the procedure of finding minimum degree nodes with the smallest index.
	 */
	NodeColoringJohnson() : m_minDegreeProcedure(MinDegreeProcedure::maxDegreeSubgraph) { }

	/**
	 * Sets the procedure of finding minimum degree nodes.
	 * @param minDegreeProcedure The desired minimum degree finding procedure
	 */
	void setMinDegreeProcedure(MinDegreeProcedure minDegreeProcedure) {
		m_minDegreeProcedure = minDegreeProcedure;
	}

	virtual NodeColor call(const Graph& graph, NodeArray<NodeColor>& colors,
			NodeColor start = 0) override {
		auto numberColorsUsed = NodeColor(0);
		// Copy the input graph
		GraphCopy graphMain = GraphCopy(graph);
		preprocessGraph(graphMain);

		// Perform the coloring while the graph is not empty
		while (!graphMain.empty()) {
			// Search for big independent sets
			GraphCopy gSubgraph = graphMain;

			// Store the degrees of the
			NodeArray<int> degreesSubGraph(gSubgraph, 0);
			for (auto v : gSubgraph.nodes) {
				degreesSubGraph[v] = v->degree();
			}

			// Color the next independent set
			while (!gSubgraph.empty()) {
				// Search for all min degree nodes
				List<node> minDegreeNodes;
				getMinimumDegreeNodes(gSubgraph, minDegreeNodes);
				auto minDegreeNode = minDegreeNodes.front();

				// Select a specific min degree node
				if (m_minDegreeProcedure != MinDegreeProcedure::smallestIndex) {
					int extremalDegree;
					if (m_minDegreeProcedure == MinDegreeProcedure::minDegreeOriginal
							|| m_minDegreeProcedure == MinDegreeProcedure::maxDegreeOriginal) {
						extremalDegree = gSubgraph.original(minDegreeNode)->degree();
					} else if (m_minDegreeProcedure == MinDegreeProcedure::minDegreeSubgraph
							|| m_minDegreeProcedure == MinDegreeProcedure::maxDegreeSubgraph) {
						extremalDegree = degreesSubGraph[minDegreeNode];
					}

					for (node v : minDegreeNodes) {
						int degreeMainGraph = gSubgraph.original(v)->degree();
						int degreeSubGraph = degreesSubGraph[v];
						if (m_minDegreeProcedure == MinDegreeProcedure::minDegreeOriginal) {
							if (degreeMainGraph < extremalDegree) {
								extremalDegree = degreeMainGraph;
								minDegreeNode = v;
							}
						} else if (m_minDegreeProcedure == MinDegreeProcedure::maxDegreeOriginal) {
							if (degreeMainGraph > extremalDegree) {
								extremalDegree = degreeMainGraph;
								minDegreeNode = v;
							}
						} else if (m_minDegreeProcedure == MinDegreeProcedure::minDegreeSubgraph) {
							if (degreeSubGraph < extremalDegree) {
								extremalDegree = degreeSubGraph;
								minDegreeNode = v;
							}
						} else if (m_minDegreeProcedure == MinDegreeProcedure::maxDegreeSubgraph) {
							if (degreeSubGraph > extremalDegree) {
								extremalDegree = degreeSubGraph;
								minDegreeNode = v;
							}
						}
					}
				}

				// Color the node
				colors[gSubgraph.original(minDegreeNode)] = start + numberColorsUsed;

				// Delete the node and its neighbors from the subgraph
				List<node> nodes_to_delete;
				for (adjEntry adj : minDegreeNode->adjEntries) {
					nodes_to_delete.emplaceBack(adj->twinNode());
				}
				for (node v : nodes_to_delete) {
					gSubgraph.delNode(v);
				}

				node minDegOrigNode = gSubgraph.original(minDegreeNode);
				gSubgraph.delNode(minDegreeNode);

				// Delete the already colored node from the graph
				graphMain.delNode(graphMain.copy(minDegOrigNode));
			}
			// Increment the color the next independent set
			numberColorsUsed++;
		}
		// Check the coloring
		OGDF_ASSERT(checkColoring(graph, colors));
		return numberColorsUsed;
	}

private:
	MinDegreeProcedure m_minDegreeProcedure;
};
}