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

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/basic.h>
#include <ogdf/graphalg/NodeColoringModule.h>
#include <ogdf/graphalg/NodeColoringSequential.h>
#include <ogdf/graphalg/NodeColoringSimple.h>

#include <cmath>

namespace ogdf {
template<class E>
class List;

/**
 * Approximation algorithms for the node coloring problem in graphs.
 * This class implements the approximation given by Wigderson which
 * colors the graph recursively until a bipartite graph is reached or
 * the given approximation ratio can be reached by a primitive coloring
 * algorithm.
 */
class OGDF_EXPORT NodeColoringWigderson : public NodeColoringModule {
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
	inline void setSearchProcedure(SearchProcedure searchProcedure) {
		m_searchProcedure = searchProcedure;
	}

	/**
	 * Sets the procedure of finding maximum degree nodes.
	 * @param maxDegreeProcedure The desired maximum degree finding procedure
	 */
	inline void setMaxDegreeProcedure(MaxDegreeProcedure maxDegreeProcedure) {
		m_maxDegreeProcedure = maxDegreeProcedure;
	}

	/**
	 * Sets the recursion anchor procedure.
	 * @param recursionAnchorProcedure The desired recursion anchor coloring procedure
	 */
	inline void setRecursionAnchorProcedure(RecursionAnchorProcedure recursionAnchorProcedure) {
		m_recursionAnchorProcedure = recursionAnchorProcedure;
	}

	/**
	 * Sets the brute force procedure.
	 * @param bruteForceProcedure The desired brute force coloring procedure
	 */
	inline void setBruteForceProcedure(BruteForceProcedure bruteForceProcedure) {
		m_bruteForceProcedure = bruteForceProcedure;
	}

	virtual NodeColor call(const Graph& graph, NodeArray<NodeColor>& colors,
			NodeColor start = 0) override;

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
	inline int wigdersonFunction(int n, int k) const {
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
	bool wigdersonCaller(const Graph& graph, NodeArray<NodeColor>& colors, NodeColor& color, int k);

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
			int k, NodeArray<int>& degreesOriginal, List<node>& nodesToBeColored);

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
