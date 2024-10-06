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

#pragma once

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/basic.h>
#include <ogdf/graphalg/NodeColoringJohnson.h>
#include <ogdf/graphalg/NodeColoringModule.h>
#include <ogdf/graphalg/NodeColoringSequential.h>
#include <ogdf/graphalg/NodeColoringSimple.h>
#include <ogdf/graphalg/NodeColoringWigderson.h>

namespace ogdf {

/**
 * Approximation algorithms for the node coloring problem in graphs.
 * This class implements the approximation given by Berger&Rompel which
 * colors the graph by finding independent sets.
 */
class OGDF_EXPORT NodeColoringBergerRompel : public NodeColoringModule {
public:
	/**
	 * Declares the procedure dealing with the brute force coloring.
	 */
	enum class BruteForceProcedure {
		trivialColoring, ///< Uses trivial coloring, i.e. colors each node with a different color
		sequentialColoringSimple, ///< Uses sequential coloring with the trivial node permutation
		sequentialColoringDegree, ///< Uses sequential coloring with the nodes sorted decreasing by their degree
		johnsonColoring, ///< Uses Johnson's coloring algorithm
		wigdersonColoring, ///< Uses Wigderson's coloring algorithm
		pooled, ///< Pools all nodes from the recursion anchor and then applies sequential coloring
	};

	/**
	 * The constructor.
	 * Initializes the search procedure with the binary-search by default.
	 * The parameter alpha is initialized with 1.0 by default.
	 */
	NodeColoringBergerRompel()
		: m_simpleColoring()
		, m_sequentialColoring()
		, m_johnsonColoring()
		, m_wigdersonColoring()
		, m_searchProcedure(SearchProcedure::wigdersonSearch)
		, m_bruteForceProcedure(BruteForceProcedure::pooled)
		, m_alpha(0.4) { }

	/**
	 * Sets the search procedure to find the smallest possible parameter k such as
	 * the graph is k-colorable.
	 * @param searchProcedure The desired search procedure
	 */
	inline void setSearchProcedure(SearchProcedure searchProcedure) {
		m_searchProcedure = searchProcedure;
	}

	/**
	 * Sets the brute force procedure.
	 * @param bruteForceProcedure The desired brute force coloring procedure
	 */
	inline void setBruteForceProcedure(BruteForceProcedure bruteForceProcedure) {
		m_bruteForceProcedure = bruteForceProcedure;
	}

	/**
	 * Sets the parameter alpha which controls the accuracy of the algorithm.
	 * A higher value of alpha results in a better approximation ratio.
	 * Also, the running time increases with alpha.
	 * @param alpha The control parameter alpha
	 */
	inline void setAlpha(double alpha) {
		OGDF_ASSERT(alpha > 0.1);
		m_alpha = alpha;
	}

	virtual NodeColor call(const Graph& graph, NodeArray<NodeColor>& colors,
			NodeColor start = 0) override;

private:
	NodeColoringSimple m_simpleColoring;
	NodeColoringSequential m_sequentialColoring;
	NodeColoringJohnson m_johnsonColoring;
	NodeColoringWigderson m_wigdersonColoring;
	SearchProcedure m_searchProcedure;
	BruteForceProcedure m_bruteForceProcedure;
	double m_alpha;

	bool bergerRompelParameterized(const Graph& graph, NodeArray<NodeColor>& colors,
			NodeColor& color, int k, double alpha);

	/**
	 * Wraps the parameterized Berger&Rompel algorithm
	 */
	struct SearchWrapperBergerRompel : public SearchWrapper {
		/**
		 * Creates the wrapper.
		 * @param coloringBergerRompel Reference to the NodeColoringBergerRompel
		 * @param graph The graph to color
		 * @param colors The array of colors to be assigned
		 * @param alpha The alpha control parameter of the Berger&Rompel algorithm
		 */
		SearchWrapperBergerRompel(NodeColoringBergerRompel& coloringBergerRompel,
				const Graph& graph, NodeArray<NodeColor>& colors, double alpha)
			: m_coloring(coloringBergerRompel), m_graph(graph), m_colors(colors), m_alpha(alpha) { }

		bool step(int k) override {
			NodeColor start = NodeColor(1);
			return m_coloring.bergerRompelParameterized(m_graph, m_colors, start, k, m_alpha);
		}

		NodeColoringBergerRompel& m_coloring;
		const Graph& m_graph;
		NodeArray<NodeColor>& m_colors;
		double m_alpha;
	};
};
}
