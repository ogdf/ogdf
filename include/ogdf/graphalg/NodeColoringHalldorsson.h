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

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/basic.h>
#include <ogdf/graphalg/NodeColoringModule.h>

namespace ogdf {
template<class E>
class List;

/**
 * Approximation algorithms for the node coloring problem in graphs.
 * This class implements the approximation given by Halldorsson which
 * colors the graph by finding independent sets.
 */
class OGDF_EXPORT NodeColoringHalldorsson : public NodeColoringModule {
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
	inline void setSearchProcedure(SearchProcedure searchProcedure) {
		m_searchProcedure = searchProcedure;
	}

	virtual NodeColor call(const Graph& graph, NodeArray<NodeColor>& colors,
			NodeColor start = 0) override;

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
	void performHalldorsson(const Graph& graph, List<node>& independentSet);

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
	bool halldorssonRecursive(const Graph& graph, List<node>& independentSet, int k, double alpha);

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
