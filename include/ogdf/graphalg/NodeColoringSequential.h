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

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/comparer.h>
#include <ogdf/graphalg/NodeColoringModule.h>

namespace ogdf {
template<class E>
class List;

/**
 * Approximation algorithms for the node coloring problem in graphs.
 * This class applies the sequential coloring algorithm which greedily
 * assigns to each node the first available color.
 */
class OGDF_EXPORT NodeColoringSequential : public NodeColoringModule {
public:
	virtual NodeColor call(const Graph& graph, NodeArray<NodeColor>& colors,
			NodeColor start = 0) override {
		return colorByDegree(graph, colors, start);
	}

	/**
	 * Sorts a list of nodes decreasingly by the degree.
	 * @param nodes List of nodes to be sorted.
	 */
	virtual void sortByDegree(List<node>& nodes);

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
			NodeColor start = 0);

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
			NodeColor start = 0);

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
			List<node>& nodePermutation, NodeColor start = 0, bool lookForNeighbors = false);

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
