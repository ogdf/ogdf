/** \file
 * \brief Applies the node coloring heuristic "Recursive Largest First"
 * specified by Leighton (1979).
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

/**
 * A simple greedy node coloring heuristic in graphs. It colors one node
 * after another by preferring nodes with a large degree.
 */
class OGDF_EXPORT NodeColoringRecursiveLargestFirst : public NodeColoringModule {
public:
	virtual NodeColor call(const Graph& graph, NodeArray<NodeColor>& colors,
			NodeColor start = 0) override;

private:
	/**
	 * Searches a node candidate for the next coloring step
	 * @param candidate The resulting node candidate
	 * @param degreesUnavailable The degrees to the unavailable nodes
	 * @param graph The graph of the nodes
	 */
	void getCandidate(node& candidate, NodeArray<int>& degreesUnavailable, Graph& graph);
};
}
