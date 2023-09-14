/** \file
 * \brief Implementation of class SeparatorDualFC.
 *
 * \author Thomas Klein
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

#include <ogdf/graphalg/SeparatorDualFC.h>

namespace ogdf {

void SeparatorDualFC::makeTree() {
	if (useTriangulatingBFS) {
		// construct triangulating bfs tree
		tree = std::make_shared<TriangulatingBFSTree>(*graph, getStartNode(*graph));
	} else {
		// embed and triangulate the planar embedding of the graph
		triangulate(*graph);
		// construct BFS Tree
		tree = std::make_shared<BFSTreeFC>(*graph, getStartNode(*graph));
	}
#ifdef OGDF_HEAVY_DEBUG
	tree->assertTreeConsistency();
#endif
}

bool SeparatorDualFC::doSeparate(const Graph& G, List<node>& separator, List<node>& first,
		List<node>& second) {
	makeTree();

	return findCycle(separator, first, second);
}

bool SeparatorDualFC::findCycle(List<node>& separator, List<node>& first, List<node>& second) {
	// start depth first search over dual graph
	planar_separator::SeparatorDualHelper helper(graph, tree);
	planar_separator::SeparatorDualHelper::CycleData cycle = helper.dfs();

	for (node n : cycle.cycle) {
		separator.pushBack(graph->original(n));
		graph->delNode(n);
	}
	exitPoint = "cycle_expansion";
	return separateComponents(*graph, separator, first, second, true);
}

}
