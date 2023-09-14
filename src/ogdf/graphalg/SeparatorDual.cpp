/** \file
 * \brief Implementation of class SeparatorDual.
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

#include <ogdf/graphalg/SeparatorDual.h>

namespace ogdf {

void SeparatorDual::makeTree() {
	tree = std::make_shared<BFSTreeClassical>(*graph, getStartNode(*graph), treeHeightIterations,
			true);
}

bool SeparatorDual::doSeparate(const Graph& G, List<node>& separator, List<node>& first,
		List<node>& second) {
	// construct BFS Tree
	makeTree();

	// if a level is small enough and works as separator, the tree is aborted and the level is returned here
	// separatorLevel contains the lowest level that fulfilled the criteria or -1 if none did
	int separatorLevel = tree->getSeparatorLevel();

	if (separatorLevel != -1) {
		exitPoint = "bfs_level";
		for (node no : tree->getLevel(separatorLevel)) {
			separator.pushBack(graph->original(no));
			graph->delNode(no);
		}
		return separateComponents(*graph, separator, first, second, true);
	}

	// check if the two levels work as a separator already
	int sum = 0;
	for (int i = tree->get_t0() + 1; i < tree->get_t2(); i++) {
		sum += tree->getSizeOfLevel(i);
	}

	// if successful, return separator
	if (sum <= 2.0 / 3.0 * graph->nodes.size()) {
		exitPoint = "two_bfs_levels";
		fillLists(separator, first, second);
		return true;
	}

	// restructure the tree: collapse t0 and below, remove t2 and above, create new root
	// if tri-BFS is used, this already produces a planarly embedded, triangulated graph
	tree->restructure(separator, second, useTriBFS);

	if (!useTriBFS) {
		// embed and triangulate the planar embedding of the graph
		planarEmbedPlanarGraph(*graph); // re-embedding because restructuring destroyed planarity
		triangulate(*graph);
	}

	// start depth first search over dual graph
	planar_separator::SeparatorDualHelper helper(graph, tree);
	planar_separator::SeparatorDualHelper::CycleData cycle = helper.dfs();

	if (cycle.isInCycle(tree->getRoot())) {
		cycle.cycle.removeFirst(tree->getRoot());
	}

	graph->delNode(tree->getRoot());

	for (node n : cycle.cycle) {
		separator.pushBack(graph->original(n));
		graph->delNode(n);
	}
	exitPoint = "cycle_expansion";
	return separateComponents(*graph, separator, first, second, true);
}

}
