/** \file
 * \brief Implementation of class SeparatorLiptonTarjan.
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

#include <ogdf/graphalg/SeparatorLiptonTarjan.h>

namespace ogdf {

// connects this instance to a graph by initializing arrays etc.
void SeparatorLiptonTarjan::makeTree() {
	tree = std::make_shared<BFSTreeClassical>(*graph, getStartNode(*graph), treeHeightIterations,
			false);
}

bool SeparatorLiptonTarjan::doSeparate(const Graph& G, List<node>& separator, List<node>& first,
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
	tree->restructure(separator, second, useTriBFS);

	if (!useTriBFS) {
		// embed and triangulate the planar embedding of the graph
		planarEmbedPlanarGraph(
				*graph); // re-planarization necessary because restructuring broke embedding
		triangulate(*graph);
	}

	// choose any non-tree edge as start-Edge for new cycle
	Cycle cycle(tree.get(), chooseEdge());

	while (cycle.getInsideCost() > 2.0 / 3.0 * graph->nodes.size()) {
		cycle = cycle.expandCycle();
	}

	exitPoint = "cycle_expansion";
	cycle.fillLists(separator, first, second);

	return true;
}

void SeparatorLiptonTarjan::fillLists(List<node>& separator, List<node>& first,
		List<node>& second) const {
	// fill separator with levels t0 and t2
	List<node> t0_level = tree->getLevel(tree->get_t0());
	List<node> t2_level = tree->getLevel(tree->get_t2());
	t0_level.conc(t2_level);
	for (const node no : t0_level) {
		separator.pushBack(graph->original(no));
	}

	// construct the two parts of the separation: inner and outer
	List<node> outer = tree->getNodesFromTo(0, tree->get_t0());
	List<node> inner = tree->getNodesFromTo(tree->get_t0() + 1, tree->get_t2());
	List<node> upper = tree->getNodesFrom(tree->get_t2() + 1);
	outer.conc(upper);

	// find the biggest chunk and add it to first
	if (outer.size() > inner.size()) {
		std::swap(inner, outer);
	}

	for (const node no : inner) {
		first.pushBack(graph->original(no));
	}
	for (const node no : outer) {
		second.pushBack(graph->original(no));
	}
}

edge SeparatorLiptonTarjan::chooseEdge() const {
	for (edge e : graph->edges) {
		if (!tree->isInTree(e)) {
			return e;
		}
	}
	// Unable to find start edge!
	OGDF_ASSERT(false);
	return graph->firstEdge();
}

}
