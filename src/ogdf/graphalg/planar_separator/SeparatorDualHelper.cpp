/** \file
 * \brief Implementation of class SeparatorDualHelper.
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

#include <ogdf/graphalg/planar_separator/SeparatorDualHelper.h>

namespace ogdf {

namespace planar_separator {

SeparatorDualHelper::CycleData SeparatorDualHelper::dfs() {
	embedding.init(*graph);
	marked.init(embedding, false);

	// DFS over faces of the embedding = dual of the graph
	// choose face right next to root node: has only one neighbour, as both adJEntries should be tree edges
	adjEntry startAdj = *tree->getRoot()->adjEntries.begin();
	while (!tree->isInTree(startAdj->theEdge())) {
		startAdj = startAdj->cyclicSucc();
	}
	face rootFace = embedding.rightFace(startAdj);
	marked[rootFace] = true;

	adjEntry adj = *rootFace->entries.begin();
	while (!tree->isInTree(adj->theEdge())) {
		adj = adj->faceCycleSucc();
	}

	CycleData cycleData = process(rootFace, adj->twin());

	return cycleData;
}

SeparatorDualHelper::CycleData SeparatorDualHelper::process(face f, adjEntry adj) {
	OGDF_ASSERT(embedding.leftFace(adj) == f);
	adjEntry twin = adj->twin();

	marked[f] = true;
	List<std::pair<face, adjEntry>> neighbours = getUnmarkedNeighbours(f, twin);
	switch (neighbours.size()) {
	case 0: {
		// f is leaf node of dual tree, simple scenario
		OGDF_ASSERT(embedding.rightFace(twin) == f);
		return CycleData(*graph, f, twin);
	}
	case 1: {
		// f has only one neighbour, so recurse on that
		auto neigh = neighbours.popFrontRet();
		OGDF_ASSERT(embedding.leftFace(neigh.second) == neigh.first);
		CycleData cycle = process(neigh.first, neigh.second);

		if (cycle.checkSize(graph->numberOfNodes())) {
			return cycle;
		}

		// check if the two nodes of the edge defined by adj are on the cycle
		if (!(cycle.isInCycle(adj->theNode()) && cycle.isInCycle(adj->twinNode()))) {
			cycle.addTriangle(twin);
		} else { // both of the nodes of the edge were on the cycle, so opportunity to make the cycle smaller
			cycle.removeTriangle(twin);
		}

		return cycle;
	}
	case 2: {
		// current face has two unprocessed neighbours, recurse on both, if one of them is big enough return, else join them
		auto neigh1 = neighbours.popFrontRet();
		OGDF_ASSERT(embedding.leftFace(neigh1.second) == neigh1.first);
		CycleData cycle1 = process(neigh1.first, neigh1.second);

		if (cycle1.checkSize(graph->numberOfNodes())) {
			return cycle1;
		}

		auto neigh2 = neighbours.popFrontRet();
		OGDF_ASSERT(embedding.leftFace(neigh2.second) == neigh2.first);
		CycleData cycle2 = process(neigh2.first, neigh2.second);

		if (cycle2.checkSize(graph->numberOfNodes())) {
			return cycle2;
		}

		// we need to join them
		return CycleData(*graph, cycle1, cycle2);
	}
	default:
		OGDF_ASSERT(false); // This should have been impossible, a triangle can not have more than three neighbours
		OGDF_THROW(ogdf::AlgorithmFailureException);
		return CycleData();
	}
}

// counts the unprocessed neighbours of a face
List<std::pair<face, adjEntry>> SeparatorDualHelper::getUnmarkedNeighbours(face f, adjEntry adj) {
	OGDF_ASSERT(embedding.rightFace(adj) == f);

	List<std::pair<face, adjEntry>> res;
	adjEntry nextAdj = adj->faceCycleSucc();

	while (nextAdj != adj) {
		if (tree->isInTree(nextAdj->theEdge())) {
			nextAdj = nextAdj->faceCycleSucc();
			continue;
		}

		face next = embedding.leftFace(nextAdj);
		if (!marked[next]) {
			res.pushBack(std::make_pair(next, nextAdj));
			marked[next] = true;
		}
		nextAdj = nextAdj->faceCycleSucc();
	}
	return res;
}

}

}
