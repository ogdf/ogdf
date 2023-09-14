/** \file
 * \brief Implementation of class SeparatorLiptonTarjanFC.
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

#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/graphalg/SeparatorLiptonTarjanFC.h>

#include <stdexcept>

namespace ogdf {

void BFSTreeFC::construct() {
	List<List<node>> levels;

	SListPure<node> bfs; // current queue of BFS
	bfs.pushBack(root);

	int currentLevel = 0;

	while (!bfs.empty()) {
		// store all nodes of the current level
		List<node> next;
		for (const node no : bfs) {
			next.pushBack(no);
			levelOfNode[no] = currentLevel;
		}
		levels.pushBack(next);

		bfs.clear();

		// collect all nodes of the next level, record their parents
		for (const node w : next) {
			for (const adjEntry adj : w->adjEntries) {
				node v = adj->twinNode();
				if (!mark[v]) {
					mark[v] = true;
					bfs.pushBack(v);
					parentOfNode[v] = w;
					edgeToParent[v] = adj->twin();
					inTree[adj->theEdge()] = true;
				}
			}
		}
		currentLevel++;
	}

	// fill up descendant information from the bottom
	for (const auto& lvl : reverse(levels)) {
		for (const node x : lvl) {
			if (x != root) {
				descendantsOfNode[parentOfNode[x]] += descendantsOfNode[x];
				childrenOfNode[parentOfNode[x]].pushBack(x);
			}
		}
	}
}

void TriangulatingBFSTree::visit(node v, node parent, adjEntry adj, SListPure<node>& bfs) {
	if (!mark[v]) {
		mark[v] = true;
		bfs.pushBack(v);
		parentOfNode[v] = parent;
		edgeToParent[v] = adj->twin();
		inTree[adj->theEdge()] = true;
	}
}

void TriangulatingBFSTree::construct() {
	List<List<node>> levels;

	SListPure<node> bfs; // current queue of BFS
	bfs.pushBack(root);

	int currentLevel = 0;

	while (!bfs.empty()) {
		// store all nodes of the current level
		List<node> next;
		for (const node no : bfs) {
			next.pushBack(no);
			levelOfNode[no] = currentLevel;
		}
		levels.pushBack(next);

		bfs.clear();

		// collect all nodes of the next level, record their parents
		for (const node w : next) {
			std::set<node> neighbors; // storing neighbors of w for quick access
			List<adjEntry> oldAdjEntries;
			for (const adjEntry adj : w->adjEntries) {
				neighbors.insert(adj->twinNode());
				oldAdjEntries.pushBack(adj);
			}

			for (const adjEntry adj : oldAdjEntries) {
				node v = adj->twinNode();
				visit(v, w, adj, bfs);

				adjEntry adjBefore = adj;
				adjEntry nextInFace = adj->faceCycleSucc();

				while (nextInFace->twinNode() != w) {
					if (neighbors.find(nextInFace->twinNode()) == neighbors.end()) {
						// w is not connected to nextInFace->twinNode, but they are in the same face, so make connection
						nextInFace = nextInFace->faceCycleSucc();
						edge newEdge = pGraph->newEdge(adjBefore, nextInFace, Direction::after);
						neighbors.insert(nextInFace->theNode());
						adjBefore = newEdge->adjSource();
						visit(nextInFace->theNode(), w, newEdge->adjSource(), bfs);
					} else {
						break; // next one would end loop anyway
					}
				}
			}
		}
		currentLevel++;
	}

	// fill up descendant information from the bottom
	for (const auto& lvl : reverse(levels)) {
		for (const node x : lvl) {
			if (x != root) {
				descendantsOfNode[parentOfNode[x]] += descendantsOfNode[x];
				childrenOfNode[parentOfNode[x]].pushBack(x);
			}
		}
	}
}

bool SeparatorLiptonTarjanFC::doSeparate(const Graph& G, List<node>& separator, List<node>& first,
		List<node>& second) {
	if (useTriangulatingBFS) {
		tree = std::make_shared<TriangulatingBFSTree>(*graph, getStartNode(*graph));

		// if any of these asserts fail, I messed up tri-BFS
		OGDF_ASSERT(isPlanar(*graph));
		OGDF_ASSERT(isSimple(*graph));
		OGDF_ASSERT(graph->representsCombEmbedding());
		OGDF_ASSERT(graph->numberOfEdges() == 3 * graph->numberOfNodes() - 6);
	} else {
		// 0. embed and triangulate the planar embedding of the graph
		triangulate(*graph);

		// 1. construct tree
		tree = std::make_shared<BFSTreeFC>(*graph, getStartNode(*graph));
	}


	return findCycle(separator, first, second);
}

bool SeparatorLiptonTarjanFC::findCycle(List<node>& separator, List<node>& first, List<node>& second) {
	// 3. choose any non-tree edge as start-Edge for new cycle
	Cycle cycle(tree.get(), chooseEdge());

	// 4. iterate on cycle until sufficiently good cycle is found
	while (cycle.getInsideCost() > 2.0 / 3.0 * graph->nodes.size()) {
		cycle = cycle.expandCycle();
	}
	exitPoint = "cycle_expansion";
	cycle.fillLists(separator, first, second, true);

	return true;
}

edge SeparatorLiptonTarjanFC::chooseEdge() const {
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
