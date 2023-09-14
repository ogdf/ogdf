/** \file
 * \brief Declaration of class SeparatorDualHelper.
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

#pragma once

#include <ogdf/basic/FaceArray.h>
#include <ogdf/graphalg/PlanarSeparatorModule.h>

#include <unordered_set>

namespace ogdf {

namespace planar_separator {

//! Helper class for SeparatorDual and SeparatorDualFC.
/**
 * @ingroup ga-plansep
 */
class OGDF_EXPORT SeparatorDualHelper {
public:
	SeparatorDualHelper(std::shared_ptr<GraphCopy> pGraph, std::shared_ptr<BFSTree> pTree)
		: graph {pGraph}, tree {pTree} { }

	/**
	 * Auxiliary lightweight data structure to represent cycles.
	 */
	struct CycleData {
		unsigned int inside; // contains number of nodes on the inside

		/** contains nodes on the cycle, starting with v, ending with u, where e = (u,v) is the initial edge */
		List<node> cycle;

		std::unordered_set<node> nodes;

		/**
		 * Checks if a node lies on the cycle.
		 *
		 * @param n the node to be checked
		 * @return true if node is on the cycle
		 */
		bool isInCycle(node n) { return nodes.find(n) != nodes.end(); }

		/* Methods to add and remove nodes from the cycle. */

		void pushBack(node n) {
			OGDF_ASSERT(!isInCycle(n));
			cycle.pushBack(n);
			nodes.insert(n);
		}

		void pushFront(node n) {
			OGDF_ASSERT(!isInCycle(n));
			cycle.pushFront(n);
			nodes.insert(n);
		}

		void popBack() {
			node x = cycle.popBackRet();
			nodes.erase(x);
		}

		void popFront() {
			node x = cycle.popFrontRet();
			nodes.erase(x);
		}

		/**
		 * Empty Constructor - only used to be able to throw empty CD in case of algorithm failure.
		 * Do not use this anywhere.
		 */
		CycleData() { }

		// case 1
		/**
		 * Constructor. Constructs a CycleData from a single face.
		 *
		 * @param G the graph
		 * @param f the face
		 * @param adj adjEntry via which we entered the face
		 */
		CycleData(const Graph& G, const face f, const adjEntry adj) : inside {0} {
			pushBack(adj->twinNode());
			pushBack(adj->faceCycleSucc()->twinNode());
			pushBack(adj->theNode());
		}

		// case 4
		/**
		 * Constructor. Constructs a CycleData by merging two Cycles.
		 *
		 * @param G the graph
		 * @param first the first cycle
		 * @param second the second cycle
		 */
		CycleData(const Graph& G, CycleData& first, CycleData& second) {
			// find path
			List<node> path;

			auto firstIt = first.cycle.crbegin();
			auto secondIt = second.cycle.cbegin();

			OGDF_ASSERT(*firstIt == *secondIt);

			while (firstIt != first.cycle.crend() && secondIt != second.cycle.cend()
					&& *firstIt == *secondIt) {
				path.pushBack(*firstIt);
				++firstIt;
				++secondIt;
			}
			node x = path.back();

			// number of nodes on the inside
			inside = first.inside + second.inside + path.size() - 1;

			// join cycles
			for (auto it = first.cycle.cbegin(); *it != x; ++it) {
				pushBack(*it);
			}

			auto it = second.cycle.cbegin();
			while (*it != x) {
				++it;
			}
			for (; it != second.cycle.cend(); ++it) {
				pushBack(*it);
			}
		}

		// case 2
		/**
		 * Expands the cycle by adding a triangle.
		 *
		 * @param adj the adjEntry of the triangle that is added
		 */
		void addTriangle(adjEntry adj) {
			OGDF_ASSERT(!(isInCycle(adj->theNode()) && isInCycle(adj->twinNode())));
			if (isInCycle(adj->theNode())) {
				pushFront(adj->twinNode());
			} else {
				pushBack(adj->theNode());
			}
		}

		// case 3
		/**
		 * Expands the cycle by removing a triangle.
		 *
		 * @param adj the adjEntry of the triangle that is removed
		 */
		void removeTriangle(adjEntry adj) {
			// two different scenarios, depending on which adjacent triangle contains the cycle
			OGDF_ASSERT(isInCycle(adj->theNode()) && isInCycle(adj->twinNode()));

			if (cycle.back() == adj->theNode()) {
				popFront();
			} else if (cycle.front() == adj->twinNode()) {
				popBack();
			} else {
				OGDF_ASSERT(false); // Removing triangle was impossible because graph layout was illegitimate.
			}

			inside++;
		}

		/**
		 * Checks the size of this cycle.
		 *
		 * @param n the size of the graph to be checked against
		 * @return true if this cycle is large enough
		 */
		bool checkSize(int n) const { return inside + cycle.size() > 1.0 / 3.0 * n; }
	};

	std::shared_ptr<GraphCopy> graph;
	std::shared_ptr<BFSTree> tree;

	// components needed for DFS
	CombinatorialEmbedding embedding;
	FaceArray<bool> marked;

	/**
	 * Recursive Depth First Search over the faces of the dual of the graph.
	 *
	 * @return a data structure that holds information on the cycle
	 */
	CycleData dfs();

	/**
	 * Processes a face: One step of the recursion in the DFS.
	 *
	 * @param f the face to be processed
	 * @param adj the adjEntry via which the face was entered
	 * @return the next CycleData
	 */
	CycleData process(face f, adjEntry adj);

	/**
	 * Finds all unmarked neighbours of a face.
	 *
	 * @param f the face to be examined
	 * @param adj the adjEntry via which the face was entered
	 * @return a list containing the unmarked neighbouring faces and the entries that lead to them
	 */
	List<std::pair<face, adjEntry>> getUnmarkedNeighbours(face f, adjEntry adj);
};

}

}
