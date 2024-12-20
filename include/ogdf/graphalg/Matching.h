/** \file
 * \brief Declares simple matching functions
 *
 * \author Stephan Beyer, Thomas Klein
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
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/basic.h>

#include <array>

namespace ogdf {
template<class E, class INDEX>
class ArrayBuffer;
template<class E>
class List;

//! Simple algorithms for matchings
namespace Matching {

//! Checks in time O(|V| + size of \p matching) if the given set of edges represents a matching.
//! @ingroup ga-matching
template<typename EdgeContainer>
inline bool isMatching(const Graph& graph, const EdgeContainer& matching) {
	NodeArray<bool> covered {graph, false};

	for (edge e : matching) {
		for (node v : e->nodes()) {
			if (covered[v]) {
				return false;
			}
			covered[v] = true;
			if (e->isSelfLoop()) {
				break;
			}
		}
	}

	return true;
}

//! @copydoc isMaximal(const Graph&, const EdgeContainer&)
//! Sets \p addable to such an edge or \c nullptr if none is found.
template<typename EdgeContainer>
bool isMaximal(const Graph& graph, const EdgeContainer& matching, edge& addable) {
	addable = nullptr;

	EdgeArray<bool> covered {graph, false};

	for (edge e : matching) {
		for (node v : e->nodes()) {
			for (adjEntry adj : v->adjEntries) {
				covered[adj->theEdge()] = true;
			}
		}
	}

	for (edge e : graph.edges) {
		if (!covered[e]) {
			addable = e;
			return false;
		}
	}

	return true;
}

//! Checks in time O(|E|) if there are edges that could be added to \p matching.
//! @ingroup ga-matching
template<typename EdgeContainer>
inline bool isMaximal(const Graph& graph, const EdgeContainer& matching) {
	edge ignored;
	return isMaximal(graph, matching, ignored);
}

//! Checks in O(|V| + |E|) time if \p matching is a maximal matching.
//! @ingroup ga-matching
template<typename EdgeContainer>
inline bool isMaximalMatching(const Graph& graph, const EdgeContainer& matching) {
	return isMatching(graph, matching) && isMaximal(graph, matching);
}

//! Checks in O(1) if \p matching (assuming it is a matching and the graph is simple and connected) is perfect.
//! @ingroup ga-matching
template<typename EdgeContainer>
inline bool isPerfect(const Graph& graph, const EdgeContainer& matching) {
	return 2 * int(matching.size()) == graph.numberOfNodes();
}

//! Checks in O(|V| + size of \p matching) if \p matching is a perfect matching.
//! @ingroup ga-matching
template<typename EdgeContainer>
inline bool isPerfectMatching(const Graph& graph, const EdgeContainer& matching) {
	return isMatching(graph, matching) && isPerfect(graph, matching);
}

//! Obtains a maximal matching in O(|E|) time
//! @ingroup ga-matching
OGDF_EXPORT void findMaximalMatching(const Graph& graph, ArrayBuffer<edge>& matching);

/**
 * Finds a maximum cardinality matching in the bipartite graph \p G = (U+V, E)
 * in O(sqrt(|U+V|) * |E|) time by using the Hopcroft-Karp-Karzanov algorithm.
 *
 * @param G the bipartite graph
 * @param U all nodes in the first half of \p G
 * @param V all nodes in the second half of \p G
 * @param matching will hold the matching
 * @return the cardinality of the matching
 */
//! @ingroup ga-matching
OGDF_EXPORT int findMaximumCardinalityMatching(const Graph& G, const List<node>& U,
		const List<node>& V, EdgeArray<bool>& matching);

}
}
