/** \file
 *
 * \author Marcel Radermacher
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

#ifdef OGDF_INCLUDE_CGAL

#	include <ogdf/geometric/cr_min/datastructure/TimestampFlags.h>
#	include <ogdf/geometric/cr_min/graph/OGDFGraphWrapper.h>

#	include <iostream>
#	include <queue>

namespace ogdf {
namespace internal {
namespace gcm {
namespace graph {

template<typename Graph = OGDFGraphWrapper, typename Flags = datastructure::TimestampFlags>
class BFS {
private:
	using Node = typename Graph::Node;
	using Edge = typename Graph::Edge;

	const Graph& graph;
	std::queue<Node> queue;
	Flags visited;

	static void settle_nothing(Node) { /* nothing to do*/
	}

	static bool expand_all(Node, typename Graph::Edge) { return true; }

	static void traverse_nothing(Node, Edge) { /* nothing to do*/
	}

public:
	BFS(const Graph& graph_) : graph(graph_), visited(graph.max_node_index() + 1) {
		/*nothing to do*/
	}

	template<typename Range, typename FSettle, typename FExpandEdge, typename FTraverseEdge,
			typename InitVisit, typename AlreadyVisited, typename MarkVisited>
	void traverse2(Range& sources, InitVisit&& init_visit, AlreadyVisited&& already_visited,
			MarkVisited&& mark_visited, FSettle&& settle_node = settle_nothing,
			FExpandEdge&& expand_edge = expand_all, FTraverseEdge&& traverse_edge = traverse_nothing) {
		init_visit();
		for (auto v : sources) {
			queue.push(v);
			mark_visited(v);
		}
		unsigned int iterations = 0;
		while (!queue.empty()) {
			OGDF_ASSERT(iterations <= graph.number_of_nodes());
			++iterations;
			Node current_node = queue.front();
			queue.pop();
			settle_node(current_node);
			for (const Edge& edge : graph.edges(current_node)) {
				Node opposite = edge->opposite(current_node);
				if (expand_edge(current_node, edge) && !already_visited(opposite)) {
					traverse_edge(current_node, edge);
					queue.push(opposite);
					mark_visited(opposite);
				}
			}
		}
	}

	template<typename Range, typename FSettle, typename FExpandEdge, typename FTraverseEdge>
	void traverse(Range sources, FSettle&& settle_node = settle_nothing,
			FExpandEdge&& expand_edge = expand_all, FTraverseEdge&& traverse_edge = traverse_nothing) {
		traverse2(
				sources, [&]() { visited.clear(); },
				[&](Node u) { return visited.is_set(u->index()); },
				[&](Node u) { visited.set(u->index()); }, settle_node, expand_edge, traverse_edge);
	}

	template<typename Range, typename FSettle>
	void traverse(Range& source, FSettle&& settle_node) {
		traverse(source, settle_node, expand_all, traverse_nothing);
	}

	template<typename FSettle, typename FExpandEdge, typename FTraverseEdge>
	void traverse(Node& source, FSettle&& settle_node = settle_nothing,
			FExpandEdge&& expand_edge = expand_all, FTraverseEdge&& traverse_edge = traverse_nothing) {
		std::vector<Node> sources = {source};
		traverse(sources, settle_node, expand_edge, traverse_edge);
	}
};
}
}
}
}

#endif
