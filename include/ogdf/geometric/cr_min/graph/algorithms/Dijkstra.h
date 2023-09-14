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

#include <ogdf/basic/PriorityQueue.h>
#include <ogdf/geometric/cr_min/datastructure/TimestampFlags.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace graph {
template<typename Graph, typename Flags = datastructure::TimestampFlags>
class Dijkstra {
private:
	using Node = typename Graph::Node;
	using Edge = typename Graph::Edge;
	using Heap = PrioritizedQueue<Node, double>;
	using Element = typename Heap::Handle;

	const Graph& graph;
	Flags visited;
	Heap heap;
	std::vector<Element> reference;
	std::vector<double> distances;

public:
	static void settle_nothing(typename Graph::Node, double weight) { /* nothing to do*/
	}

	static bool expand_all(typename Graph::Node, typename Graph::Edge) { return true; }

	static void traverse_nothing(typename Graph::Node, typename Graph::Edge) { /* nothing to do*/
	}

	Node cycle_vertex = nullptr;

	Dijkstra(const Graph& _graph)
		: graph(_graph)
		, visited(graph.max_node_index() + 1)
		, heap()
		, reference(graph.max_node_index() + 1, nullptr)
		, distances(graph.max_node_index() + 1, 0) {
		// nothing to do
	}

	//return true on success, false on negative cycle
	template<typename Range, typename FWeight, typename FSettle, typename FExpand, typename FTraverse>
	bool traverse(const Range& sources, FWeight&& weight, FSettle&& settle, FExpand&& expand,
			FTraverse&& f_traverse) {
		visited.clear();
		heap.clear();

		for (auto v : sources) {
			reference[v->index()] = heap.push(v, 0);
			distances[v->index()] = 0;
			visited.set(v->index());
		}

		while (!heap.empty()) {
			const Node v = heap.topElement();
			const double current_weight = distances[v->index()];
			heap.pop();
			reference[v->index()] = nullptr;
			settle(v, current_weight);
			for (const Edge e : graph.edges(v)) {
				const Node w = e->opposite(v);
				if (expand(v, e)) {
					if (!visited.is_set(w->index())) {
						f_traverse(v, e);
						distances[w->index()] = current_weight + weight(v, e);
						reference[w->index()] = heap.push(w, current_weight + weight(v, e));
						visited.set(w->index());
					} else if (current_weight + weight(v, e) < distances[w->index()]) {
						distances[w->index()] = current_weight + weight(v, e);
						f_traverse(v, e);
						if (!reference[w->index()]) {
							cycle_vertex = w;
							return false;
						}

						heap.decrease(reference[w->index()], distances[w->index()]);
					}
					auto r = reference[w->index()];
				}
			}
		}
		return true;
	}

	template<typename FWeight, typename FSettle, typename FExpand, typename FTraverse>
	bool traverse_single(const Node& source, FWeight&& weight, FSettle&& settle, FExpand&& expand,
			FTraverse&& f_traverse) {
		return traverse(std::vector<Node>(1, source), weight, settle, expand, f_traverse);
	}
};

}
}
}
}
