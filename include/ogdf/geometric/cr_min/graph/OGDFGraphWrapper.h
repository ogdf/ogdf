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

#	include <ogdf/basic/Graph.h>
#	include <ogdf/fileformats/GraphIO.h>
#	include <ogdf/geometric/cr_min/datastructure/Iterators.h>
#	include <ogdf/geometric/cr_min/datastructure/OGDFVector.h>
#	include <ogdf/geometric/cr_min/graph/ogdf_iterator.h>
#	include <ogdf/geometric/cr_min/tools/ogdf/Universal.h>

#	include <algorithm>
#	include <iostream>
#	include <map>
#	include <vector>

namespace ogdf {
namespace internal {
namespace gcm {
namespace graph {

/*! a warpper for ogdf graphs
 */
class OGDFGraphWrapper {
public:
	using Node = node;
	using Edge = edge;
	using AdjEntry = adjEntry;

protected:
	std::shared_ptr<Graph> self_generated_graph;
	const Graph& graph;

	static void id_sg_handle_node(Node, Node) { }

	static void id_sg_handle_edge(Edge, Edge) { }


public:
	OGDFGraphWrapper() : self_generated_graph(new Graph()), graph(*self_generated_graph) {
		/*nothing to do*/
	}

	OGDFGraphWrapper(const Graph& _graph) : self_generated_graph(nullptr), graph(_graph) {
		/*nothing to do*/
	}

	~OGDFGraphWrapper() { }

	inline const Graph& get_graph() { return graph; }

	inline const Graph& get_ogdf_graph() { return graph; }

	inline size_t number_of_nodes() const { return graph.numberOfNodes(); }

	inline size_t number_of_edges() const { return graph.numberOfEdges(); }

	inline size_t max_node_index() const { return graph.maxNodeIndex(); }

	inline size_t max_edge_index() const { return graph.maxEdgeIndex(); }

	inline Node get_node(unsigned int i) const {
		for (auto w : nodes()) {
			if (w->index() == (int)i) {
				return w;
			}
		}
		return nullptr;
	}

	/**
	 * Look up the the  vw
	 * \param v a node
	 * \param w a node
	 * \return the edge e=(v,w)  or e = (w, v) if it exists, otherwise null
	 */
	inline Edge search_edge(Node v, Node w) const {
		if (v->degree() < w->degree()) {
			for (auto e : edges(v)) {
				if (e->opposite(v) == w) {
					return e;
				}
			}
			return nullptr;
		} else {
			for (auto e : edges(w)) {
				if (e->opposite(w) == v) {
					return e;
				}
			}
			return nullptr;
		}
	}

	/** checks whether there is an edge vw or wv
	 * \param v a node of the graph
	 * \param w a node of the graph
	 * \return true if and only if there is an edge vw or wv
	 */
	inline bool are_adjacent(Node v, Node w) const { return search_edge(v, w) != nullptr; }

	/**
	 * Range over all nodes of the graph
	 */
	inline datastructure::IteratorRange<NodeIterator> nodes() const {
		if (number_of_nodes() > 0) {
			return datastructure::IteratorRange<NodeIterator>(NodeIterator(graph.firstNode()),
					NodeIterator(graph.lastNode()->succ()));
		} else {
			return datastructure::IteratorRange<NodeIterator>(Node(0), Node(0));
		}
	}

	/**
	 * Range over all edges of the graph
	 */
	inline datastructure::IteratorRange<EdgeIterator> edges() const {
		if (number_of_edges() > 0) {
			return datastructure::IteratorRange<EdgeIterator>(EdgeIterator(graph.firstEdge()),
					EdgeIterator(graph.lastEdge()->succ()));
		} else {
			return datastructure::IteratorRange<EdgeIterator>(Edge(0), Edge(0));
		}
	}

	/** Range of all edges incident to w.
	 * \param w a node
	 * \return range of edges incident to w
	 */
	inline datastructure::IteratorRange<IncidentEdgeIterator> edges(const Node w) const {
		if (w->degree() == 0) {
			return datastructure::IteratorRange<IncidentEdgeIterator>(
					IncidentEdgeIterator(adjEntry(0)), IncidentEdgeIterator(adjEntry(0)));
		} else {
			return datastructure::IteratorRange<IncidentEdgeIterator>(
					IncidentEdgeIterator(w->firstAdj()), IncidentEdgeIterator(w->lastAdj()->succ()));
		}
	}

	/**
	 * Range of adjacency entries of a node w
	 * \param w a node
	 * \return range of adjacency entries
	 */
	inline datastructure::IteratorRange<AdjEntryIterator> adj(const Node w) const {
		if (w->degree() == 0) {
			return datastructure::IteratorRange<AdjEntryIterator>(AdjEntryIterator(adjEntry(0)),
					AdjEntryIterator(adjEntry(0)));
		} else {
			return datastructure::IteratorRange<AdjEntryIterator>(AdjEntryIterator(w->firstAdj()),
					AdjEntryIterator(w->lastAdj()->succ()));
		}
	}

	/** Range of neighboring nodes of a node w. A neighbor can occur twice.
	 * \param w
	 * \return range of neighbors of w
	 */
	inline datastructure::IteratorRange<AdjacentNodeIterator> neighbors(const Node w) const {
		if (w->degree() == 0) {
			return datastructure::IteratorRange<AdjacentNodeIterator>(
					AdjacentNodeIterator(adjEntry(0)), AdjacentNodeIterator(adjEntry(0)));

		} else {
			return datastructure::IteratorRange<AdjacentNodeIterator>(
					AdjacentNodeIterator(w->firstAdj()), AdjacentNodeIterator(w->lastAdj()->succ()));
		}
	}
};

}
}
}
}


#endif
