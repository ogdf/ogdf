/** \file
 * \brief Declaration of the OnePlanarization class.
 *
 * \author Matthias Pfretzschner
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
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/GraphSets.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/extended_graph_alg.h>

namespace ogdf::oneplan_backtracking {
class EdgePairPartition;

//! The graph corresponding to a partial solution for the 1-Planarity problem.
/**
 * Represents the graph corresponding to a \ref EdgePairPartition using a \ref GraphCopy of the
 * original input graph. Crossed edge pairs are replaced with crossing vertices and kite edges
 * around crossing vertices are inserted. If this graph is planar, it represents a solution for
 * the 1-Planarity problem.
 */
class OGDF_EXPORT OnePlanarization : public GraphCopy {
private:
	const EdgePairPartition* m_epp = nullptr;
	NodeSet m_saturatedVertices; // vertices that are incident to an uncrossable edge
	NodeSet m_crossingVertices;
	EdgeSet m_kiteEdges; // uncrossable
	EdgeSet m_crossingEdges; // edges incident to crossing vertices (uncrossable)
	EdgeSet m_freeEdges; // edges that still have at least one other edge they are allowed to cross
	EdgeSet m_remainingEdges; // uncrossable edges that are not kite edges or crossing edges

public:
	OnePlanarization() = default;

	//! Associates the graph with \p ep and computes the planarization
	explicit OnePlanarization(const EdgePairPartition* ep) { init(ep); }

	//! Re-initializes the graph with \p ep.
	void init(const EdgePairPartition* ep);

	//! Returns the associated \ref EdgePairPartition.
	const EdgePairPartition* getEdgePairPartition() const { return m_epp; }

	//! Returns true if the graph is planar and thus represents a solution.
	bool isPlanar() const { return ogdf::isPlanar(*this); }

	//! Returns the set of all crossing vertices in the graph.
	const NodeSet& crossingVertices() const { return m_crossingVertices; }

	//! Returns the set of all kite edges in the graph.
	/**
	 *	Kite edges are non-edges of the original graph between the endpoints of crossed edge pairs.
	 *	Note that these edges are not associated with edges of the original graph.
	 */
	const EdgeSet& kiteEdges() const { return m_kiteEdges; }

	//! Returns the set of all edges that are incident to a crossing vertex.
	/**
	 * The four edges incident to each crossing vertex are associated with the two corresponding
	 * edges of the original graph.
	 */
	const EdgeSet& crossingEdges() const { return m_crossingEdges; }

	//! Returns the set of all edges that may still cross some other edge.
	const EdgeSet& freeEdges() const { return m_freeEdges; }

	//! Returns the set of all uncrossable edges that are neither crossing edges nor kite edges.
	const EdgeSet& remainingEdges() const { return m_remainingEdges; }

	//! Computes the subgraph induced by uncrossable edges.
	void getSaturatedSubgraph(GraphCopy& cp) const {
		cp.init(*reinterpret_cast<const Graph*>(this));
		for (edge e : m_freeEdges) {
			cp.delEdge(cp.copy(e));
		}
	}

private:
	using GraphCopy::init;
};
}
