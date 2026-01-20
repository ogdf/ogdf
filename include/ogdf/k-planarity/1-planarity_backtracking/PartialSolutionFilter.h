/** \file
 * \brief Filters for partial solutions in backtracking for the 1-Planarity problem.
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
#include <ogdf/k-planarity/1-planarity_backtracking/OnePlanarization.h>

namespace ogdf::oneplan_backtracking {

//! Interface of filters for partial solutions for 1-Planarity
/**
 *	A filter takes as input a partial solution representing a node in the backtracking tree and
 *	rejects if it detects that the partial solution cannot be extended to a solution.
 */
class OGDF_EXPORT PartialSolutionFilter {
public:
	//! Returns true if the partial solution represented by \p pl can be classified as non-realizable, false otherwise.
	virtual bool canCut(OnePlanarization& pl) = 0;

	virtual ~PartialSolutionFilter() = default;
};

//! Rejects partial solutions for 1-Planarity if the crossable edges violate the density for 1-Planarity.
class OGDF_EXPORT FreeEdgesDensityFilter : public PartialSolutionFilter {
public:
	//! Returns true if the subgraph of \p induced by its free edges violates the density for 1-Planarity, false otherwise.
	bool canCut(OnePlanarization& pl) override {
		ogdf::NodeSet seenVertices(pl);
		for (ogdf::edge e : pl.freeEdges()) {
			seenVertices.insert(e->target());
			seenVertices.insert(e->source());
		}

		return pl.freeEdges().size() > 4 * seenVertices.size() - 8;
	}
};

//! Rejects partial solutions for 1-Planarity if they are too dense.
class OGDF_EXPORT PlanarizationDensityFilter : public PartialSolutionFilter {
public:
	//! Returns true if the subgraph of \p pl induced by its free edges violates the density formula for 1-Planarity, false otherwise.
	bool canCut(OnePlanarization& pl) override {
		return pl.numberOfEdges() > 4 * pl.numberOfNodes() - 8;
	}
};

//! Rejects partial solutions for 1-Planarity if the uncrossable edges induce a non-planar graph.
class OGDF_EXPORT SaturatedSubgraphPlanarityFilter : public PartialSolutionFilter {
	//! Returns true if the subgraph of \p pl induced by its uncrossable edges is non-planar, false otherwise.
	bool canCut(OnePlanarization& pl) override {
		ogdf::GraphCopy saturatedSubgraph;
		pl.getSaturatedSubgraph(saturatedSubgraph);
		return !isPlanar(saturatedSubgraph);
	}
};

//! Rejects partial solutions for 1-Planarity if a separating cycle with too few crossable edges is found.
/**
 * For every pair (e, f) of crossed edges in a given partial solution, this filter searches a cycle c
 * containing e and tests whether the maximum number of edge-disjoint paths between the endpoints of f
 * that are vertex-disjoint from c exceeds the number of crossable edges contained in c.
 * If this is the case, the partial solution is not realizable and is rejected.
 */
class OGDF_EXPORT SeparatingCycleFilter : public PartialSolutionFilter {
public:
	//! Returns true if \p pl can be rejected due to a found separating cycle, false otherwise.
	bool canCut(OnePlanarization& pl) override;

private:
	bool test(OnePlanarization& pl, ogdf::node crossingVertex, ogdf::NodePair cyclePair,
			ogdf::NodePair otherPair);

	//! Searches for a path between \p cyclePair that is disjoint from \p otherPair.
	static int minFreeEdgeCycle(const OnePlanarization& pl, node crossingVertex,
			NodePair& cyclePair, NodePair& otherPair, EdgeSet& cycleEdges, NodeSet& cycleNodes);

	//! Contracts all edges of \p cp that cannot cross any edge of the cycle.
	static void contractUncrossableEdges(const OnePlanarization& pl, GraphCopy& cp,
			const EdgeSet& cycleEdges, const NodeSet& cycleNodes, node& v1, node& v2);
};
}
