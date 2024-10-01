/** \file
 * \brief Internal class used to embed a biconnected component with Q-vertices.
 *
 * \author Simon D. Fink <ogdf@niko.fink.bayern>
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
#include <ogdf/basic/List.h>
#include <ogdf/cluster/sync_plan/basic/RegisteredMultiArray.h>
#include <ogdf/cluster/sync_plan/basic/TwoSAT.h>
#include <ogdf/decomposition/StaticPlanarSPQRTree.h>

namespace ogdf::sync_plan {

class SyncPlanComponents;

namespace internal {
struct BlockEmbedding;

template<typename V>
using NA = NodeArray<V>;

using GnMultiArray = RegisteredMultiArray<node, BlockEmbedding*, node, NA>;

//! Internal class used to embed a biconnected component with Q-vertices.
struct BlockEmbedding {
	Graph subgraph;
	StaticPlanarSPQRTree* spqr = nullptr;
	List<node> q_vertices;
	NodeArray<twosat_var> rigid_vars;

	GnMultiArray& Gn_to_subgraph;
	EdgeArray<edge> subgraph_to_Ge;

	explicit BlockEmbedding(GnMultiArray& gnToSubgraph)
		: rigid_vars(nullptr, TwoSAT_Var_Undefined), Gn_to_subgraph(gnToSubgraph) { }

	virtual ~BlockEmbedding() { delete spqr; }

	void init(Graph& G, SyncPlanComponents& components, node bc, EdgeArray<edge>& Ge_to_subgraph,
			EdgeArray<BlockEmbedding*>& Ge_to_block);

	bool addQVertex(node q, EdgeArray<edge>& Ge_to_subgraph, TwoSAT& sat, twosat_var part_var);
};

}
}
