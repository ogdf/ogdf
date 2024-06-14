#pragma once

#include <ogdf/basic/Graph.h>
#include <ogdf/decomposition/StaticPlanarSPQRTree.h>

#include "PQPlanarityComponents.h"
#include "utils/RegisteredMultiArray.h"
#include "utils/TwoSAT.h"

using namespace ogdf;

struct BlockEmbedding;

using GnMultiArray = RegisteredMultiArray<node, BlockEmbedding*, node, NodeArray>;

struct BlockEmbedding {
	Graph subgraph;
	StaticPlanarSPQRTree* spqr = nullptr;
	List<node> q_vertices;
	NodeArray<twosat_var> rigid_vars;

	GnMultiArray& Gn_to_subgraph;
	EdgeArray<edge> subgraph_to_Ge;

	explicit BlockEmbedding(GnMultiArray& gnToSubgraph) : Gn_to_subgraph(gnToSubgraph) { }

	virtual ~BlockEmbedding() { delete spqr; }

	void init(Graph& G, PQPlanarityComponents& components, node bc, EdgeArray<edge>& Ge_to_subgraph,
			EdgeArray<BlockEmbedding*>& Ge_to_block);

	bool addQVertex(node q, EdgeArray<edge>& Ge_to_subgraph, TwoSAT& sat, twosat_var part_var);
};