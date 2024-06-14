#pragma once

#include <ogdf/basic/Graph.h>

#include "Bijection.h"
#include "Iterators.h"

using namespace ogdf;

void moveEnd(Graph& G, edge e, node keep_end, node new_end);

void moveEnd(Graph& G, adjEntry keep_adj, adjEntry new_adj, Direction dir = Direction::after);

edge splitEdge(Graph& G, edge old_edge, node new_adj_to_source, node new_adj_to_target,
		int new_edge_idx = -1);

adjEntry splitEdge(Graph& G, adjEntry adj, node new_adj_to_node, node new_adj_to_twin,
		int new_edge_idx = -1);

bool joinEdge(Graph& G, edge u_e, edge v_e, node u, node v);

bool joinEdge(Graph& G, adjEntry u_adj, adjEntry v_adj, node u, node v);

std::pair<node, node> split(Graph& G, PipeBij& bij, const EdgeArray<int>* split_idcs = nullptr,
		const EdgeArray<bool>* split_reverse = nullptr, int src_idx = -1, int tgt_idx = -1);

void join(Graph& G, node u, node v, PipeBij& bij, List<bool>* reverse_v = nullptr);

void assertStarCentreAndRay(node centre, node ray);

node getCentreOfStar(node g_n);

enum OrderComp { SAME, REVERSED, DIFFERENT };

OrderComp compareCyclicOrder(node n, List<adjEntry>& o, bool full_check = false);

void moveAdjToFront(Graph& G, adjEntry f);

void moveAdjToBack(Graph& G, adjEntry b);
