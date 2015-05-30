/** \file
 * \brief Implements the class NonPlanarCore.
 *
 * \author Carsten Gutwenger
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.txt in the root directory of the OGDF installation for details.
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
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/


#include <ogdf/planarity/NonPlanarCore.h>
#include <ogdf/decomposition/StaticSPQRTree.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/Queue.h>
#include <ogdf/basic/CombinatorialEmbedding.h>
#include <ogdf/basic/FaceArray.h>
#include <ogdf/basic/simple_graph_alg.h>


namespace ogdf {


NonPlanarCore::NonPlanarCore(const Graph &G) : m_pOriginal(&G), m_orig(m_graph),
	m_real(m_graph,nullptr), m_mincut(m_graph), m_cost(m_graph)
{
	if(G.numberOfNodes() <= 4)
		return; // nothing to do; planar graph => empty core

	OGDF_ASSERT(isBiconnected(G));

	// Build SPQR-tree of graph
	StaticSPQRTree T(G);

	// mark tree nodes in the core
	NodeArray<bool> mark;
	markCore(T,mark);

	NodeArray<node> map(G,nullptr);
	NodeArray<node> mapAux(G,nullptr);
	const Graph &tree = T.tree();

	for (node v : tree.nodes) {
		if(mark[v] == false)
			continue;

		Skeleton &S = T.skeleton(v);
		for (edge e : S.getGraph().edges) {
			node src = S.original(e->source());
			node tgt = S.original(e->target());

			if(map[src] == nullptr) {
				m_orig[map[src] = m_graph.newNode()] = S.original(e->source());
			}
			if(map[tgt] == nullptr) {
				m_orig[map[tgt] = m_graph.newNode()] = S.original(e->target());
			}

			if(S.isVirtual(e)) {
				node w = S.twinTreeNode(e);
				if(mark[w] == false) {
					// new virtual edge in core graph
					edge eC = m_graph.newEdge(map[src],map[tgt]);
					traversingPath(S,e,m_mincut[eC],mapAux);
				}

			} else {
				// new real edge in core graph
				edge eC = m_graph.newEdge(map[src],map[tgt]);
				m_real[eC] = S.realEdge(e);
				m_mincut[eC].pushBack(S.realEdge(e));
			}
		}
	}

	for (edge e : m_graph.edges) {
		m_cost[e] = m_mincut[e].size();
	}
}


// This function marks all nodes in the SPQR-tree which induce the
// non-planar core.
void NonPlanarCore::markCore(const SPQRTree &T, NodeArray<bool> &mark)
{
	const Graph &tree = T.tree();

	// We mark every tree node that belongs to the core
	mark.init(tree,true);  // start with all nodes and unmark planar leaves
	NodeArray<int> degree(tree);

	Queue<node> Q;

	for (node v : tree.nodes) {
		degree[v] = v->degree();
		if(degree[v] <= 1) // also append deg-0 node (T has only one node)
			Q.append(v);
	}

	while(!Q.empty())
	{
		node v = Q.pop();

		// if v has a planar skeleton
		if(T.typeOf(v) != SPQRTree::RNode ||
			isPlanar(T.skeleton(v).getGraph()) == true)
		{
			mark[v] = false; // unmark this leaf

			node w = nullptr;
			for (adjEntry adj : v->adjEdges) {
				node x = adj->twinNode();
				if(mark[x] == true) {
					w = x; break;
				}
			}

			if(w != nullptr) {
				--degree[w];
				if(degree[w] == 1)
					Q.append(w);
			}
		}
	}
}

struct OGDF_EXPORT QueueEntry
{
	QueueEntry(node p, node v) : m_parent(p), m_current(v) { }

	node m_parent;
	node m_current;
};

void NonPlanarCore::traversingPath(Skeleton &Sv, edge eS, List<edge> &path, NodeArray<node> &mapV)
{
	const SPQRTree &T = Sv.owner();

	//-----------------------------------------------------
	// Build the graph representing the planar st-component
	Graph H;
	EdgeArray<edge> mapE(H,nullptr);
	SListPure<node> nodes;

	Queue<QueueEntry> Q;
	Q.append(QueueEntry(Sv.treeNode(),Sv.twinTreeNode(eS)));

	while(!Q.empty())
	{
		QueueEntry x = Q.pop();
		node parent = x.m_parent;
		node current = x.m_current;

		const Skeleton &S = T.skeleton(current);

		for (edge e : S.getGraph().edges) {
			if(S.isVirtual(e) == true)
				continue;

			node src = S.original(e->source());
			node tgt = S.original(e->target());

			if(mapV[src] == nullptr) {
				nodes.pushBack(src);
				mapV[src] = H.newNode();
			}
			if(mapV[tgt] == nullptr) {
				nodes.pushBack(tgt);
				mapV[tgt] = H.newNode();
			}

			mapE[H.newEdge(mapV[src],mapV[tgt])] = S.realEdge(e);
		}

		for (adjEntry adj : current->adjEdges) {
			node w = adj->twinNode();
			if(w != parent)
				Q.append(QueueEntry(current,w));
		}
	}

	// add st-edge
	edge e_st = H.newEdge(mapV[Sv.original(eS->source())],mapV[Sv.original(eS->target())]);

	// Compute planar embedding of H
#ifdef OGDF_DEBUG
	bool ok =
#endif
		planarEmbed(H);
	OGDF_ASSERT(ok)
	CombinatorialEmbedding E(H);

	//---------------------------------
	// Compute corresponding dual graph
	Graph dual;
	FaceArray<node> nodeOf(E);
	EdgeArray<adjEntry> primalAdj(dual);

	// insert a node in the dual graph for each face in E
	for (face f : E.faces)
		nodeOf[f] = dual.newNode();


	node s = nodeOf[E.rightFace(e_st->adjSource())];
	node t = nodeOf[E.rightFace(e_st->adjTarget())];

	// Insert an edge into the dual graph for each adjacency entry in E.
	// The edges are directed from the left face to the right face.
	for (node v : H.nodes)
	{
		for (adjEntry adj : v->adjEdges)
		{
			// do not insert edges crossing e_st
			if(adj->theEdge() == e_st)
				continue;

			node vLeft  = nodeOf[E.leftFace (adj)];
			node vRight = nodeOf[E.rightFace(adj)];

			primalAdj[dual.newEdge(vLeft,vRight)] = adj;
		}
	}

	//---------------------------
	// Find shortest path in dual
	NodeArray<edge> spPred(dual,nullptr);
	QueuePure<edge> queue;

	edge eDual;
	forall_adj_edges(eDual,s) {
		if(s == eDual->source())
			queue.append(eDual);
	}

	// actual search (using bfs on directed dual)
	for( ; ; )
	{
		// next candidate edge
		edge eCand = queue.pop();
		node v = eCand->target();

		// leads to an unvisited node?
		if (spPred[v] == nullptr)
		{
			// yes, then we set v's predecessor in search tree
			spPred[v] = eCand;

			// have we reached t ...
			if (v == t)
			{
				// ... then search is done.
				// constructed list of used edges (translated to crossed
				// edges entries in G) from t back to s (including first
				// and last!)

				do {
					edge eDual = spPred[v];
					edge eG = mapE[primalAdj[eDual]->theEdge()];
					OGDF_ASSERT(eG != 0)
					path.pushFront(eG);
					v = eDual->source();
				} while(v != s);

				break;
			}

			// append next candidate edges to queue
			// (all edges leaving v)
			edge e;
			forall_adj_edges(e,v) {
				if (v == e->source())
					queue.append(e);
			}
		}
	}


	//---------
	// Clean-up
	for(node v : nodes)
		mapV[v] = nullptr;
}


} // end namespace ogdf
