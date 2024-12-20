/** \file
 * \brief Computes a connected subgraph G' of G containing node n.
 *
 * \author Thorsten Kerkhof
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

namespace ogdf {
namespace embedder {

template<class T>
class ConnectedSubgraph {
public:
	//constructor
	ConnectedSubgraph() { }

	/**
	 * \brief Computes a connected subgraph SG of G containing node nG.
	 * \param G is the original graph.
	 * \param SG is assigned the connected subgraph containing nG.
	 * \param nG is a node in G.
	 * \param nSG is assigned the corresponding node of nG in SG.
	 * \param nodeLengthG stores for each node of G its length.
	 * \param nodeLengthSG is assigned for each node of SG its length.
	 */
	static void call(const Graph& G, Graph& SG, const node& nG, node& nSG,
			const NodeArray<T>& nodeLengthG, NodeArray<T>& nodeLengthSG) {
		EdgeArray<T> edgeLengthG(G, 1);
		EdgeArray<T> edgeLengthSG;
		call(G, SG, nG, nSG, nodeLengthG, nodeLengthSG, edgeLengthG, edgeLengthSG);
	}

	/**
	 * \brief Computes a connected subgraph SG of G containing node nG.
	 * \param G is the original graph.
	 * \param SG is assigned the connected subgraph containing nG.
	 * \param nG is a node in G.
	 * \param nodeLengthG stores for each node of G its length.
	 * \param nodeLengthSG is assigned for each node of SG its length.
	 * \param nG_to_nSG is assigned a mapping of nodes in G to nodes in SG.
	 */
	static void call(const Graph& G, Graph& SG, const node& nG, const NodeArray<T>& nodeLengthG,
			NodeArray<T>& nodeLengthSG, NodeArray<node>& nG_to_nSG) {
		node nSG;
		NodeArray<node> nSG_to_nG;
		EdgeArray<edge> eSG_to_eG;
		EdgeArray<edge> eG_to_eSG;
		EdgeArray<T> edgeLengthG(G, 1);
		EdgeArray<T> edgeLengthSG;
		call(G, SG, nG, nSG, nSG_to_nG, eSG_to_eG, nG_to_nSG, eG_to_eSG, nodeLengthG, nodeLengthSG,
				edgeLengthG, edgeLengthSG);
	}

	/**
	 * \brief Computes a connected subgraph SG of G containing node nG.
	 * \param G is the original graph.
	 * \param SG is assigned the connected subgraph containing nG.
	 * \param nG is a node in G.
	 * \param nSG is assigned the corresponding node of nG in SG.
	 * \param nodeLengthG is saving for each node of G its length.
	 * \param nodeLengthSG is assigned for each node of SG its length.
	 * \param edgeLengthG is saving for each edge of G its length.
	 * \param edgeLengthSG is assigned for each edge of SG its length.
	 */
	static void call(const Graph& G, Graph& SG, const node& nG, node& nSG,
			const NodeArray<T>& nodeLengthG, NodeArray<T>& nodeLengthSG,
			const EdgeArray<T>& edgeLengthG, EdgeArray<T>& edgeLengthSG) {
		NodeArray<node> nSG_to_nG(SG);
		EdgeArray<edge> eSG_to_eG(SG);
		call(G, SG, nG, nSG, nSG_to_nG, eSG_to_eG, nodeLengthG, nodeLengthSG, edgeLengthG,
				edgeLengthSG);
	}

	/**
	 * \brief Computes a connected subgraph SG of G containing node nG.
	 * \param G is the original graph.
	 * \param SG is assigned the connected subgraph containing nG.
	 * \param nG is a node in G.
	 * \param nodeLengthG stores for each node of G its length.
	 * \param nodeLengthSG is assigned for each node of SG its length.
	 */
	static void call(const Graph& G, Graph& SG, const node& nG, const NodeArray<T>& nodeLengthG,
			NodeArray<T>& nodeLengthSG) {
		node nSG;
		call(G, SG, nG, nSG, nodeLengthG, nodeLengthSG);
	}

	/**
	 * \brief Computes a connected subgraph SG of G containing node nG.
	 * \param G is the original graph.
	 * \param SG is assigned the connected subgraph containing nG.
	 * \param nG is a node in G.
	 * \param nodeLengthG stores for each node of G its length.
	 * \param nodeLengthSG is assigned for each node of SG its length.
	 * \param edgeLengthG stores for each edge of G its length.
	 * \param edgeLengthSG is assigned for each edge of SG its length.
	 */
	static void call(const Graph& G, Graph& SG, const node& nG, const NodeArray<T>& nodeLengthG,
			NodeArray<T>& nodeLengthSG, const EdgeArray<T>& edgeLengthG, EdgeArray<T>& edgeLengthSG) {
		node nSG = nullptr;
		call(G, SG, nG, nSG, nodeLengthG, nodeLengthSG, edgeLengthG, edgeLengthSG);
	}

	/**
	 * \brief Computes a connected subgraph SG of G containing node nG.
	 * \param G is the original graph.
	 * \param SG is assigned the connected subgraph containing nG.
	 * \param nG is a node in G.
	 * \param nSG
	 * \param nSG_to_nG is mapping nodes in SG to nodes in G
	 * \param eSG_to_eG is mapping edges in SG to edges in G
	 * \param nodeLengthG is saving for each node of G its length.
	 * \param nodeLengthSG is assigned for each node of SG its length.
	 * \param edgeLengthG stores for each edge of G its length.
	 * \param edgeLengthSG is assigned for each edge of SG its length.
	 */
	static void call(const Graph& G, Graph& SG, const node& nG, node& nSG,
			NodeArray<node>& nSG_to_nG, EdgeArray<edge>& eSG_to_eG, const NodeArray<T>& nodeLengthG,
			NodeArray<T>& nodeLengthSG, const EdgeArray<T>& edgeLengthG, EdgeArray<T>& edgeLengthSG) {
		NodeArray<node> nG_to_nSG;
		EdgeArray<edge> eG_to_eSG;
		call(G, SG, nG, nSG, nSG_to_nG, eSG_to_eG, nG_to_nSG, eG_to_eSG, nodeLengthG, nodeLengthSG,
				edgeLengthG, edgeLengthSG);
	}

	/**
	 * \brief Computes a connected subgraph SG of G containing node nG.
	 * \param G is the original graph.
	 * \param SG is assigned the connected subgraph containing nG.
	 * \param nG is a node in G.
	 * \param nSG
	 * \param nSG_to_nG is mapping nodes in SG to nodes in G
	 * \param eSG_to_eG is mapping edges in SG to edges in G
	 * \param nG_to_nSG is mapping nodes in G to nodes in SG
	 * \param eG_to_eSG is mapping edges in G to edges in SG
	 * \param nodeLengthG stores for each node of G its length.
	 * \param nodeLengthSG is assigned for each node of SG its length.
	 * \param edgeLengthG stores for each edge of G its length.
	 * \param edgeLengthSG is assigned for each edge of SG its length.
	 */
	static void call(const Graph& G, Graph& SG, const node& nG, node& nSG,
			NodeArray<node>& nSG_to_nG, EdgeArray<edge>& eSG_to_eG, NodeArray<node>& nG_to_nSG,
			EdgeArray<edge>& eG_to_eSG, const NodeArray<T>& nodeLengthG, NodeArray<T>& nodeLengthSG,
			const EdgeArray<T>& edgeLengthG, EdgeArray<T>& edgeLengthSG);

	/**
	 * \brief Computes a connected subgraph SG of G containing node nG.
	 * \param G is the original graph.
	 * \param SG is assigned the connected subgraph containing nG.
	 * \param nG is a node in G.
	 * \param nSG_to_nG is mapping nodes in SG to nodes in G
	 * \param eSG_to_eG is mapping edges in SG to edges in G
	 * \param nG_to_nSG is mapping nodes in G to nodes in SG
	 * \param eG_to_eSG is mapping edges in G to edges in SG
	 */
	static void call(const Graph& G, Graph& SG, const node& nG, NodeArray<node>& nSG_to_nG,
			EdgeArray<edge>& eSG_to_eG, NodeArray<node>& nG_to_nSG, EdgeArray<edge>& eG_to_eSG);

	/**
	 * \brief Computes a connected subgraph SG of G containing node nG.
	 * \param G is the original graph.
	 * \param SG is assigned the connected subgraph containing nG.
	 * \param nG is a node in G.
	 * \param nSG_to_nG is mapping nodes in SG to nodes in G
	 */
	static void call(const Graph& G, Graph& SG, const node& nG, NodeArray<node>& nSG_to_nG) {
		NodeArray<T> nodeLengthG(G, 0);
		NodeArray<T> nodeLengthSG(SG);
		EdgeArray<T> edgeLengthG(G, 0);
		EdgeArray<T> edgeLengthSG(SG);
		node nSG;
		EdgeArray<edge> eSG_to_eG;
		call(G, SG, nG, nSG, nSG_to_nG, eSG_to_eG, nodeLengthG, nodeLengthSG, edgeLengthG,
				edgeLengthSG);
	}

private:
	/**
	 * \brief Copies to SG node nG and recursively all adjacent edges
	 *  and nodes.
	 *
	 * \param SG is the connected subgraph.
	 * \param nodeVisited saves for all nodes in G, if they were already
	 *    treated.
	 * \param edgeVisited saves for all edges in G, if they were already
	 *    treated.
	 * \param nG is a node in G.
	 * \param nodeLengthG is saving for each node of G its length.
	 * \param nodeLengthSG is saving for each node of SG its length.
	 * \param edgeLengthG is saving for each edge of G its length.
	 * \param edgeLengthSG is saving for each edge of SG its length.
	 * \param nSG_to_nG is mapping nodes in SG to nodes in G
	 * \param eSG_to_eG is mapping edges in SG to edges in G
	 * \param nG_to_nSG is mapping nodes in G to nodes in SG
	 * \param eG_to_eSG is mapping edges in G to edges in SG
	 */
	static void recursion(Graph& SG, NodeArray<bool>& nodeVisited, EdgeArray<bool>& edgeVisited,
			const node& nG, const NodeArray<T>& nodeLengthG, NodeArray<T>& nodeLengthSG,
			const EdgeArray<T>& edgeLengthG, EdgeArray<T>& edgeLengthSG, NodeArray<node>& nSG_to_nG,
			EdgeArray<edge>& eSG_to_eG, NodeArray<node>& nG_to_nSG, EdgeArray<edge>& eG_to_eSG);
};

template<class T>
void ConnectedSubgraph<T>::recursion(Graph& SG, NodeArray<bool>& nodeVisited,
		EdgeArray<bool>& edgeVisited, const node& nG, const NodeArray<T>& nodeLengthG,
		NodeArray<T>& nodeLengthSG, const EdgeArray<T>& edgeLengthG, EdgeArray<T>& edgeLengthSG,
		NodeArray<node>& nSG_to_nG, EdgeArray<edge>& eSG_to_eG, NodeArray<node>& nG_to_nSG,
		EdgeArray<edge>& eG_to_eSG) {
	node nSG = SG.newNode();
	nodeLengthSG[nSG] = nodeLengthG[nG];
	nG_to_nSG[nG] = nSG;
	nSG_to_nG[nSG] = nG;
	nodeVisited[nG] = true;

	for (adjEntry adj : nG->adjEntries) {
		edge eG = adj->theEdge();
		if (!nodeVisited[eG->source()]) {
			recursion(SG, nodeVisited, edgeVisited, eG->source(), nodeLengthG, nodeLengthSG,
					edgeLengthG, edgeLengthSG, nSG_to_nG, eSG_to_eG, nG_to_nSG, eG_to_eSG);
		} else if (!nodeVisited[eG->target()]) {
			recursion(SG, nodeVisited, edgeVisited, eG->target(), nodeLengthG, nodeLengthSG,
					edgeLengthG, edgeLengthSG, nSG_to_nG, eSG_to_eG, nG_to_nSG, eG_to_eSG);
		}
		if (!edgeVisited[eG]) {
			edge eSG = SG.newEdge(nG_to_nSG[eG->source()], nG_to_nSG[eG->target()]);
			edgeLengthSG[eSG] = edgeLengthG[eG];
			eG_to_eSG[eG] = eSG;
			eSG_to_eG[eSG] = eG;
			edgeVisited[eG] = true;
		}
	}
}

template<class T>
void ConnectedSubgraph<T>::call(const Graph& G, Graph& SG, const node& nG, node& nSG,
		NodeArray<node>& nSG_to_nG, EdgeArray<edge>& eSG_to_eG, NodeArray<node>& nG_to_nSG,
		EdgeArray<edge>& eG_to_eSG, const NodeArray<T>& nodeLengthG, NodeArray<T>& nodeLengthSG,
		const EdgeArray<T>& edgeLengthG, EdgeArray<T>& edgeLengthSG) {
	SG.clear();

	NodeArray<bool> nodeVisited(G, false);
	EdgeArray<bool> edgeVisited(G, false);

#if 0
	const int n = G.numberOfNodes();
	const int m = G.numberOfEdges();

	bool* nodeVisited = new bool[n];
	bool* edgeVisited = new bool[m];
	for (int i = 0; i < n; i++)
		nodeVisited[i] = false;
	for (int i = 0; i < m; i++)
		edgeVisited[i] = false;
#endif
	nSG_to_nG.init(SG);
	eSG_to_eG.init(SG);
	nodeLengthSG.init(SG);
	edgeLengthSG.init(SG);
	nG_to_nSG.init(G);
	eG_to_eSG.init(G);

	recursion(SG, nodeVisited, edgeVisited, nG, nodeLengthG, nodeLengthSG, edgeLengthG,
			edgeLengthSG, nSG_to_nG, eSG_to_eG, nG_to_nSG, eG_to_eSG);
	nSG = nG_to_nSG[nG];

#if 0
	delete[] nodeVisited;
	delete[] edgeVisited;
#endif
}

template<class T>
void ConnectedSubgraph<T>::call(const Graph& G, Graph& SG, const node& nG, NodeArray<node>& nSG_to_nG,
		EdgeArray<edge>& eSG_to_eG, NodeArray<node>& nG_to_nSG, EdgeArray<edge>& eG_to_eSG) {
	SG.clear();

	NodeArray<bool> nodeVisited(G, false);
	EdgeArray<bool> edgeVisited(G, false);

#if 0
	bool* nodeVisited = new bool[G.numberOfNodes()];
	bool* edgeVisited = new bool[G.numberOfEdges()];
	for (int i = 0; i < G.numberOfNodes(); i++)
		nodeVisited[i] = false;
	for (int i = 0; i < G.numberOfEdges(); i++)
		edgeVisited[i] = false;
#endif
	nSG_to_nG.init(SG);
	eSG_to_eG.init(SG);
	NodeArray<T> nodeLengthG(G, 0);
	NodeArray<T> nodeLengthSG(SG);
	EdgeArray<T> edgeLengthG(G, 1);
	EdgeArray<T> edgeLengthSG(SG);
	nG_to_nSG.init(G);
	eG_to_eSG.init(G);

	recursion(SG, nodeVisited, edgeVisited, nG, nodeLengthG, nodeLengthSG, edgeLengthG,
			edgeLengthSG, nSG_to_nG, eSG_to_eG, nG_to_nSG, eG_to_eSG);

#if 0
	delete[] nodeVisited;
	delete[] edgeVisited;
#endif
}

}
}
