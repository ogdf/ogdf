/*
 * $Revision: 3396 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-15 14:49:03 +0200 (Mo, 15. Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of several shortest path algorithms.
 *
 * \author Mark Ortmann, University of Konstanz
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

#include <ogdf/graphalg/ShortestPathAlgorithms.h>
#include <ogdf/graphalg/Dijkstra.h>

namespace ogdf {

void bfs_SPAP(const Graph& G, NodeArray<NodeArray<double> >& shortestPathMatrix,
		double edgeCosts)
{
	node v;
	forall_nodes(v, G) {
		bfs_SPSS(v, G, shortestPathMatrix[v], edgeCosts);
	}
}

void bfs_SPSS(const node& v, const Graph& G, NodeArray<double>& distanceArray,
		double edgeCosts)
{
	NodeArray<bool> mark(G, false);
	SListPure<node> bfs;
	bfs.pushBack(v);
	// mark v and set distance to itself 0
	mark[v] = true;
	distanceArray[v] = 0;
	node w;
	node adj;
	edge e;
	while (!bfs.empty()) {
		w = bfs.popFrontRet();
		double d = distanceArray[w] + edgeCosts;
		forall_adj_edges(e,w){
		adj = e->opposite(w);
		if (!mark[adj]) {
			mark[adj] = true;
			bfs.pushBack(adj);
			distanceArray[adj] = d;
		}
	}
}
}

double dijkstra_SPAP(const GraphAttributes& GA,
		NodeArray<NodeArray<double> >& shortestPathMatrix)
{
	const Graph& G = GA.constGraph();
	EdgeArray<double> edgeCosts(G);
	edge e;
	double avgCosts = 0;
	forall_edges(e,G) {
		edgeCosts[e] = GA.doubleWeight(e);
		avgCosts += edgeCosts[e];
	}
	dijkstra_SPAP(G, shortestPathMatrix, edgeCosts);
	return avgCosts / G.numberOfEdges();
}

void dijkstra_SPAP(const Graph& G,
		NodeArray<NodeArray<double> >& shortestPathMatrix,
		const EdgeArray<double>& edgeCosts)
{
	node v;
	forall_nodes(v, G) {
		dijkstra_SPSS(v, G, shortestPathMatrix[v], edgeCosts);
	}
}

void dijkstra_SPSS(node s, const Graph& G, NodeArray<double>& distance,
		const EdgeArray<double>& edgeCosts)
{
	NodeArray<edge> predecessor(G);
	Dijkstra<double> sssp;
	sssp.call(G, edgeCosts, s, predecessor, distance);
}

void floydWarshall_SPAP(NodeArray<NodeArray<double> >& shortestPathMatrix,
		const Graph& G)
{
	node u;
	forall_nodes(u, G) {
		node v;
		forall_nodes(v, G) {
			node w;
			forall_nodes(w, G) {
				shortestPathMatrix[v][w] = min(shortestPathMatrix[v][w],
						shortestPathMatrix[u][v] + shortestPathMatrix[u][w]);
			}
		}
	}
}

}
