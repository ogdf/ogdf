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
	for (node v : G.nodes) {
		bfs_SPSS(v, G, shortestPathMatrix[v], edgeCosts);
	}
}


void bfs_SPSS(node s, const Graph& G, NodeArray<double>& distanceArray, double edgeCosts)
{
	NodeArray<bool> mark(G, false);
	SListPure<node> bfs;
	bfs.pushBack(s);
	// mark s and set distance to itself 0
	mark[s] = true;
	distanceArray[s] = 0;
	while (!bfs.empty()) {
		node w = bfs.popFrontRet();
		double d = distanceArray[w] + edgeCosts;
		for(adjEntry adj : w->adjEntries) {
			node v = adj->twinNode();
			if (!mark[v]) {
				mark[v] = true;
				bfs.pushBack(v);
				distanceArray[v] = d;
			}
		}
	}
}


double dijkstra_SPAP(const GraphAttributes& GA,
	NodeArray<NodeArray<double> >& shortestPathMatrix)
{
	const Graph& G = GA.constGraph();
	EdgeArray<double> edgeCosts(G);
	double avgCosts = 0;
	for (edge e : G.edges) {
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
	for (node v : G.nodes) {
		dijkstra_SPSS(v, G, shortestPathMatrix[v], edgeCosts);
	}
}


void dijkstra_SPSS(node s, const Graph& G, NodeArray<double>& distance,
	const EdgeArray<double>& edgeCosts)
{
	NodeArray<edge> predecessor;
	Dijkstra<double> sssp;
	sssp.call(G, edgeCosts, s, predecessor, distance);
}


void floydWarshall_SPAP(NodeArray<NodeArray<double> >& shortestPathMatrix,
	const Graph& G)
{
	for (node u : G.nodes) {
		for (node v : G.nodes) {
			for (node w : G.nodes) {
				shortestPathMatrix[v][w] = min(shortestPathMatrix[v][w],
					shortestPathMatrix[u][v] + shortestPathMatrix[u][w]);
			}
		}
	}
}

}
