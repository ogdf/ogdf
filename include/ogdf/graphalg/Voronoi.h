/*
 * $Revision: 3624 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-07-02 13:57:06 +0200 (Di, 02. Jul 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of Voronoi regions in an EdgeWeightedGraph
 *
 * \author Stephan Beyer
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

#ifdef _MSC_VER
#pragma once
#endif

#ifndef OGDF_VORONOI_H
#define OGDF_VORONOI_H


#include <ogdf/basic/List.h>
#include <ogdf/graphalg/Dijkstra.h>
#include <map>

namespace ogdf {

template<typename T>
class Voronoi {
protected:
	NodeArray<edge> m_predecessor;
	NodeArray<T> m_distance;
	NodeArray<node> m_seedOfNode;
	std::map< node, List<node> > m_nodeList;

	void computeVoronoiRegions(const Graph &G, const EdgeArray<T> &weights, const List<node> &seeds);

public:
	//! Build data structure to query Voronoi regions of edge-weighted graph G with given Voronoi seeds
	//! @param G the input graph
	//! @param weights edge weights
	//! @param seeds a list of Voronoi seed nodes, the centers of the Voronoi regions
	Voronoi(const Graph &G, const EdgeArray<T> &weights, const List<node> &seeds)
	 : m_predecessor(G)
	 , m_distance(G)
	 , m_seedOfNode(G)
	 , m_nodeList()
	{
		computeVoronoiRegions(G, weights, seeds);
	}

	//! Returns the edge incident to v and its predecessor.  Note that the predecessor of a terminal is NULL.
	edge predecessorEdge(node v) const
	{
		return m_predecessor[v];
	}

	//! Returns the nearest node to v on the shortest path to its Voronoi seed.
	node predecessor(node v) const
	{
		edge tmp = predecessorEdge(v);
		return (tmp ? tmp->opposite(v) : NULL);
	}

	//! Returns the distance between v and its Voronoi seed.
	T distance(node v) const
	{
		return m_distance[v];
	}

	//! Returns the Voronoi seed of node v.
	node seed(node v) const
	{
		return m_seedOfNode[v];
	}

	//! Returns the list of nodes in the Voronoi region of node v.
	const List<node> &nodesInRegion(node v) const
	{
		return m_nodeList.find(seed(v))->second;
	}
};

//// implementation

template<typename T>
void Voronoi<T>::computeVoronoiRegions(const Graph &G, const EdgeArray<T> &weights, const List<node> &seeds)
{
	Dijkstra<T> sssp;
	sssp.call(G, weights, seeds, m_predecessor, m_distance);

	// extract Voronoi seeds for each node and Voronoi regions for each seed
	NodeArray<bool> processed(G, false);
	forall_listiterators(node, seed, seeds) {
		processed[*seed] = true;
		m_seedOfNode[*seed] = *seed;
		m_nodeList[*seed].pushBack(*seed);
	}

	node u;
	forall_nodes(u, G) {
		List<node> foundNodes;
		node v;
		for (v = u; !processed[v]; v = predecessor(v)) {
			processed[v] = true;
			foundNodes.pushBack(v);
		}
		OGDF_ASSERT(processed[v]);
		forall_listiterators(node, passedNode, foundNodes) {
			m_seedOfNode[*passedNode] = m_seedOfNode[v];
			m_nodeList[m_seedOfNode[v]].pushBack(*passedNode);
		}
	}
}

} // end namespace ogdf

#endif // OGDF_VORONOI_H
