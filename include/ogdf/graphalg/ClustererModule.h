/** \file
 * \brief Declaration of interface for clustering algorithms that
 * 		  compute a clustering for a given graph based on some
 *        structural or semantical properties.
 *
 * Precondition:
 *    Input graph has to be connected
 *
 * \author Karsten Klein
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
#include <ogdf/basic/SList.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/memory.h>
#include <ogdf/basic/simple_graph_alg.h>

namespace ogdf {
class ClusterGraph;

//Helper classes to code a clustering, used as an interface to applications that
//need to build the clustergraph structure themselves
class SimpleCluster {
public:
	SimpleCluster(SimpleCluster* parent = nullptr) : m_size(0), m_parent(parent), m_index(-1) { }

	//insert vertices and children
	void pushBackVertex(node v) { m_nodes.pushBack(v); }

	void pushBackChild(SimpleCluster* s) { m_children.pushBack(s); }

	void setParent(SimpleCluster* parent) { m_parent = parent; }

	SimpleCluster* getParent() { return m_parent; }

	void setIndex(int index) { m_index = index; }

	int getIndex() { return m_index; }

	SList<node>& nodes() { return m_nodes; }

	SList<SimpleCluster*>& children() { return m_children; }

	int m_size; //preliminary: allowed to be inconsistent with real vertex number to
	//allow lazy vertex addition (should therefore be local Array?)

private:
	SList<node> m_nodes;
	SList<SimpleCluster*> m_children;
	SimpleCluster* m_parent;
	int m_index; //index of the constructed cluster
};

/**
 * \brief Interface for algorithms that compute a clustering for a
 *      given graph
 *
 * The class ClustererModule is the base class for clustering
 * classes that allow  to compute some hierarchical clustering
 */
class OGDF_EXPORT ClustererModule {
public:
	//Constructor taking a graph G to be clustered
	explicit ClustererModule(const Graph& G) : m_pGraph(&G) { OGDF_ASSERT(isConnected(G)); }

	//! Default constructor, initializes a clustering module.
	// Allows to cluster multiple
	// graphs with the same instance of the Clusterer
	ClustererModule() { }

	/**
	 * \brief Sets the graph to be clustered
	 *
	 * @param G is the input graph
	 *
	 */
	void setGraph(const Graph& G) {
		OGDF_ASSERT(isConnected(G));
		m_pGraph = &G;
	}

	//! Returns the graph to be clustered
	const Graph& getGraph() const { return *m_pGraph; }

	/**
	 * \brief compute some kind of clustering on the graph m_pGraph
	 *
	 * This is the algorithm call that has to be implemented by derived classes
	 *
	 * @param sl is the resulting list of clusters
	 */
	virtual void computeClustering(SList<SimpleCluster*>& sl) = 0;

	//! translate computed clustering into cluster hierarchy in cluster graph C
	virtual void createClusterGraph(ClusterGraph& C) = 0;

	//! compute a clustering index for each vertex
	virtual double computeCIndex(const Graph& G, node v) = 0;

	//! compute a clustering index for each vertex
	virtual double computeCIndex(node v) = 0;

	//! compute the average clustering index for the given graph
	virtual double averageCIndex() { return averageCIndex(*m_pGraph); }

	virtual double averageCIndex(const Graph& G) {
		double ciSum = 0.0;
		for (node v : G.nodes) {
			ciSum += computeCIndex(G, v);
		}
		return ciSum / (G.numberOfNodes());
	}

protected:
	const Graph* m_pGraph; //the graph to be clustered

	OGDF_MALLOC_NEW_DELETE
};

}
