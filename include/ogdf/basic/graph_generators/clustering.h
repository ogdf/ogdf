/** \file
 * \brief Declaration of randomized clustering generators.
 *
 * \author Carsten Gutwenger, Markus Chimani, Jöran Schierbaum, Simon D. Fink
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
#include <ogdf/basic/basic.h>

#include <cstdint>
#include <iosfwd>

namespace ogdf {

class ClusterGraph;

namespace sync_plan {
class SyncPlan;
}

//! @addtogroup graph-generators
//! @{

//! @name Randomized clustering generators
//! @{

//! Creates a random c-connected clustering for a given graph \p G.
/**
 * The resulting cluster graph is always c-connected but may or may not be c-planar
 * (e.g. when a connected cluster encloses some outside vertices).
 * @param C is a cluster graph for \p G.
 * @param cNum is the maximal number of Clusters introduced.
 * \pre \p G is connected and not empty.
 */
OGDF_EXPORT void randomCConnectedClustering(ClusterGraph& C, int cNum);

//! Creates a random clustering for a given graph \p G.
/**
 * @param C is a cluster graph for \p G.
 * @param cNum is the maximal number of clusters introduced.
 * \pre \p G is connected and not empty.
 */
OGDF_EXPORT void randomClustering(ClusterGraph& C, int cNum);

//! Creates a specified cluster structure for a given graph \p G, and assigns vertices to clusters.
/**
 * This function is called with a graph \p G and the root of a second graph, resembling a tree,
 * that gives the cluster structure. Then, the vertices of G are randomly assigned to the clusters,
 * where we can guarantee that any leaf-cluster has (on average) <i>moreInLeaves</i>-times more vertices
 * than a non-leaf cluster. (E.g. if \p moreInLeaves = 5, any leaf will contain roughly 5 times more vertices than
 * an inner cluster)
 * @param C is a cluster graph for \p G, to be assigned the solution.
 * @param root is a node in some other graph (say \a T). \a T is a tree that we will consider rooted at \p root.
 *        \a T is the pattern for the cluster hierarchy.
 * @param moreInLeaves is a factor such that leaf-clusters have on average <i>moreInLeaves</i>-times more
 *        vertices than inner clusters
 * \pre \p G contains at least twice as many nodes as \a T has leaves.
 */
OGDF_EXPORT void randomClustering(ClusterGraph& C, const node root, int moreInLeaves);

//! Parameters for the randomPlanarClustering() method
struct OGDF_EXPORT RandomClusterConfig {
	//! The hard limit of nodes that can be direct children of a single cluster.
	int max_nodes_in_cluster = 0;
	//! The probability of not adding another node to the current cluster.
	double prob_no_further_node = 0.1;
	//! The probability of not creating a further cluster.
	double prob_no_further_cluster = 0.0;
	//! The hard limit of clusters to create.
	int max_clusters = 0;
	//! The hard minimum of nodes that must remain directly within the root cluster.
	int min_root_nodes = 0;
	//! Whether the resulting clustering must be c-connected.
	bool cconnected = false;
	//! Timeout for the generation in seconds.
	int timeout = 0;

	//! Get the expected number of nodes per cluster, i.e., 1.0 / prob_no_further_node.
	double expected_nodes() const { return 1.0 / prob_no_further_node; }

	//! Set the expected number of nodes per cluster, i.e., prob_no_further_node = 1.0 / \p n.
	void expected_nodes(double n) { prob_no_further_node = 1.0 / n; }

	friend std::ostream& operator<<(std::ostream& os, const RandomClusterConfig& config);
};

//! Creates a random c-planar clustering for a given planar graph \p G.
/**
 * The clusters are created by working on a copy of \p G and making it connected as well as triangulated and then:
 * 1. selecting a random node \p u and putting it into a new cluster
 * 2. selecting a random incident edge \p e such that its other endpoint \p v and \p u have at most
 *    two common neighbors (this ensures c-planarity)
 * 3. adding \p v into the cluster of \p u by contracting \p e in the copy of \p G
 * 4. going to step 2. if the cluster may grow further or otherwise step 1. if we may another cluster
 *
 * The clustering is c-connected if we only select edges that were already present in \p G instead
 * of being added by the triangulation.
 *
 * @param CG is a cluster graph for some \p G.
 * @param config configuration for the random generation
 * @return false if the clustering ran into the configured timeout, true otherwise
 */
OGDF_EXPORT bool randomPlanarClustering(ClusterGraph& CG, const RandomClusterConfig& config);

//! Create a random planar graph with a c-planar clustering.
/**
 * The graph is iteratively created by starting with a random planar connected graph and then,
 * for each cluster, choosing an arbitrary vertex and replacing it with a cluster containing another
 * new random planar graph. The replacement is done by identifying pairs of incident edges between
 * the chosen vertex and another random vertex from the new graph. Thus, when a cut-vertex is joined
 * in this way, the resulting graph will no longer be c-connected.
 *
 * Note that the given parameters are upper bounds that might not be always reached.
 * Especially, replacing a vertex (which needs degree at least 4) with a cluster makes the parent
 * cluster contain one vertex less, and no further clusters can be added if all vertices have degree
 * 3 or less.
 *
 * @param G will be assigned the planar graph.
 * @param CG will be assigned the c-planar clustering of \p G.
 * @param clusters how many clusters to generate.
 * @param node_per_cluster how many nodes each cluster should directly contain.
 * @param edges_per_cluster how many edges to add between each clusters' nodes.
 */
OGDF_EXPORT void randomClusterPlanarGraph(Graph& G, ClusterGraph& CG, int clusters,
		int node_per_cluster, int edges_per_cluster);

//! Create a random SynchronizedPlanarity instance by introducing \p pipe_count pipes between vertices of degree at least \p min_deg.
OGDF_EXPORT void randomSyncPlanInstance(sync_plan::SyncPlan& pq, int pipe_count, int min_deg = 3);

//! Create a (simultaneously planar) 2-SEFE instance with a given shared graph.
/**
 * This works by randomly subdividing the faces of the embedded \p sefe, separately for the two
 * exclusive graphs.
 *
 * @param sefe contains the shared graph and will be modified to also contain the exclusive edges.
 * @param edge_types will be assigned 3 for all shared edges and 1 or 2 for edges in either of the exclusive graphs.
 * @param edges1 the number of edges to create for the first exclusive graph.
 * @param edges2 the number of edges to create for the second exclusive graph.
 * @pre \p sefe must contain a connected, planarly embedded graph
 */
OGDF_EXPORT void randomSEFEInstanceBySharedGraph(Graph* sefe, EdgeArray<uint8_t>& edge_types,
		int edges1, int edges2);

//! Create a (simultaneously planar) 2-SEFE instance with a given union graph.
/**
 *This works by randomly assigning edges of the embedded \p sefe to either an exclusive or the shared graph.
 *
 * @param sefe contains the desired union graph.
 * @param edge_types will be assigned 3 for all shared edges and 1 or 2 for edges in either of the exclusive graphs.
 * @param frac_shared the (expected) fraction of edges to make shared (i.e. type 3).
 * @param frac_g1 the (expected) fraction of edges to assign to exclusive graph 1.
 */
OGDF_EXPORT void randomSEFEInstanceByUnionGraph(const Graph* sefe, EdgeArray<uint8_t>& edge_types,
		double frac_shared = 0.34, double frac_g1 = 0.33);

//! @}
//! @}

}
