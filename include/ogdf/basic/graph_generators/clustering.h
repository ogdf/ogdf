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

#include <ogdf/cluster/ClusterGraph.h>

namespace ogdf {

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
 * Creates a cluster-planar clustering for a given planar graph. The clusters are created by working
 * on a copy of \p G and making it connected as well as triangulated and then:
 * 1. selecting a random node \p u and putting it into a new cluster
 * 2. selecting a random incident edge \p e such that its other endpoint \p v and \p u have at most
 *    two common neighbors (this ensures c-planarity)
 * 3. adding \p v into the cluster of \p u by contracting \p e in the copy of \p G
 * 4. going to step 2. if the cluster may grow further or otherwise step 1. if we may another cluster
 *
 * The clustering is c-connected if we only select edges that were already present in \p G instead
 * of being added by the triangulation.
 *
 * @param C is a cluster graph for some \p G.
 * @param config configuration for the random generation
 * @return false if the clustering ran into the configured timeout, true otherwise
 */
OGDF_EXPORT bool randomPlanarClustering(ClusterGraph& CG, const RandomClusterConfig& config);

//! @}
//! @}

}
