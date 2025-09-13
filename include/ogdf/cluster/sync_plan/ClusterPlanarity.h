/** \file
 * \brief Utilities for reducing from Cluster Planarity to SyncPlan.
 *
 * \author Simon D. Fink <ogdf@niko.fink.bayern>
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
#include <ogdf/basic/GraphSets.h> // IWYU pragma: keep
#include <ogdf/basic/basic.h>
#include <ogdf/cluster/ClusterGraph.h>
#include <ogdf/cluster/ClusterPlanarityModule.h>

#include <utility>
#include <vector>

namespace ogdf {

//! ClusterPlanarity testing in quadratic time using the \ref ogdf::sync_plan::SyncPlan "Synchronized Planarity" approach.
class OGDF_EXPORT SyncPlanClusterPlanarityModule : public ClusterPlanarityModule {
	std::vector<std::pair<adjEntry, adjEntry>>* m_augmentation = nullptr;

public:
	bool isClusterPlanar(const ClusterGraph& CG) override;
	bool isClusterPlanarDestructive(ClusterGraph& CG, Graph& G) override;
	bool clusterPlanarEmbedClusterPlanarGraph(ClusterGraph& CG, Graph& G) override;

	//! When set to a non-null pointer, will contain the augmentation edges to make the graph c-connected c-plane after calling clusterPlanarEmbed().
	/**
	 * @sa insertAugmentationEdges()
	 * @sa SyncPlan::SyncPlan(Graph*, ClusterGraph*, std::vector<std::pair<adjEntry, adjEntry>>*, ClusterGraphAttributes*)
	 */
	void setStoreAugmentation(std::vector<std::pair<adjEntry, adjEntry>>* augmentation) {
		m_augmentation = augmentation;
	}

	std::vector<std::pair<adjEntry, adjEntry>>* getStoreAugmentation() const {
		return m_augmentation;
	}

protected:
	void copyBackEmbedding(ClusterGraph& CG, Graph& G, const ClusterGraph& CGcopy, const Graph& Gcopy,
			const ClusterArray<cluster, true>& copyC, const NodeArray<node, true>& copyN,
			const EdgeArray<edge, true>& copyE, const EdgeArray<edge, true>& origE) const override;
};

//! Perform the reduction from level- to cluster planarity.
/**
 * See Section 6.4.4 of https://doi.org/10.15475/cpatp.2024
 *
 * @param LG the graph that should be tested for level planarity.
 * @param emb the level assignment, containing a list of its contained nodes for each level.
 * @param G will be assigned the graph resulting from the reduction.
 * @param CG will be assigned the clustering resulting from the reduction.
 * @param embMap maps from some edges in \p G to the nodes of \p LG,
 * 		  the order in which they cross the borders of \p CG will induce a cluster planar embedding
 */
OGDF_EXPORT void reduceLevelPlanarityToClusterPlanarity(const Graph& LG,
		const std::vector<std::vector<node>>& emb, Graph& G, ClusterGraph& CG,
		EdgeArray<node>& embMap);

//! Inserts augmentation edges to make a c-plane graph c-connected while maintaining the combinatorial embedding.
/**
 * @param CG the ClusterGraph.
 * @param G the corresponding Graph.
 * @param augmentation a list of adjEntry pairs pointing to (before) the start and (after) the end point of the augmentation edges.
 * @param added if non-null, will be assigned all newly added edges.
 * @param embedded whether CG represents a c-plane embedding that need to be maintained throughout the insertion.
 * @param assert_minimal whether we should assert that the set of augmentation edges is minimal for c-connectivity.
 * @sa SyncPlan::SyncPlan(Graph*, ClusterGraph*, std::vector<std::pair<adjEntry, adjEntry>>*, ClusterGraphAttributes*)
 */
OGDF_EXPORT void insertAugmentationEdges(const ClusterGraph& CG, Graph& G,
		std::vector<std::pair<adjEntry, adjEntry>>& augmentation, EdgeSet* added = nullptr,
		bool embedded = true, bool assert_minimal = true);

}
