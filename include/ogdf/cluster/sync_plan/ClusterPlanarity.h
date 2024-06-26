/** \file
* \brief TODO Document
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

#include <ogdf/cluster/ClusterPlanarityModule.h>

namespace ogdf {
/**
 * Warning: the destructive methods invalidate all node and edge objects even if the instance is c-planar
 */
class OGDF_EXPORT SyncPlanClusterPlanarityModule : public ClusterPlanarityModule {
	std::vector<std::pair<adjEntry, adjEntry>>* m_augmentation;

public:
	bool isClusterPlanarDestructive(ClusterGraph& CG, Graph& G) override;
	bool clusterPlanarEmbedClusterPlanarGraph(ClusterGraph& CG, Graph& G) override;
	bool clusterPlanarEmbed(ClusterGraph& CG, Graph& G)  override;;

	//! When set to a non-null pointer, will contain the augmentation edges to make the graph c-connected c-plane after calling clusterPlanarEmbed().
	void setStoredAugmentation(std::vector<std::pair<adjEntry, adjEntry>>* augmentation) {
		m_augmentation = augmentation;
	}

	std::vector<std::pair<adjEntry, adjEntry>>* getStoredAugmentation() const {
		return m_augmentation;
	}
};

//! Perform the reduction from level- to cluster planarity.
/**
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
}
