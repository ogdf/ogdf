/** \file
 * \brief Declaration of ClusterPlanarityModule which implements a cluster-planarity
 *        test and, optionally, an embedder.
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

#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/Module.h>
#include <ogdf/cluster/ClusterGraph.h>

namespace ogdf {

class OGDF_EXPORT ClusterPlanarityModule : public Module {
public:
	ClusterPlanarityModule() = default;

	virtual ~ClusterPlanarityModule() = default;

	//! Returns true, if \p CG is cluster-planar, false otherwise.
	virtual bool isClusterPlanar(const ClusterGraph& CG) {
		Graph Gcopy;
		ClusterGraph CGcopy(CG, Gcopy);
		return isClusterPlanarDestructive(CGcopy, Gcopy);
	};

	//! Returns true, if \p CG is cluster-planar, false otherwise. In it is non-cluster-planar, the (Cluster)Graph may be arbitrarily changed after the call.
	/**
	 * This variant may be slightly faster than the default isClusterPlanar
	 */
	virtual bool isClusterPlanarDestructive(ClusterGraph& CG, Graph& G) = 0;

	//! Returns true, if \p CG is cluster-planar, false otherwise. If true, \p CG contains a cluster-planar embedding.
	virtual bool clusterPlanarEmbed(ClusterGraph& CG, Graph& G) {
		GraphCopy Gcopy;
		Gcopy.setOriginalGraph(G);

		ClusterArray<cluster> origC;
		NodeArray<node> origN;
		EdgeArray<edge> origE;
		ClusterGraph CGcopy(CG, Gcopy, origC, origN, origE);

		if (!clusterPlanarEmbedClusterPlanarGraph(CGcopy, Gcopy)) {
			return false;
		}
		CG.adjAvailable(true);
		Gcopy.setOriginalEmbedding();
		for (cluster c : CGcopy.clusters) {
			origC[c]->adjEntries.clear();
			for (adjEntry adj : c->adjEntries) {
				origC[c]->adjEntries.pushBack(Gcopy.original(adj));
			}
		}
		return true;
	};

	//! Constructs a cluster-planar embedding of \p CG. \p CG \b has to be cluster-planar!
	/**
	 * Returns true if the embedding was successful.
	 * Returns false if the given graph was non-cluster-planar
	 * (and leaves the (Cluster)Graph in an at least partially invalidated state).
	 *
	 * This routine may be slightly faster than clusterPlanarEmbed, but requires \p CG to be cluster-planar.
	 * If \p CG is not cluster-planar, the (Cluster)Graph will be (partially) destroyed while trying to embed it!
	 */
	virtual bool clusterPlanarEmbedClusterPlanarGraph(ClusterGraph& CG, Graph& G) {
		throw std::runtime_error(
				"Embedding is (currently) not implemented by this ClusterPlanarityModule!");
	}
};

}