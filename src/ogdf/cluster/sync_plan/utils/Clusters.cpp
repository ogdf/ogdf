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


#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/Logger.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/cluster/ClusterGraph.h>
#include <ogdf/cluster/sync_plan/utils/Clusters.h>

#include <functional>
#include <ostream>
#include <utility>

// #include <ogdf/planarlayout/FPPLayout.h>

using namespace std;

namespace ogdf::sync_plan {


bool isClusterPlanarEmbedding(const ClusterGraph& CG) {
	if (!CG.adjAvailable()) {
		return false;
	}
#ifdef OGDF_DEBUG
	CG.consistencyCheck();
#endif
	GraphCopy gcopy(CG.constGraph());
	gcopy.setOriginalEmbedding();
	clusterBorderToEdges(CG, gcopy, nullptr, [&gcopy](edge e) -> edge { return gcopy.copy(e); });
	bool copy_valid = gcopy.representsCombEmbedding();
	//bool cluster_valid = CG.representsCombEmbedding();
	//OGDF_ASSERT(copy_valid == cluster_valid);
	return copy_valid;
}

void clusterBorderToEdges(const ClusterGraph& CG, Graph& G,
		EdgeArray<List<std::pair<adjEntry, cluster>>>* subdivisions,
		const std::function<edge(edge)>& translate) {
	for (cluster c = CG.firstPostOrderCluster(); c != CG.rootCluster(); c = c->pSucc()) {
		adjEntry prev_ray = nullptr, first_ray = nullptr;
		for (adjEntry adj : c->adjEntries) {
			bool reverse = adj->isSource();
			edge the_edge = translate(adj->theEdge());
			if (reverse) {
				G.reverseEdge(the_edge);
			}
			edge new_edge = G.split(the_edge);
			adjEntry spike_to_src =
					the_edge->adjTarget(); // adjacency of the new node toward the source of the_edge
			if (subdivisions != nullptr) {
				if (reverse) {
					(*subdivisions)[adj->theEdge()].emplaceBack(spike_to_src, c);
				} else {
					(*subdivisions)[adj->theEdge()].emplaceFront(spike_to_src, c);
				}
			}
			if (reverse) {
				G.reverseEdge(the_edge);
				G.reverseEdge(new_edge);
			}

			if (prev_ray != nullptr) {
				G.newEdge(prev_ray, spike_to_src->cyclicPred());
			}
			prev_ray = spike_to_src;
			if (first_ray == nullptr) {
				first_ray = spike_to_src;
			}
		}
		if (first_ray != nullptr) {
			G.newEdge(prev_ray, first_ray->cyclicPred());
		}
	}
}

}
