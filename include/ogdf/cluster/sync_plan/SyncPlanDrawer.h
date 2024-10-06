/** \file
 * \brief Utilities by dumping a drawing of the current state of a SyncPlan instance.
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
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/LayoutModule.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/basic.h>
#include <ogdf/cluster/ClusterGraph.h>
#include <ogdf/fileformats/GraphIO.h>

#include <functional>
#include <memory>
#include <utility>

namespace ogdf {
class PlanarGridLayoutModule;
} // namespace ogdf

namespace ogdf::sync_plan {
class SyncPlan;

//! Simple util for apply a default style to nodes, including a group-based coloring.
OGDF_EXPORT void formatNode(node n, GraphAttributes* ga, int group);

//! Properly style a cluster border inserted by planarizeClusterBorderCrossings()
OGDF_EXPORT void styleClusterBorder(
		const ClusterGraph& CG, const EdgeArray<List<std::pair<adjEntry, cluster>>>& subdivisions,
		GraphAttributes& GA,
		const std::function<edge(edge)>& translate = [](edge e) -> edge { return e; });

//! Draw a graph using a given planar layout by turning its cluster borders into edges via planarizeClusterBorderCrossings()
OGDF_EXPORT std::unique_ptr<std::pair<GraphCopy, GraphAttributes>> drawClusterGraph(ClusterGraph& CG,
		GraphAttributes& GA, PlanarGridLayoutModule& layout, adjEntry adjExternal = nullptr);

//! Utilities by dumping a drawing of the current state of a SyncPlan instance.
class OGDF_EXPORT SyncPlanDrawer {
	std::unique_ptr<LayoutModule> planar_layout;
	std::unique_ptr<LayoutModule> non_planar_layout;
	GraphIO::SVGSettings svg;

	SyncPlan* PQ;
	GraphAttributes BC_GA;
	List<edge> g_edges;
	List<edge> bc_edges;
	List<int> reuse_g_edge_idx;
	List<int> reuse_bc_edge_idx;
	std::unique_ptr<GraphAttributes> own_GA;

public:
	SyncPlanDrawer(SyncPlan* pq);

	virtual ~SyncPlanDrawer() { cleanUp(); }

	GraphAttributes& ensureGraphAttributes();

	[[nodiscard]] const GraphAttributes& getBC_GA() const { return BC_GA; }

	[[nodiscard]] const GraphIO::SVGSettings& getSvg() const { return svg; }

	void layout(bool format = true, bool components = true);

	void cleanUp();
};
}
