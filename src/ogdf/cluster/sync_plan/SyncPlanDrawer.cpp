/** \file
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
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/LayoutModule.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/Logger.h>
#include <ogdf/basic/PreprocessorLayout.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/geometry.h>
#include <ogdf/basic/graphics.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/cluster/ClusterGraph.h>
#include <ogdf/cluster/sync_plan/PMatching.h>
#include <ogdf/cluster/sync_plan/SyncPlan.h>
#include <ogdf/cluster/sync_plan/SyncPlanComponents.h>
#include <ogdf/cluster/sync_plan/SyncPlanDrawer.h>
#include <ogdf/cluster/sync_plan/basic/Drawing.h>
#include <ogdf/cluster/sync_plan/utils/Logging.h>
#include <ogdf/layered/OptimalHierarchyLayout.h>
#include <ogdf/layered/SugiyamaLayout.h>
#include <ogdf/packing/ComponentSplitterLayout.h>
#include <ogdf/planarlayout/FPPLayout.h>
#include <ogdf/planarlayout/GridLayoutModule.h>

#include <algorithm>
#include <array>
#include <functional>
#include <memory>
#include <sstream>
#include <string>
#include <utility>

using namespace ogdf::sync_plan::internal;

namespace ogdf {

const std::array<Color, 63> colors = {Color("#00FF00"), Color("#0000FF"), Color("#FF0000"),
		Color("#01FFFE"), Color("#FFA6FE"), Color("#FFDB66"), Color("#006401"), Color("#010067"),
		Color("#95003A"), Color("#007DB5"), Color("#FF00F6"), Color("#FFEEE8"), Color("#774D00"),
		Color("#90FB92"), Color("#0076FF"), Color("#D5FF00"), Color("#FF937E"), Color("#6A826C"),
		Color("#FF029D"), Color("#FE8900"), Color("#7A4782"), Color("#7E2DD2"), Color("#85A900"),
		Color("#FF0056"), Color("#A42400"), Color("#00AE7E"), Color("#683D3B"), Color("#BDC6FF"),
		Color("#263400"), Color("#BDD393"), Color("#00B917"), Color("#9E008E"), Color("#001544"),
		Color("#C28C9F"), Color("#FF74A3"), Color("#01D0FF"), Color("#004754"), Color("#E56FFE"),
		Color("#788231"), Color("#0E4CA1"), Color("#91D0CB"), Color("#BE9970"), Color("#968AE8"),
		Color("#BB8800"), Color("#43002C"), Color("#DEFF74"), Color("#00FFC6"), Color("#FFE502"),
		Color("#620E00"), Color("#008F9C"), Color("#98FF52"), Color("#7544B1"), Color("#B500FF"),
		Color("#00FF78"), Color("#FF6E41"), Color("#005F39"), Color("#6B6882"), Color("#5FAD4E"),
		Color("#A75740"), Color("#A5FFD2"), Color("#FFB167"), Color("#009BFF"), Color("#E85EBE")};

void spreadParallels(GraphAttributes& GA, double min_spread, double max_spread, double max_abs) {
	max_abs /= 2;
	GA.clearAllBends();
	EdgeArray<List<edge>> parallels(GA.constGraph());
	getParallelFreeUndirected(GA.constGraph(), parallels);
	for (auto e : GA.constGraph().edges) {
		if (parallels[e].empty()) {
			continue;
		}
		double spread = std::min(min_spread, max_spread / parallels[e].size());
		double width = spread * parallels[e].size();
		parallels[e].pushFront(e);

		DPoint vec = GA.point(e->target()) - GA.point(e->source());
		DPoint mid = GA.point(e->source()) + (vec / 2);
		DPoint ort = vec.orthogonal();
		ort *= vec.norm() / ort.norm();
		if (ort.norm() > max_abs) {
			ort *= max_abs / ort.norm();
		}

		int idx = 0;
		for (auto pe : parallels[e]) {
			OGDF_ASSERT(GA.bends(pe).empty());
			GA.bends(pe).pushBack(mid + ort * (-width / 2 + idx * spread));
			idx++;
		}
	}
}

void fixLoops(Graph& G, const std::function<void(edge, edge)>& cb) {
	List<edge> edges;
	G.allEdges(edges);
	for (edge e : edges) {
		if (e->isSelfLoop()) {
			cb(e, G.split(e));
		}
	}
}

void fixParallels(Graph& G, const std::function<void(edge, edge)>& cb) {
	EdgeArray<SListPure<edge>> parallels(G);
	getParallelFreeUndirected(G, parallels);
	for (const auto& list : parallels) {
		for (auto e : list) {
			cb(e, G.split(e));
		}
	}
}

void bendEdge(GraphAttributes& GA, edge e, double bend) {
	DPoint vec = GA.point(e->target()) - GA.point(e->source());
	DPoint mid = GA.point(e->source()) + (vec / 2);
	DPoint ort = vec.orthogonal();
	ort *= vec.norm() / ort.norm();
	GA.bends(e).pushBack(mid + ort * bend);
}

namespace sync_plan {

void formatNode(node n, GraphAttributes* ga, int group) {
	ga->shape(n) = Shape::Ellipse;
	ga->width(n) = 5;
	ga->height(n) = 5;
	ga->fillColor(n) = Color::Name::Lightgrey;
	ga->strokeColor(n) = colors[group % colors.size()];
	ga->strokeWidth(n) = 2;
	if (ga->label(n).empty()) {
		std::stringstream ss;
		ss << n->index();
		ga->label(n) = ss.str();
	}
}

void styleClusterBorder(const ClusterGraph& CG,
		const EdgeArray<List<std::pair<adjEntry, cluster>>>& subdivisions, GraphAttributes& GA,
		const std::function<edge(edge)>& translate) {
	OGDF_ASSERT(subdivisions.graphOf() == CG.getGraph());
	ClusterArray<bool> seen(CG, false);
	for (edge e : CG.constGraph().edges) {
		if (subdivisions[e].empty()) {
			continue;
		}
		node src = translate(e)->source();
		node tgt = translate(e)->target();
		int i = 1, m = subdivisions[e].size() + 1;
		double dx = (GA.x(tgt) - GA.x(src)) / m, dy = (GA.y(tgt) - GA.y(src)) / m;
		for (std::pair<adjEntry, cluster> subdiv : subdivisions[e]) {
			node bn = subdiv.first->theNode();
			OGDF_ASSERT(bn->degree() == 4);
			edge be = subdiv.first->cyclicSucc()->theEdge();
			cluster c = subdiv.second;
			GA.x(bn) = GA.x(src) + dx * i;
			GA.y(bn) = GA.y(src) + dy * i;
			GA.width(bn) = 1;
			GA.height(bn) = 1;
			GA.strokeColor(be) = colors[c->index() % colors.size()];
			GA.strokeWidth(be) = 2;
			if (!seen[c]) {
				std::stringstream ss;
				ss << "cluster " << c->index();
				GA.label(be) = ss.str();
				seen[c] = true;
			}
			i++;
		}
	}
}

std::unique_ptr<std::pair<GraphCopy, GraphAttributes>> drawClusterGraph(ClusterGraph& CG,
		GraphAttributes& GA, PlanarGridLayoutModule& layout, adjEntry adjExternal) {
	std::unique_ptr<std::pair<GraphCopy, GraphAttributes>> pair(
			new std::pair<GraphCopy, GraphAttributes>());
	GraphCopy& GC = pair->first;
	GraphAttributes& GCA = pair->second;

	GC.init(CG.constGraph());
	GC.setOriginalEmbedding();
	GCA.init(GC, GA.attributes());
	GA.transferToCopy(GCA);

	EdgeArray<List<std::pair<adjEntry, cluster>>> subdivisions(CG.constGraph());
	std::function<edge(edge)> translate = [&GC](edge e) -> edge { return GC.copy(e); };
	planarizeClusterBorderCrossings(CG, GC, &subdivisions, translate);
	styleClusterBorder(CG, subdivisions, GCA, translate);
	auto formatSplitEdge = [&GCA](edge e, edge e2) {
		GCA.width(e2->source()) = GCA.width(e2->target());
		GCA.height(e2->source()) = GCA.height(e2->target());
		GCA.strokeColor(e2) = GCA.strokeColor(e);
		GCA.strokeWidth(e2) = GCA.strokeWidth(e);
	};
	fixLoops(GC, formatSplitEdge);
	fixParallels(GC, formatSplitEdge);
	if (GC.representsCombEmbedding()) {
		if (adjExternal) {
			layout.callFixEmbed(GCA, GC.copy(adjExternal));
		} else {
			layout.callFixEmbed(GCA);
		}
	} else {
		layout.call(GCA);
	}

	return pair;
}

SyncPlanDrawer::SyncPlanDrawer(SyncPlan* pq) : PQ(pq) {
	auto* pre = new PreprocessorLayout();
	planar_layout.reset(pre);
	pre->setRandomizePositions(false);
	auto* split = new ComponentSplitterLayout();
	pre->setLayoutModule(split);
	split->setBorder(150);
	auto* act = new FPPLayout();
	split->setLayoutModule(act);
	act->separation(50);

	auto* sugi = new SugiyamaLayout();
	non_planar_layout.reset(sugi);
	sugi->minDistCC(150);
	sugi->runs(1);
	auto* ohl = new OptimalHierarchyLayout();
	sugi->setLayout(ohl);
	ohl->nodeDistance(50);

	svg.margin(50);
	svg.bezierInterpolation(true);
	svg.curviness(0.3);
}

void SyncPlanDrawer::layout(bool format, bool components) {
	if (PQ->GA != nullptr) {
		if (isPlanar(*PQ->G)) {
			planar_layout->call(*PQ->GA);
		} else {
			PQ->log.lout(Logger::Level::Alarm)
					<< "Graph is non-planar, still trying to generate layout" << std::endl;
			non_planar_layout->call(*PQ->GA);
		}
		ogdf::spreadParallels(*PQ->GA);
		if (format) {
			for (node n : PQ->G->nodes) {
				PQ->formatNode(n);
			}
		}
	}

	if (components) {
		BC_GA.init(PQ->components.bcTree(), GraphAttributes::all);
		planar_layout->call(BC_GA);
		for (node bc : PQ->components.bcTree().nodes) {
			BC_GA.label(bc) = internal::to_string(PQ->components.fmtBCNode(bc));
			formatNode(bc, &BC_GA, bc->index());
			if (PQ->components.isCutComponent(bc)) {
				BC_GA.shape(bc) = Shape::Ellipse;
			} else {
				BC_GA.shape(bc) = Shape::Rect;
			}
		}
	}

	for (const auto& pipe : PQ->matchings) {
		edge e;
		if (!reuse_g_edge_idx.empty()) {
			e = PQ->G->newEdge(pipe.node1, pipe.node2, reuse_g_edge_idx.popFrontRet());
		} else {
			e = PQ->G->newEdge(pipe.node1, pipe.node2);
		}
		if (PQ->GA) {
			PQ->GA->strokeColor(e) = Color::Name::Darkgray;
			PQ->GA->strokeWidth(e) = 3;
			PQ->GA->strokeType(e) = StrokeType::Dash;
			PQ->GA->bends(e).clear();
			bendEdge(*PQ->GA, e, 0.1);
		}
		g_edges.pushBack(e);

		if (components) {
			edge e2;
			node e2f = PQ->components.biconnectedComponent(pipe.node1);
			node e2t = PQ->components.biconnectedComponent(pipe.node2);
			if (!reuse_bc_edge_idx.empty()) {
				e2 = PQ->components.BC.newEdge(e2f, e2t, reuse_bc_edge_idx.popFrontRet());
			} else {
				e2 = PQ->components.BC.newEdge(e2f, e2t);
			}
			BC_GA.strokeColor(e2) = Color::Name::Darkgray;
			BC_GA.strokeType(e2) = StrokeType::Dash;
			bendEdge(BC_GA, e2, 0.1);
			bc_edges.pushBack(e2);
		}
	}
}

void SyncPlanDrawer::cleanUp() {
	for (edge e : g_edges) {
		reuse_g_edge_idx.pushBack(e->index());
		PQ->G->delEdge(e);
	}
	for (edge e : bc_edges) {
		reuse_bc_edge_idx.pushBack(e->index());
		PQ->components.BC.delEdge(e);
	}
	g_edges.clear();
	bc_edges.clear();
	BC_GA.init(0);
}

GraphAttributes& SyncPlanDrawer::ensureGraphAttributes() {
	if (PQ->GA == nullptr) {
		own_GA = std::make_unique<GraphAttributes>(*(PQ->G), GraphAttributes::all);
		PQ->GA = own_GA.get();
	}
	return *PQ->GA;
}

}
}
