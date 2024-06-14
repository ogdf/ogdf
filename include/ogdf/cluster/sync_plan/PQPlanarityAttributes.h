#pragma once

#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/LayoutModule.h>
#include <ogdf/fileformats/GraphIO.h>

#include <memory>

using namespace ogdf;

extern const std::array<Color, 63> colors;

void formatNode(node n, GraphAttributes* ga, int group);

void styleClusterBorder(
		const ClusterGraph& CG, const EdgeArray<List<std::pair<adjEntry, cluster>>>& subdivisions,
		GraphAttributes& GA,
		const std::function<edge(edge)>& translate = [](edge e) -> edge { return e; });

void spreadParallels(GraphAttributes& GA, double min_spread = 0.1, double max_spread = 0.6,
		double max_abs = 100);

std::unique_ptr<std::pair<GraphCopy, GraphAttributes>> drawClusterGraph(ClusterGraph& CG,
		GraphAttributes& GA, adjEntry adjExternal = nullptr);

void fixLoops(Graph& G, const std::function<void(edge, edge)>& cb);

void fixParallels(Graph& G, const std::function<void(edge, edge)>& cb);

void bendEdge(GraphAttributes& GA, edge e, double bend);

class PQPlanarity;

class PQPlanarityDrawer {
	std::unique_ptr<LayoutModule> planar_layout;
	std::unique_ptr<LayoutModule> non_planar_layout;
	GraphIO::SVGSettings svg;

	PQPlanarity* PQ;
	GraphAttributes BC_GA;
	List<edge> g_edges;
	List<edge> bc_edges;
	List<int> reuse_g_edge_idx;
	List<int> reuse_bc_edge_idx;
	std::unique_ptr<GraphAttributes> own_GA;

public:
	PQPlanarityDrawer(PQPlanarity* pq);

	virtual ~PQPlanarityDrawer() { cleanUp(); }

	GraphAttributes& ensureGraphAttributes();

	[[nodiscard]] const GraphAttributes& getBC_GA() const { return BC_GA; }

	[[nodiscard]] const GraphIO::SVGSettings& getSvg() const { return svg; }

	void layout(bool format = true, bool components = true);

	void cleanUp();
};
