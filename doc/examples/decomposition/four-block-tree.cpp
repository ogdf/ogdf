#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/graph_generators/randomized.h>
#include <ogdf/basic/graphics.h>
#include <ogdf/decomposition/FourBlockTree.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/planarlayout/PlanarStraightLayout.h>
#include <memory>
#include <string>

using namespace ogdf;

int main(void) {
	Graph g;
	constexpr int n = 16;
	randomPlanarConnectedGraph(g, n, 3 * n - 6);
	const adjEntry externalFace = g.firstNode()->firstAdj()->cyclicSucc();

	const auto fbt = FourBlockTree::construct(g, externalFace);

	PlanarStraightLayout layout;
	{
		GraphAttributes ga(g,
				ogdf::GraphAttributes::nodeLabel | ogdf::GraphAttributes::nodeGraphics
						| ogdf::GraphAttributes::nodeStyle | ogdf::GraphAttributes::edgeGraphics
						| ogdf::GraphAttributes::edgeStyle);
		ga.directed() = false;
		layout.callFixEmbed(ga, externalFace);
		for (const node v : g.nodes) {
			ga.label(v) = std::to_string(v->index());
			ga.fillColor(v) = Color::Name::White;
		}
		GraphIO::write(ga, "output-g.svg", GraphIO::drawSVG);
	}
	int i = 0;
	fbt->preorder([&](const FourBlockTree& treeNode) -> void {
		GraphAttributes ga(*treeNode.g,
				ogdf::GraphAttributes::nodeLabel | ogdf::GraphAttributes::nodeGraphics
						| ogdf::GraphAttributes::nodeStyle | ogdf::GraphAttributes::edgeGraphics
						| ogdf::GraphAttributes::edgeStyle);
		ga.directed() = false;
		layout.callFixEmbed(ga, treeNode.externalFace);
		for (const node v : treeNode.g->nodes) {
			ga.label(v) = std::to_string(treeNode.originalNodes[v]->index());
			ga.fillColor(v) = Color::Name::White;
		}
		GraphIO::write(ga, std::string("output-node-") + std::to_string(i++) + ".svg",
				GraphIO::drawSVG);
	});
	return 0;
}
