#include <ogdf/planarity/PlanarizationLayout.h>
#include <ogdf/planarity/SubgraphPlanarizer.h>
#include <ogdf/planarity/VariableEmbeddingInserter.h>
#include <ogdf/planarity/FastPlanarSubgraph.h>
#include <ogdf/orthogonal/OrthoLayout.h>
#include <ogdf/planarity/EmbedderMinDepthMaxFaceLayers.h>
#include <ogdf/fileformats/GraphIO.h>


using namespace ogdf;

int main()
{
	Graph G;
	GraphAttributes GA(G,
		GraphAttributes::nodeGraphics |
		GraphAttributes::edgeGraphics |
		GraphAttributes::nodeLabel |
		GraphAttributes::edgeStyle |
		GraphAttributes::nodeStyle |
		GraphAttributes::nodeTemplate);
	GraphIO::readGML(GA, G, "ERDiagram.gml");

	PlanarizationLayout pl;

	SubgraphPlanarizer *crossMin = new SubgraphPlanarizer;

	FastPlanarSubgraph *ps = new FastPlanarSubgraph;
	ps->runs(100);
	VariableEmbeddingInserter *ves = new VariableEmbeddingInserter;
	ves->removeReinsert(rrAll);

	crossMin->setSubgraph(ps);
	crossMin->setInserter(ves);

	EmbedderMinDepthMaxFaceLayers *emb = new EmbedderMinDepthMaxFaceLayers;
	pl.setEmbedder(emb);

	OrthoLayout *ol = new OrthoLayout;
	ol->separation(20.0);
	ol->cOverhang(0.4);
	pl.setPlanarLayouter(ol);

	pl.call(GA);

	GraphIO::writeGML(GA, "ERDiagram-layout.gml");

	return 0;
}
