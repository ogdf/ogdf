#include <ogdf/basic/graph_generators.h>
#include <ogdf/planarity/SubgraphPlanarizer.h>
#include <ogdf/planarity/PlanarSubgraphFast.h>
#include <ogdf/planarity/VariableEmbeddingInserter.h>
#include <ogdf/fileformats/GraphIO.h>

using namespace ogdf;

int main()
{
	Graph G;
	randomSimpleGraph(G, 100, 150);

	SubgraphPlanarizer SP;
	SP.setSubgraph(new PlanarSubgraphFast);
	SP.setInserter(new VariableEmbeddingInserter);

	int crossNum;
	PlanRep PR(G);
	SP.call(PR, 0, crossNum);

	cout << crossNum << " crossings" << endl;
	GraphIO::write(PR, "output-plan.gml", GraphIO::writeGML);

	return 0;
}
