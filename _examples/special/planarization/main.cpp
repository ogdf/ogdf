#include <ogdf/basic/graph_generators.h>
#include <ogdf/planarity/SubgraphPlanarizer.h>
#include <ogdf/planarity/FastPlanarSubgraph.h>
#include <ogdf/planarity/VariableEmbeddingInserter.h>
#include <ogdf/fileformats/GraphIO.h>

using namespace ogdf;

int main()
{
	Graph G;
	randomSimpleGraph(G, 100, 150);
	PlanRep PR(G);

	int crossNum;
	SubgraphPlanarizer SP;
	SP.setSubgraph(new FastPlanarSubgraph);
	SP.setInserter(new VariableEmbeddingInserter);
	SP.call(PR,0,crossNum);

	cout << crossNum << " crossings" << endl;
	GraphIO::writeGML(PR, "plan.gml");

	return 0;
}
