#include <ogdf/energybased/FMMMLayout.h>
#include <ogdf/fileformats/GraphIO.h>

using namespace ogdf;

int main()
{
	Graph G;
	GraphAttributes GA(G);
	if (!GraphIO::readGML(G, "sierpinski_04.gml") ) {
		cerr << "Could not load sierpinski_04.gml" << endl;
		return 1;
	}

	node v;
	forall_nodes(v,G)
		GA.width(v) = GA.height(v) = 10.0;

	FMMMLayout fmmm;

	fmmm.useHighLevelOptions(true);
	fmmm.unitEdgeLength(15.0);
	fmmm.newInitialPlacement(true);
	fmmm.qualityVersusSpeed(FMMMLayout::qvsGorgeousAndEfficient);

	fmmm.call(GA);
	GraphIO::writeGML(GA, "sierpinski_04-layout.gml");

	return 0;
}
