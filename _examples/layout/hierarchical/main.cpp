#include <ogdf/layered/SugiyamaLayout.h>
#include <ogdf/layered/OptimalRanking.h>
#include <ogdf/layered/MedianHeuristic.h>
#include <ogdf/layered/OptimalHierarchyLayout.h>
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
	if (!GraphIO::readGML(GA, G, "unix-history.gml") ) {
		cerr << "Could not load unix-history.gml" << endl;
		return 1;
	}

	SugiyamaLayout SL;
	SL.setRanking(new OptimalRanking);
	SL.setCrossMin(new MedianHeuristic);

	OptimalHierarchyLayout *ohl = new OptimalHierarchyLayout;
	ohl->layerDistance(30.0);
	ohl->nodeDistance(25.0);
	ohl->weightBalancing(0.8);
	SL.setLayout(ohl);

	SL.call(GA);
	GraphIO::writeGML(GA, "unix-history-layout.gml");

	return 0;
}
