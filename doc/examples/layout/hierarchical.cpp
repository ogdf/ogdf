#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/layered/MedianHeuristic.h>
#include <ogdf/layered/OptimalHierarchyLayout.h>
#include <ogdf/layered/OptimalRanking.h>
#include <ogdf/layered/SugiyamaLayout.h>
#include <iostream>
#include <string>

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
	if (!GraphIO::read(GA, G, "unix-history.gml")) {
		std::cerr << "Could not load unix-history.gml" << std::endl;
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
	GraphIO::write(GA, "output-unix-history-hierarchical.gml");
	GraphIO::write(GA, "output-unix-history-hierarchical.svg");

	return 0;
}
