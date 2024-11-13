#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/hypergraph/EdgeStandardRep.h>
#include <ogdf/hypergraph/Hypergraph.h>
#include <ogdf/hypergraph/HypergraphAttributes.h>
#include <ogdf/hypergraph/HypergraphLayout.h>
#include <string>

using namespace ogdf;

int main()
{
	Hypergraph H;

	H.readBenchHypergraph("c17.bench");

	HypergraphAttributesES HA(H, EdgeStandardType::tree);
	HypergraphLayoutES hlES;

	hlES.setProfile(HypergraphLayoutES::Profile::Normal);
	hlES.call(HA);

	GraphIO::write(HA.repGA(), "output-c17.gml", GraphIO::writeGML);
	GraphIO::write(HA.repGA(), "output-c17.svg", GraphIO::drawSVG);

	return 0;
}
