#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/hypergraph/Hypergraph.h>
#include <ogdf/hypergraph/HypergraphLayout.h>
#include <ogdf/hypergraph/HypergraphAttributes.h>

using namespace ogdf;

int main()
{
	Hypergraph H;

	H.readBenchHypergraph("c17.bench");

	HypergraphAttributesES HA(H, EdgeStandardType::tree);
	HypergraphLayoutES hlES;

	hlES.setProfile(HypergraphLayoutES::Normal);
	hlES.call(HA);

	GraphIO::writeGML(HA.repGA(), "c17.gml");

	return 0;
}
