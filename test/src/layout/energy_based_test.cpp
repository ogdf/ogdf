//*********************************************************
// Tested classes:
//    - FMMMLayout
//    - SpringEmbedderGridVariant
//    - GEMLayout
//    - DavidsonHarelLayout
//    - PivotMDS
//
//  Author: Carsten Gutwenger, Tilo Wiedera
//*********************************************************

#include <bandit/bandit.h>

#include <ogdf/energybased/FMMMLayout.h>
#include <ogdf/energybased/SpringEmbedderGridVariant.h>
#include <ogdf/energybased/GEMLayout.h>
#include <ogdf/energybased/DavidsonHarelLayout.h>
#include <ogdf/energybased/PivotMDS.h>
#include <ogdf/basic/graph_generators.h>
#include <ogdf/basic/simple_graph_alg.h>

#include "layout_helpers.h"

using namespace ogdf;

go_bandit([](){ bandit::describe("Energy-based layouts", [](){
	FMMMLayout                fmmm, fmmmNice, fmmmHQ;
	SpringEmbedderGridVariant frl, frlHQ;
	GEMLayout                 gem;
	DavidsonHarelLayout       dhl;
	PivotMDS                  pmds;

	fmmmHQ.qualityVersusSpeed(FMMMLayout::qvsGorgeousAndEfficient);
	fmmmHQ.useHighLevelOptions(true);
	fmmmNice.qualityVersusSpeed(FMMMLayout::qvsNiceAndIncredibleSpeed);
	frlHQ.iterations(1000);

	describeLayoutModule("Fast Multipole Multilevel Embedder", fmmm);
	describeLayoutModule("Fast Multipole Multilevel Embedder with high quality settings", fmmmHQ);
	describeLayoutModule("Fast Multipole Multilevel Embedder with nice quality and incredible speed", fmmmHQ);
	describeLayoutModule("Spring Embedder grid variant", frl);
	describeLayoutModule("Spring Embedder grid variant with high quality settings", fmmmHQ);
	describeLayoutModule("GEM layout", gem);
	describeLayoutModule("Davidson-Harel layout", dhl);
	describeLayoutModule("PivotMDS layout", pmds, 0, GR_CONNECTED);
}); });
