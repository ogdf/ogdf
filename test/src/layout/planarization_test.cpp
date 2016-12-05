//*********************************************************
// Tested classes:
//    - PlanarizationLayout
//    - PlanarizationGridLayout
//
//  Author: Carsten Gutwenger, Tilo Wiedera
//*********************************************************

#include <bandit/bandit.h>

#include <ogdf/planarity/PlanarizationLayout.h>
#include <ogdf/planarity/PlanarizationGridLayout.h>
#include <ogdf/planarity/SubgraphPlanarizer.h>
#include <ogdf/planarity/PlanarSubgraphFast.h>
#include <ogdf/planarity/VariableEmbeddingInserter.h>
#include <ogdf/planarity/FixedEmbeddingInserter.h>
#include <ogdf/planarlayout/MixedModelLayout.h>
#include <ogdf/planarlayout/MMCBLocalStretch.h>

#include "layout_helpers.h"

using namespace ogdf;

go_bandit([](){ bandit::describe("Planarization layouts", [](){
	PlanarizationLayout pl, plFixed;
	PlanarizationGridLayout pgl, pglMM;

	VariableEmbeddingInserter *pVarInserter = new VariableEmbeddingInserter;
	FixedEmbeddingInserter *pFixInserter = new FixedEmbeddingInserter;

	SubgraphPlanarizer *pCrossMin = new SubgraphPlanarizer;
	pCrossMin->setSubgraph(new PlanarSubgraphFast);
	pCrossMin->setInserter(pVarInserter);
	pCrossMin->permutations(4);

	pl.setCrossMin(pCrossMin->clone());
	pgl.setCrossMin(pCrossMin->clone());
	pglMM.setCrossMin(pCrossMin->clone());

	pCrossMin->setInserter(pFixInserter);
	pCrossMin->permutations(1);
	plFixed.setCrossMin(pCrossMin);

	MixedModelLayout *pMml = new MixedModelLayout;
	pMml->setCrossingsBeautifier(new MMCBLocalStretch);
	pglMM.setPlanarLayouter(pMml);

	describeLayoutModule("planarization layout", pl, GraphAttributes::edgeType | GraphAttributes::nodeType, GR_ALL, 50);
	describeLayoutModule("planarization layout with fixed inserter", plFixed, GraphAttributes::edgeType | GraphAttributes::nodeType, GR_ALL, 50);
	describeGridLayoutModule("planarization grid layout", pgl, GR_ALL, 50);
	describeGridLayoutModule("planarization grid layout with mixed model", pglMM, GR_ALL, 50);
}); });
