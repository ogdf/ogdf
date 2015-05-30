//*********************************************************
//  Tested classes:
//    - PlanarStraightLayout
//    - PlanarDrawLayout
//    - MixedModelLayout
//
//  Author: Carsten Gutwenger, Tilo Wiedera
//*********************************************************

#include <bandit/bandit.h>

#include <ogdf/planarlayout/PlanarStraightLayout.h>
#include <ogdf/planarlayout/PlanarDrawLayout.h>
#include <ogdf/planarlayout/MixedModelLayout.h>
#include <ogdf/planarity/EmbedderMaxFace.h>
#include <ogdf/planarity/SimpleEmbedder.h>
#include <ogdf/planarlayout/TriconnectedShellingOrder.h>
#include <ogdf/planarlayout/BiconnectedShellingOrder.h>

#include "layout_helpers.h"

using namespace ogdf;

go_bandit([](){ bandit::describe("Planar layouts", [](){
	PlanarStraightLayout psl, pslMF, pslBi, pslTri;
	PlanarDrawLayout     pdl, pdlMF, pdlBi, pdlTri;
	MixedModelLayout     mml, mmlMF, mmlBi, mmlTri;

	describeGridLayoutModule("planar straight layout", psl, GR_PLANAR);
	describeGridLayoutModule("planar draw layout", pdl, GR_PLANAR);
	describeGridLayoutModule("mixed mode layout", mml, GR_PLANAR);

	pslMF.setEmbedder(new EmbedderMaxFace);
	pdlMF.setEmbedder(new EmbedderMaxFace);
	mmlMF.setEmbedder(new EmbedderMaxFace);

	describeGridLayoutModule("planar straight layout with embedder max face", pslMF, GR_PLANAR);
	describeGridLayoutModule("planar draw layout with embedder max face", pdlMF, GR_PLANAR);
	describeGridLayoutModule("mixed mode layout with embedder max face", mmlMF, GR_PLANAR);

	pslBi.setEmbedder(new EmbedderMaxFace);
	pslBi.setShellingOrder(new BiconnectedShellingOrder);
	pdlBi.setEmbedder(new EmbedderMaxFace);
	pdlBi.setShellingOrder(new BiconnectedShellingOrder);
	mmlBi.setEmbedder(new SimpleEmbedder);
	mmlBi.setShellingOrder(new BiconnectedShellingOrder);

	describeGridLayoutModule("planar straight layout with embedder max face and bic-order", pslBi, GR_PLANAR);
	describeGridLayoutModule("planar draw layout with embedder max face and bic-order", pdlBi, GR_PLANAR);
	describeGridLayoutModule("mixed mode layout with simple embedder and bic-order", mmlBi, GR_PLANAR);


	pslTri.setEmbedder(new EmbedderMaxFace);
	pslTri.setShellingOrder(new TriconnectedShellingOrder);
	pdlTri.setEmbedder(new EmbedderMaxFace);
	pdlTri.setShellingOrder(new TriconnectedShellingOrder);
	mmlTri.setEmbedder(new EmbedderMaxFace);
	mmlTri.setShellingOrder(new TriconnectedShellingOrder);

	describeGridLayoutModule("planar straight layout with embedder max face and tric-order", pslTri, GR_PLANAR | GR_TRIPLE_CONNECTED);
	describeGridLayoutModule("planar draw layout with embedder max face and tric-order", pdlTri, GR_PLANAR | GR_TRIPLE_CONNECTED);
	describeGridLayoutModule("mixed mode layout with embedder max face and tric-order", mmlTri, GR_PLANAR | GR_TRIPLE_CONNECTED);
}); });
