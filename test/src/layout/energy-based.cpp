/** \file
 * \brief Tests for several energy-based layout classes
 *
 * \author Carsten Gutwenger, Tilo Wiedera
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.md in the OGDF root directory for details.
 *
 * \par
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * Version 2 or 3 as published by the Free Software Foundation;
 * see the file LICENSE.txt included in the packaging of this file
 * for details.
 *
 * \par
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * \par
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, see
 * http://www.gnu.org/copyleft/gpl.html
 */

#include <bandit/bandit.h>

#include <ogdf/energybased/FMMMLayout.h>
#include <ogdf/energybased/SpringEmbedderGridVariant.h>
#include <ogdf/energybased/GEMLayout.h>
#include <ogdf/energybased/DavidsonHarelLayout.h>
#include <ogdf/energybased/PivotMDS.h>

#include "layout_helpers.h"

using namespace ogdf;

go_bandit([](){ bandit::describe("Energy-based layouts", [](){
	FMMMLayout                fmmm, fmmmNice, fmmmHQ;
	SpringEmbedderGridVariant frl, frlHQ;
	GEMLayout                 gem;
	DavidsonHarelLayout       dhl;
	PivotMDS                  pmds;

	fmmmHQ.qualityVersusSpeed(FMMMOptions::QualityVsSpeed::GorgeousAndEfficient);
	fmmmHQ.useHighLevelOptions(true);
	fmmmNice.qualityVersusSpeed(FMMMOptions::QualityVsSpeed::NiceAndIncredibleSpeed);
	frlHQ.iterations(1000);

	// call              name                                                                         module    extraAttributes req        maxNodes isGridLayout skip
	describeLayoutModule("Fast Multipole Multilevel Embedder"                                       , fmmm);
	describeLayoutModule("Fast Multipole Multilevel Embedder with high quality settings"            , fmmmHQ);
	describeLayoutModule("Fast Multipole Multilevel Embedder with nice quality and incredible speed", fmmmNice);
	describeLayoutModule("Spring Embedder grid variant"                                             , frl     , 0              , {}      , 200    , false);
	describeLayoutModule("Spring Embedder grid variant with high quality settings"                  , frlHQ   , 0              , {}      , 200    , false);
	describeLayoutModule("GEM layout"                                                               , gem);
	describeLayoutModule("Davidson-Harel layout"                                                    , dhl);
	describeLayoutModule("PivotMDS layout"                                                          , pmds    , 0              , {GraphRequirement::connected});
}); });
