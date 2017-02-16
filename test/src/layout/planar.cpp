/** \file
 * \brief Tests for several planar layouts
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

#include <ogdf/planarlayout/PlanarStraightLayout.h>
#include <ogdf/planarlayout/PlanarDrawLayout.h>
#include <ogdf/planarlayout/MixedModelLayout.h>
#include <ogdf/planarity/EmbedderMaxFace.h>
#include <ogdf/planarity/EmbedderMaxFaceLayers.h>
#include <ogdf/planarity/EmbedderMinDepth.h>
#include <ogdf/planarity/EmbedderMinDepthMaxFace.h>
#include <ogdf/planarity/EmbedderMinDepthMaxFaceLayers.h>
#include <ogdf/planarity/EmbedderMinDepthPiTa.h>
#include <ogdf/planarity/EmbedderOptimalFlexDraw.h>
#include <ogdf/planarity/SimpleEmbedder.h>
#include <ogdf/planarlayout/TriconnectedShellingOrder.h>
#include <ogdf/planarlayout/BiconnectedShellingOrder.h>

#include "layout_helpers.h"

using namespace ogdf;

template<typename Layout>
static void describeForAllEmbedders(Layout &layout, string name, std::initializer_list<GraphRequirement> requirements) {
	name += " and ";

	layout.setEmbedder(new SimpleEmbedder);
	describeGridLayoutModule(name + " SimpleEmbedder", layout, requirements);

	layout.setEmbedder(new EmbedderMaxFace);
	describeGridLayoutModule(name + " EmbedderMaxFace", layout, requirements);

	layout.setEmbedder(new EmbedderMaxFaceLayers);
	describeGridLayoutModule(name + "EmbedderMaxFaceLayers", layout, requirements);

	layout.setEmbedder(new EmbedderMinDepth);
	describeGridLayoutModule(name + "EmbedderMinDepth", layout, requirements);

	layout.setEmbedder(new EmbedderMinDepthMaxFace);
	describeGridLayoutModule(name + "EmbedderMinDepthMaxFace", layout, requirements);

	layout.setEmbedder(new EmbedderMinDepthMaxFaceLayers);
	describeGridLayoutModule(name + "EmbedderMinDepthMaxFaceLayers", layout, requirements);

	// the two embedders below are currently disabled since they cause failures

	layout.setEmbedder(new EmbedderMinDepthPiTa);
	describeGridLayoutModule(name + "EmbedderMinDepthPiTa", layout, requirements, MAX_NODES, true);

	layout.setEmbedder(new EmbedderOptimalFlexDraw);
	describeGridLayoutModule(name + "EmbedderOptimalFlexDraw", layout, requirements, MAX_NODES, true);
}

template<typename Layout>
static void describePlanarLayout(string name) {
	Layout layout;
	name += " with ";

	layout.setShellingOrder(new BiconnectedShellingOrder);
	describeForAllEmbedders(layout, name + "BiconnectedShellingOrder", {GraphRequirement::planar});

	layout.setShellingOrder(new TriconnectedShellingOrder);
	describeForAllEmbedders(layout, name + "TriconnectedShellingOrder", {GraphRequirement::planar, GraphRequirement::triconnected});
}

go_bandit([] { bandit::describe("Planar layouts", [] {
	describePlanarLayout<PlanarStraightLayout>("PlanarStraightLayout");
	describePlanarLayout<PlanarDrawLayout>("PlanarDrawLayout");
	describePlanarLayout<MixedModelLayout>("MixedModelLayout");
}); });
