/** \file
 * \brief Tests for node coloring algorithms
 *
 * \author Max Ilsen
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

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/Math.h>
#include <ogdf/graphalg/NodeColoringBergerRompel.h>
#include <ogdf/graphalg/NodeColoringBoppanaHalldorsson.h>
#include <ogdf/graphalg/NodeColoringHalldorsson.h>
#include <ogdf/graphalg/NodeColoringJohnson.h>
#include <ogdf/graphalg/NodeColoringModule.h>
#include <ogdf/graphalg/NodeColoringRecursiveLargestFirst.h>
#include <ogdf/graphalg/NodeColoringSequential.h>
#include <ogdf/graphalg/NodeColoringSimple.h>
#include <ogdf/graphalg/NodeColoringWigderson.h>

#include <functional>
#include <set>
#include <string>
#include <vector>

#include <graphs.h>

#include <testing.h>

using NColor = NodeColoringModule::NodeColor;

static NColor getNumberOfUsedColors(const Graph& graph, const NodeArray<NColor>& colors) {
	NColor numColor {0};
	NColor maxColor {0};
	for (node v : graph.nodes) {
		Math::updateMax(maxColor, colors[v]);
	}
	std::vector<bool> colorUsed(maxColor + 1, false);
	for (node v : graph.nodes) {
		colorUsed[colors[v]] = true;
	}
	for (bool isUsed : colorUsed) {
		if (isUsed) {
			numColor++;
		}
	}
	return numColor;
}

static void describeNodeColoringModule(NodeColoringModule& module, bool isSlow) {
	forEachGraphItWorks(
			{},
			[&](Graph& graph) {
				NodeArray<NColor> colors(graph);
				auto numColors = module.call(graph, colors);
				AssertThat(module.checkColoring(graph, colors), IsTrue());
				AssertThat(numColors == (NColor)0, Equals(graph.numberOfNodes() == 0));
				AssertThat(numColors, Equals(getNumberOfUsedColors(graph, colors)));
			},
			GraphSizes(16, isSlow ? 58 : 100, 42));
}

template<typename ModuleType>
static void describeModule(std::string&& moduleName, bool isSlow = false) {
	describe(moduleName, [&isSlow] {
		ModuleType module;
		describeNodeColoringModule(module, isSlow);
	});
}

go_bandit([]() {
	describeModule<NodeColoringBergerRompel>("NodeColoringBergerRompel");
	describeModule<NodeColoringBoppanaHalldorsson>("NodeColoringBoppanaHalldorsson");
	describeModule<NodeColoringHalldorsson>("NodeColoringHalldorsson", true);
	describeModule<NodeColoringJohnson>("NodeColoringJohnson");
	describeModule<NodeColoringRecursiveLargestFirst>("NodeColoringRecursiveLargestFirst");
	describeModule<NodeColoringSequential>("NodeColoringSequential");
	describeModule<NodeColoringSimple>("NodeColoringSimple");
	describeModule<NodeColoringWigderson>("NodeColoringWigderson");
});
