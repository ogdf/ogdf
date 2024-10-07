/** \file
 * \brief Applies the node coloring approximation specified by Boppana&Halldorsson.
 *
 * \author Jan-Niklas Buckow
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
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/basic.h>
#include <ogdf/graphalg/NodeColoringBoppanaHalldorsson.h>
#include <ogdf/graphalg/NodeColoringModule.h>

namespace ogdf {

using NColor = NodeColoringModule::NodeColor;

NColor NodeColoringBoppanaHalldorsson::call(const Graph& graph, NodeArray<NColor>& colors,
		NColor start) {
	NColor numberOfColorsUsed = 0;
	// Copy the input graph
	GraphCopy graphMain = GraphCopy(graph);
	preprocessGraph(graphMain);
	NodeArray<NColor> colorsMain(graphMain);

	// Color each independent set until the graph is colored
	while (!graphMain.empty()) {
		List<node> ramseyClique;
		List<node> ramseyIndependentSet;
		ramseyAlgorithm(graphMain, ramseyClique, ramseyIndependentSet);

		for (node& v : ramseyIndependentSet) {
			colors[graphMain.original(v)] = start;
			graphMain.delNode(v);
		}

		start++;
		numberOfColorsUsed++;
	}

	OGDF_ASSERT(checkColoring(graph, colors));
	return numberOfColorsUsed;
}

}
