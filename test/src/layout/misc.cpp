/** \file
 * \brief Tests for several planar layouts
 *
 * \author Tilo Wiedera
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

#include <testing.h>
// IWYU pragma: begin_keep
#include <ogdf/basic/PreprocessorLayout.h>
#include <ogdf/basic/graph_generators/randomized.h>
#include <ogdf/energybased/FMMMLayout.h>
#include <ogdf/misclayout/BalloonLayout.h>
#include <ogdf/misclayout/BertaultLayout.h>
#include <ogdf/misclayout/CircularLayout.h>
#include <ogdf/misclayout/LinearLayout.h>
#include <ogdf/misclayout/ProcrustesSubLayout.h>
#include <ogdf/packing/ComponentSplitterLayout.h>
#include <ogdf/packing/SimpleCCPacker.h>
#include <ogdf/planarlayout/MixedModelLayout.h>
#include <ogdf/tree/RadialTreeLayout.h>
#include <ogdf/tree/TreeLayout.h>
#include <ogdf/upward/DominanceLayout.h>
#include <ogdf/upward/VisibilityLayout.h>

#include <testing.h>
// IWYU pragma: end_keep

#include "layout_helpers.h" // IWYU pragma: associated

static bool edgesHaveBends(const Graph& g, const GraphAttributes& ga) {
	for (edge e : g.edges) {
		if (ga.bends(e).size() > 0) {
			return true;
		}
	}
	return false;
}

go_bandit([] {
	describe("Miscellaneous layouts", [] {
		GraphSizes smallSizes = GraphSizes(16, 32, 16);

		PreprocessorLayout preProc;
		// CircularLayout requires simple graphs
		preProc.setLayoutModule(new CircularLayout);
		describeLayout("PreprocessorLayout with CircularLayout", preProc, 0);

		TEST_LAYOUT(BalloonLayout, GraphProperty::connected);

		describeLayout<BertaultLayout>("BertaultLayout", 0,
				{GraphProperty::sparse, GraphProperty::simple}, false, smallSizes);
		TEST_LAYOUT(CircularLayout, GraphProperty::simple);
		TEST_LAYOUT(LinearLayout);

		ProcrustesSubLayout procrustesLayout(new FMMMLayout);
		describeLayout("ProcrustesSubLayout", procrustesLayout);

		TEST_LAYOUT(ComponentSplitterLayout);

		// BalloonLayout requires connectivity
		SimpleCCPacker packerLayout(new BalloonLayout);
		describeLayout("SimpleCCPacker with BalloonLayout", packerLayout);

		TEST_LAYOUT(RadialTreeLayout, GraphProperty::arborescenceForest, GraphProperty::connected);
		TEST_LAYOUT(TreeLayout, GraphProperty::arborescenceForest);

		// skip iteration with maximum number of nodes as it takes too long
		describeLayout<DominanceLayout>("DominanceLayout", 0,
				{GraphProperty::connected, GraphProperty::simple, GraphProperty::sparse}, false,
				smallSizes);
		describeLayout<VisibilityLayout>("VisibilityLayout", 0,
				{GraphProperty::connected, GraphProperty::simple, GraphProperty::sparse}, false,
				smallSizes);
	});
	describe("ComponentSplitterLayout", [] {
		it("should preserve edge bends of connected component drawings", [&]() {
			Graph graph, graph2;
			GraphAttributes graphAttr(graph), graphAttr2(graph2);

			// generate two planar graphs
			randomPlanarConnectedGraph(graph, 30, 84);
			randomPlanarConnectedGraph(graph2, 30, 84);

			// draw first graph using mixed model that generates bends to edges
			MixedModelLayout mixedModelLayout;
			mixedModelLayout.call(graphAttr);
			// check bends have been set to some edges
			AssertThat(edgesHaveBends(graph, graphAttr), Equals(true));

			// draw second graph using mixed model but wrapped by ComponentSplitterLayout
			ComponentSplitterLayout compSplitterLayout;
			compSplitterLayout.setLayoutModule(new MixedModelLayout);
			compSplitterLayout.call(graphAttr2);
			// edge bends should have been preserved atfer the drawing of connected components
			AssertThat(edgesHaveBends(graph2, graphAttr2), Equals(true));
		});
	});
	describe("SimpleCCPacker", [] {
		it("should preserve edge bends of connected component drawings", [&]() {
			Graph graph, graph2;
			GraphAttributes graphAttr(graph), graphAttr2(graph2);

			// generate two planar graphs
			randomPlanarConnectedGraph(graph, 30, 84);
			randomPlanarConnectedGraph(graph2, 30, 84);

			// add a new connected component in second graph
			node n1 = graph2.newNode();
			node n2 = graph2.newNode();
			node n3 = graph2.newNode();
			graph2.newEdge(n1, n2);
			graph2.newEdge(n2, n3);
			graph2.newEdge(n3, n1);

			// draw first graph using mixed model that generates bends to edges
			MixedModelLayout mixedModelLayout;
			mixedModelLayout.call(graphAttr);
			// check bends have been set to some edges
			AssertThat(edgesHaveBends(graph, graphAttr), Equals(true));

			// draw second graph using mixed model but wrapped by SimpleCCPacker
			SimpleCCPacker simpleCCPacker(new MixedModelLayout);
			simpleCCPacker.call(graphAttr2);
			// edge bends should have been preserved atfer the drawing of connected components
			AssertThat(edgesHaveBends(graph2, graphAttr2), Equals(true));
		});
	});
});
