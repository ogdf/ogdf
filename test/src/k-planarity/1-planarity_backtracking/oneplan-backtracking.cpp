/** \file
 * \brief Tests for OnePlanarityBacktracking.
 *
 * \author Matthias Pfretzschner
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
#include <ogdf/basic/GraphSets.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/Module.h>
#include <ogdf/basic/graph_generators/deterministic.h>
#include <ogdf/k-planarity/1-planarity_backtracking/EdgePairPartition.h>
#include <ogdf/k-planarity/1-planarity_backtracking/OnePlanarityBacktracking.h>
#include <ogdf/k-planarity/1-planarity_backtracking/OnePlanarization.h>

#include <functional>
#include <set>
#include <string>
#include <vector>

#include <graphs.h>
#include <resources.h>

#include <testing.h>

using namespace oneplan_backtracking;

void testGraph(const Graph& G, OneplanMode mode, bool yesInstance) {
	OnePlanarityBacktracking solver;
	OnePlanarization opl;
	Module::ReturnType res;
	switch (mode) {
	case OneplanMode::IC:
		res = solver.testICPlanarity(G, &opl);
		break;
	case OneplanMode::NIC:
		res = solver.testNICPlanarity(G, &opl);
		break;
	case OneplanMode::Normal:
	default:
		res = solver.testOnePlanarity(G, &opl);
		break;
	}

	if (yesInstance) {
		AssertThat(res, Equals(Module::ReturnType::Feasible));
		AssertThat(opl.isPlanar(), IsTrue());
		AssertThat(opl.numberOfNodes(), Equals(G.numberOfNodes() + opl.crossingVertices().size()));

		EdgeSet mappedEdges(G);
		for (edge e : opl.crossingEdges()) {
			AssertThat(opl.original(e), !IsNull());
			mappedEdges.insert(opl.original(e));
		}
		for (edge e : opl.freeEdges()) {
			AssertThat(opl.original(e), !IsNull());
			mappedEdges.insert(opl.original(e));
		}

		for (edge e : opl.kiteEdges()) {
			AssertThat(opl.original(e), IsNull());
		}

		for (edge e : opl.remainingEdges()) {
			AssertThat(opl.original(e), !IsNull());
			mappedEdges.insert(opl.original(e));
		}
		AssertThat(mappedEdges.size(), Equals(G.numberOfEdges()));
	} else {
		AssertThat(res, Equals(Module::ReturnType::NoFeasibleSolution));
	}
}

go_bandit([]() {
	describe("1-Planarity test", []() {
		it("works for planar graphs", []() {
			forEachGraphItWorks({GraphProperty::planar, GraphProperty::simple}, [&](Graph& G) {
				testGraph(G, OneplanMode::Normal, true);
				testGraph(G, OneplanMode::NIC, true);
				testGraph(G, OneplanMode::IC, true);
			});
		});
		it("works for trivially non-planar graphs", []() {
			Graph G;
			completeGraph(G, 10);
			testGraph(G, OneplanMode::Normal, false);
			testGraph(G, OneplanMode::NIC, false);
			testGraph(G, OneplanMode::IC, false);
		});

		for_each_graph_it("works", {"north/g.41.26.gml", "north/g.73.8.gml", "north/g.18.1.gml"},
				[&](Graph& G) {
					testGraph(G, OneplanMode::Normal, true);
					testGraph(G, OneplanMode::NIC, true);
					testGraph(G, OneplanMode::IC, false);
				});

		for_each_graph_it("works", {"north/g.13.94.gml"}, [&](Graph& G) {
			testGraph(G, OneplanMode::Normal, false);
			testGraph(G, OneplanMode::NIC, false);
			testGraph(G, OneplanMode::IC, false);
		});

		for_each_graph_it("works", {"north/g.23.108.gml"}, [&](Graph& G) {
			testGraph(G, OneplanMode::Normal, true);
			testGraph(G, OneplanMode::NIC, true);
			testGraph(G, OneplanMode::IC, true);
		});

		for_each_graph_it("works", {"north/g.10.27.gml"}, [&](Graph& G) {
			testGraph(G, OneplanMode::Normal, true);
			testGraph(G, OneplanMode::NIC, false);
			testGraph(G, OneplanMode::IC, false);
		});

		for_each_graph_it("works with fewer threads", {"north/g.10.27.gml"}, [&](Graph& G) {
			OnePlanarityBacktracking solver(1);
			AssertThat(solver.testOnePlanarity(G), Equals(Module::ReturnType::Feasible));
			AssertThat(solver.testNICPlanarity(G), Equals(Module::ReturnType::NoFeasibleSolution));
			AssertThat(solver.testICPlanarity(G), Equals(Module::ReturnType::NoFeasibleSolution));
		});
	});
});
