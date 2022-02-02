/** \file
 * \brief Tests for st-numbering algorithms
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

#include <ogdf/basic/graph_generators.h>
#include <ogdf/basic/STNumbering.h>
#include <testing.h>
#include <graphs.h>

/**
 * Assert that calling computeSTNumbering() on a Graph returns a valid st-numbering.
 *
 * @param G is the Graph to be tested.
 * @param source is the source node for the st-numbering.
 * @param target is the target node for the st-numbering.
 */
void stNumberingAssert(const Graph &G, node source = nullptr, node target = nullptr) {
	NodeArray<int> numbering(G);
	int result = computeSTNumbering(G, numbering, source, target);
	AssertThat(result, Equals(G.numberOfNodes()));
	AssertThat(isSTNumbering(G, numbering, G.numberOfNodes()), IsTrue());

	if (source != nullptr) {
		AssertThat(numbering[source], Equals(1));
	}

	if (target != nullptr) {
		AssertThat(numbering[target], Equals(G.numberOfNodes()));
	}
}


void describeComputeSTNumbering() {
	forEachGraphItWorks({GraphProperty::biconnected, GraphProperty::simple}, [](const Graph &G) {
		stNumberingAssert(G);
		stNumberingAssert(G, G.firstNode());
		stNumberingAssert(G, nullptr, G.firstNode());
		stNumberingAssert(G, G.firstNode(), G.firstNode()->firstAdj()->twinNode());
	}, GraphSizes(), 2);

	it("works on a large simple biconnected graph", []() {
		Graph G;
		randomBiconnectedGraph(G, 300000, 600000);
		makeSimple(G);
		stNumberingAssert(G);
	});
}

go_bandit([]() {
	describe("st-Numbering", []() {
		describe("computeSTNumbering", []() {
			describeComputeSTNumbering();
		});
	});
});
