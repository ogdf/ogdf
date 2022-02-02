/** \file
 * \brief Tests for induced subgraphs.
 *
 * \author Finn Stutzenstein
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
#include <ogdf/basic/graph_generators/deterministic.h>
#include <ogdf/basic/extended_graph_alg.h>

#include <graphs.h>
#include <testing.h>

void assertGraphEqual(const Graph& G, const GraphCopySimple& GC) {
    AssertThat(GC.numberOfNodes(), Equals(G.numberOfNodes()));
    AssertThat(GC.numberOfEdges(), Equals(G.numberOfEdges()));

    for (node n : G.nodes) {
        AssertThat(GC.copy(n), Is().Not().Null());
    }

    for (edge e : G.edges) {
        edge copy = GC.copy(e);
        AssertThat(copy, Is().Not().Null());

        AssertThat(GC.original(copy->source()), Equals(e->source()));
        AssertThat(GC.original(copy->target()), Equals(e->target()));
    }
}

void testFullCopy(const Graph& G) {
    GraphCopySimple GC;
    List<node> nodes;
    for (node n : G.nodes) {
        nodes.pushBack(n);
    }
    inducedSubGraph<ListIterator<node>>(G, nodes.begin(), GC);
    assertGraphEqual(G, GC);
}

go_bandit([] {
	describe("induced subgraph", [] {
        describe("preserves direction", [] {
            it("works forward", [] {
                Graph G;
                customGraph(G, 2, {{0, 1}});
                testFullCopy(G);
            });
            it("works backward", [] {
                Graph G;
                customGraph(G, 2, {{1, 0}});
                testFullCopy(G);
            });
        });
        describe("Can copy full graph", [] {
            forEachGraphItWorks({}, [](Graph& G) {
                testFullCopy(G);
            });
        });
    });
});
