/** \file
 * \brief Simple tests for generating various graphs.
 *
 * \author Christoph Schulz, Tilo Wiedera
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

#include "bandit/bandit.h"
#include "ogdf/basic/Graph.h"
#include "ogdf/basic/graph_generators.h"
#include "ogdf/basic/simple_graph_alg.h"

using namespace ogdf;
using namespace bandit;

// TODO: Test overloaded functions

go_bandit([](){

describe("Special graph classes", [](){

    it("generates two circulant graphs",[](){
        Graph G;
        circulantGraph(G, 11, Array<int>({1,2,4}));
        AssertThat(G.numberOfEdges(), Equals(33));
        AssertThat(G.numberOfNodes(), Equals(11));
        AssertThat(isConnected(G),Equals(true));

        circulantGraph(G, 12, Array<int>({2,4,6}));
        AssertThat(G.numberOfNodes(), Equals(12));
        AssertThat(isConnected(G),Equals(false));
    });

});

describe("Random generators", [](){

	describe("randomGraph", [](){
		for(int n = 0; n < 100; n++) {
			int m = randomNumber(0, (n*(n-1))/2);
			it(string("generates a graph with " + to_string(n) + " nodes and " + to_string(m) + " edges"), [&](){
				Graph G;
				randomGraph(G, n, m);
				AssertThat(G.numberOfNodes(), Equals(n));
				AssertThat(G.numberOfEdges(), Equals(m));
			});
		}
	});

	describe("randomSimpleGraph", [](){
		for(int n = 0; n < 100; n++) {
			int m = randomNumber(0, (n*(n-1))/2);
			it(string("generates a graph with " + to_string(n) + " nodes and " + to_string(m) + " edges"), [&](){
				Graph G;
				randomSimpleGraph(G, n, m);
				AssertThat(G.numberOfNodes(), Equals(n));
				AssertThat(G.numberOfEdges(), Equals(m));
				AssertThat(isSimple(G), Equals(true));
			});
		}
	});

	describe("randomBiconnectedGraph", [](){
		for(int n = 3; n < 100; n++) {
			int m = randomNumber(n, (n*(n-1))/2);
			it(string("generates a graph with " + to_string(n) + " nodes and " + to_string(m) + " edges"), [&](){
				Graph G;
				randomBiconnectedGraph(G, n, m);
				AssertThat(G.numberOfNodes(), Equals(n));
				AssertThat(G.numberOfEdges(), Equals(m));
				AssertThat(isBiconnected(G), Equals(true));
			});
		}
	});

	describe("randomTriconnectedGraph", [](){
		for(int n = 4; n < 100; n++) {
			it(string("generates a graph with " + to_string(n) + " nodes"), [&](){
				Graph G;
				randomTriconnectedGraph(G, n, .5, .5);
				AssertThat(G.numberOfNodes(), Equals(n));
				AssertThat(isTriconnected(G), Equals(true));
			});
		}
	});

	describe("randomTree", [](){
		for(int n = 1; n < 100; n++) {
			it(string("generates a graph with " + to_string(n) + " nodes"), [&](){
				Graph G;
				randomTree(G, n);
				AssertThat(G.numberOfNodes(), Equals(n));
				AssertThat(isTree(G), Equals(true));
			});
		}
	});

	// TODO: dont skip me
	describe_skip("randomHierarchy", [](){
		for(int n = 1; n < 100; n++) {
			int m = randomNumber(n-1, (n*(n-1))/2);
			it(string("generates a graph with " + to_string(n) + " nodes and " + to_string(m) + " edges"), [&](){
				Graph G;
				randomHierarchy(G, n, m, false, false, true);
				AssertThat(G.numberOfNodes(), Equals(n));
				AssertThat(G.numberOfEdges(), Equals(m));
			});
		}
	});

	describe("randomDiGraph", [](){
		for(int n = 1; n < 100; n++) {
			it(string("generates a graph with " + to_string(n) + " nodes"), [&](){
				Graph G;
				randomDiGraph(G, n, .5);
				AssertThat(G.numberOfNodes(), Equals(n));
				AssertThat(isSimple(G), Equals(true));
			});
		}
	});

	describe("randomRegularGraph", []() {
		for (int n = 10; n <= 30; n += 5) {
			for (int d = 2; d <= 6; d += 2) {
				it(string("generates a graph with degree " + to_string(d) + " and " + to_string(n) + " nodes"), [&]() {
					Graph G;
					randomRegularGraph(G, n, d);
					AssertThat(G.numberOfNodes(), Equals(n));
					AssertThat(isSimple(G), Equals(true));
					AssertThat(isRegular(G, d), Equals(true));
				});
			}
		}
	});

	describe("randomGeometricCubeGraph", [](){
		for(int d = 1; d < 4; d++){
			for(double t : {0.0, 0.1, 0.5}) {
				for(int n = 0; n < 100; n++) {
					it(string("generates a graph with " + to_string(n) +
							  " nodes in dim " + to_string(d) +
							  " and threshold " + to_string(t)), [&](){
						Graph G;
						randomGeometricCubeGraph(G,n,t,d);
						AssertThat(G.numberOfNodes(), Equals(n));
						AssertThat(isSimple(G), Equals(true));
					});
				}
			}
		}
	});

	describe("emptyGraph", [](){
		for(int n = 0; n < 100; n++) {
			it(string("generates a graph with " + to_string(n) + " nodes"), [&](){
				Graph G;
				emptyGraph(G, n);
				AssertThat(G.numberOfNodes(), Equals(n));
				AssertThat(G.numberOfEdges(), Equals(0));
			});
		}
	});
});
});
