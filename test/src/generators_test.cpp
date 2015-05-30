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
 * See README.txt in the root directory of the OGDF installation for details.
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
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/

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
			it(string("generates a graph with " + to_string(n) + " nodes and " + to_string(m) + " edges").c_str(), [&](){
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
			it(string("generates a graph with " + to_string(n) + " nodes and " + to_string(m) + " edges").c_str(), [&](){
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
			it(string("generates a graph with " + to_string(n) + " nodes and " + to_string(m) + " edges").c_str(), [&](){
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
			it(string("generates a graph with " + to_string(n) + " nodes").c_str(), [&](){
				Graph G;
				randomTriconnectedGraph(G, n, .5, .5);
				AssertThat(G.numberOfNodes(), Equals(n));
				AssertThat(isTriconnected(G), Equals(true));
			});
		}
	});

	describe("randomTree", [](){
		for(int n = 1; n < 100; n++) {
			it(string("generates a graph with " + to_string(n) + " nodes").c_str(), [&](){
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
			it(string("generates a graph with " + to_string(n) + " nodes and " + to_string(m) + " edges").c_str(), [&](){
				Graph G;
				randomHierarchy(G, n, m, false, false, true);
				AssertThat(G.numberOfNodes(), Equals(n));
				AssertThat(G.numberOfEdges(), Equals(m));
			});
		}
	});

	describe("randomDiGraph", [](){
		for(int n = 1; n < 100; n++) {
			it(string("generates a graph with " + to_string(n) + " nodes").c_str(), [&](){
				Graph G;
				randomDiGraph(G, n, .5);
				AssertThat(G.numberOfNodes(), Equals(n));
				AssertThat(isSimple(G), Equals(true));
			});
		}
	});

});
});
