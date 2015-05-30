/** \file
 * \brief Tests for LCA (lowest common ancestor) class
 *
 * \author Stephan Beyer
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

#include <sstream>
#include <bandit/bandit.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/List.h>
#include <ogdf/tree/LCA.h>
#include <ogdf/basic/graph_generators.h>

using namespace ogdf;
using namespace bandit;

// get ancestor in tree
static node
ancestor(node v)
{
	edge e;
	forall_adj_edges(e, v) {
		if (e->target() == v) {
			return e->source();
		}
	}
	return nullptr;
}

// check if a is LCA of v and u
static bool
checkLCA(node a, node v, node u)
{
	if (a == v && v == u) {
		return true;
	}
	node vp = v, up = u; // initialization needed for the case that v == a or u == a
	while (v && v != a) {
		vp = v;
		v = ancestor(v);
	}
	while (u && u != a) {
		up = u;
		u = ancestor(u);
	}
	return (v == a && u == a && vp != up);
}

static void
testTree(int q, int n, int64_t &genTime, int64_t &evalTime, int64_t &checkTime, int maxDeg = 0, int maxWidth = 0)
{
	Graph G;
	randomTree(G, n, maxDeg, maxWidth);
	node root = G.firstNode();

	int64_t time;
	System::usedRealTime(time);
	LCA lca(G, root);
	genTime += System::usedRealTime(time);

	while (q--) {
		node u = G.chooseNode();
		node v = G.chooseNode();
		System::usedRealTime(time);
		node a = lca.call(u, v);
		evalTime += System::usedRealTime(time);
		bool isValidQueryResult = checkLCA(a, v, u);
		checkTime += System::usedRealTime(time);
		AssertThat(isValidQueryResult, Equals(true));
	}
}

static void
lotsOfQueries(const int q, const int n, int nodes, int maxDeg = 0, int maxWidth = 0)
{
	std::stringstream ss;
	ss << "makes " << q << " LCA queries on " << n << " random trees with " << nodes << " nodes";
	if (maxDeg > 0) {
		ss << ", maximum degree " << maxDeg;
	}
	if (maxWidth > 0) {
		ss << " and maximum width " << maxWidth;
	}
	it(ss.str().c_str(), [&]() {
		Graph G;
		int64_t genTime(0), evalTime(0), checkTime(0), time;

		int i = n;
		System::usedRealTime(time);
		while (i--) {
			testTree(q, nodes, genTime, evalTime, checkTime, maxDeg, maxWidth);
		}
		cout << "\n"
		  "      avg LCA build time in ms: " << double(genTime) / n << "\n"
		  "      avg LCA query time in ms: " << double(evalTime) / n / q << "\n"
		  "      avg LCA check time in ms: " << double(checkTime) / n / q << " (naive)\n"
		  "              total time in  s: " << 0.001 * double(System::usedRealTime(time)) << "\n";
	});
}

static void
fewNodeTests()
{
	Graph G;
	node root = G.newNode();
	it_skip("constructs LCA data structure on a tree with one node", [&]() {
		LCA lca(G, root);
	});
	it_skip("makes LCA query on tree with one node", [&]() {
		LCA lca(G, root);
		node commonAncestor = lca.call(root, root);
		AssertThat(commonAncestor, Equals(root));
	});
	node vl = G.newNode();
	G.newEdge(root, vl);
	it("constructs LCA data structure on a tree with two nodes", [&]() {
		LCA lca(G, root);
	});
	it("makes LCA queries on a tree with two nodes", [&]() {
		LCA lca(G, root);
		AssertThat(lca.call(root, root), Equals(root));
		AssertThat(lca.call(vl, root), Equals(root));
		AssertThat(lca.call(root, vl), Equals(root));
		AssertThat(lca.call(vl, vl), Equals(vl));
	});
	node vr = G.newNode();
	G.newEdge(root, vr);
	it("makes LCA queries on a tree with three nodes", [&]() {
		LCA lca(G, root);
		AssertThat(lca.call(vr, root), Equals(root));
		AssertThat(lca.call(vl, vr), Equals(root));
		AssertThat(lca.call(vr, vr), Equals(vr));
	});
	return;
}

go_bandit([]() {
	describe("Lowest Common Ancestor algorithm", []() {
		fewNodeTests();
		describe("lots of LCA queries", []() {
			lotsOfQueries(3000, 1000,    100);
			lotsOfQueries(3000, 1000,    100, 33);
			lotsOfQueries(3000, 1000,    100,  2);
			lotsOfQueries(3000, 1000,    100,  1);
			lotsOfQueries(3000, 1000,    100,  0, 10);
			lotsOfQueries(3000,  200,   1000,  0, 15);
			lotsOfQueries(3000,  200,   1000,  2);
			lotsOfQueries(3000,   25,  10000);
			lotsOfQueries(1500,    5, 100000);
		});
	});
});
