//**************************************************************
//  regression test for LCA class
//
//  Author: Stephan Beyer
//**************************************************************

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/List.h>
#include <ogdf/tree/LCA.h>
#include <ogdf/basic/graph_generators.h>

using namespace ogdf;

// get ancestor in tree
static node
ancestor(node v) {
	edge e;
	forall_adj_edges(e, v) {
		if (e->target() == v) {
			return e->source();
		}
	}
	return NULL;
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

static bool
testTree(int q, int n, __int64 &genTime, __int64 &evalTime, __int64 &checkTime, int maxDeg = 0, int maxWidth = 0)
{
	Graph G;
	randomTree(G, n, maxDeg, maxWidth);
	node root = G.firstNode();

	__int64 time;
	System::usedRealTime(time);
	LCA lca(G, root);
	genTime += System::usedRealTime(time);

	while (q--) {
		node u = G.chooseNode();
		node v = G.chooseNode();
		System::usedRealTime(time);
		node a = lca.call(u, v);
		evalTime += System::usedRealTime(time);
		bool check = checkLCA(a, v, u);
		checkTime += System::usedRealTime(time);
		if (!check) {
			return false;
		}
	}
	return true;
}

static bool
testTrees(const int q, const int n, int nodes, int maxDeg = 0, int maxWidth = 0)
{
	Graph G;
	__int64 genTime(0), evalTime(0), checkTime(0), time;

	cout << "-> " << q << " LCA queries on " << n << " random trees with " << nodes << " nodes";
	if (maxDeg > 0) {
		cout << ", maximum degree " << maxDeg;
	}
	if (maxWidth > 0) {
		cout << " and maximum width " << maxWidth;
	}
	cout << "\n";
	int i = n;
	System::usedRealTime(time);
	while (i--) {
		if (!testTree(q, nodes, genTime, evalTime, checkTime, maxDeg, maxWidth)) {
			return false;
		}
	}
	cout
	 << "    avg LCA build time in ms: " << double(genTime) / n << "\n"
	    "    avg LCA query time in ms: " << double(evalTime) / n / q << "\n"
	    "    avg LCA check time in ms: " << double(checkTime) / n / q << " (naive)\n"
	    "            total time in  s: " << 0.001 * double(System::usedRealTime(time)) << "\n";
	return true;
}

bool
regLCA()
{
	return testTrees(3000, 1000, 100)
	    && testTrees(3000, 1000, 100,  33)
	    && testTrees(3000, 1000, 100,  2)
	    && testTrees(3000, 1000, 100,  1)
	    && testTrees(3000, 1000, 100,  0, 10)
	    && testTrees(3000, 200,  1000, 0, 15)
	    && testTrees(3000, 200,  1000, 2)
	    && testTrees(3000, 25,   10000)
	    && testTrees(1500, 5,    100000);
}
