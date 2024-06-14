#pragma once
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/NodeSet.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/layered/FastHierarchyLayout.h>
#include <ogdf/layered/SugiyamaLayout.h>

#include <vector>

using namespace ogdf;
using namespace std;

void writeLevelGraph(const Graph& G, const vector<vector<node>>& emb, const NodeArray<int>& pos,
		ostream& os);

void readLevelGraph(Graph& G, vector<vector<node>>& emb, NodeArray<int>& pos, istream& is);

void layout(const Graph& G, const NodeArray<int>& lvl, const NodeArray<int>& pos,
		GraphAttributes& GA);

void drawSVG(const Graph& G, const NodeArray<int>& lvl, const NodeArray<int>& pos,
		const string& name);

struct AdjCompLess {
	const NodeArray<int>& pos;

	explicit AdjCompLess(const NodeArray<int>& pos) : pos(pos) { }

	bool operator()(adjEntry a, adjEntry b);
};

void embedPLE(Graph& G, const vector<vector<node>>& emb, NodeArray<int>& lvl, NodeArray<int>& pos);

void checkPLE(const Graph& G, const vector<vector<node>>& emb, const NodeArray<int>& lvl,
		const NodeArray<int>& pos);

void randomProperMaximalLevelPlaneGraph(Graph& G, vector<vector<node>>& emb, int N, int K,
		bool radial);

void pruneEdges(Graph& G, int max_edges, int min_deg);

void reduceLevelToCluster(const Graph& LG, const vector<vector<node>>& emb, Graph& G,
		ClusterGraph& CG);
