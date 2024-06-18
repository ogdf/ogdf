/** \file
 * \brief TODO Document
 *
 * \author Simon D. Fink <ogdf@niko.fink.bayern>
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
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/basic.h>
#include <ogdf/cluster/ClusterGraph.h>
#include <ogdf/cluster/sync_plan/utils/Levels.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/layered/FastHierarchyLayout.h>
#include <ogdf/layered/Hierarchy.h>
#include <ogdf/layered/HierarchyLevels.h>

#include <algorithm>
#include <istream>
#include <random>
#include <string>
#include <utility>
#include <vector>

using namespace std;

namespace ogdf::sync_plan {

void writeLevelGraph(const Graph& G, const vector<vector<node>>& emb, const NodeArray<int>& pos,
		ostream& os) {
	os << G.numberOfNodes() << " " << G.numberOfEdges() << " " << emb.size() << "\n";
	for (int l = 0; l < emb.size(); ++l) {
		os << emb[l].size() << " ";
	}
	os << "\n";

	for (int l = 0; l + 1 < emb.size(); ++l) {
		for (int i = 0; i < emb[l].size(); ++i) {
			if (emb[l][i]->outdeg() == 0) {
				continue;
			}
			os << l << " " << i << " : ";
			for (adjEntry a : emb[l][i]->adjEntries) {
				if (a->isSource()) {
					os << pos[a->twinNode()] << " ";
				}
			}
			os << ";\n";
		}
	}
}

void readLevelGraph(Graph& G, vector<vector<node>>& emb, NodeArray<int>& pos, istream& is) {
	int N, M, L;
	is >> N >> M >> L;
	emb.resize(L);
	for (int l = 0; l < L; ++l) {
		int c;
		is >> c;
		for (int i = 0; i < c; ++i) {
			node n = G.newNode();
			pos[n] = i;
			emb[l].push_back(n);
		}
	}
	while (G.numberOfEdges() < M) {
		int l, lp, up;
		char c;
		is >> l >> lp >> c;
		OGDF_ASSERT(c == ':');
		node n = emb[l][lp];
		while (is >> c && c != ';') {
			is.putback(c);
			is >> up;
			G.newEdge(n, emb[l + 1][up]);
		}
	}
	OGDF_ASSERT(G.numberOfNodes() == N);
	OGDF_ASSERT(G.numberOfEdges() == M);
}

void layout(const Graph& G, const NodeArray<int>& lvl, const NodeArray<int>& pos,
		GraphAttributes& GA) {
	Hierarchy H {G, lvl};
	HierarchyLevels HL {H};
	const GraphCopy& GC = (const GraphCopy&)H;
	NodeArray<int> HLpos {GC, -1};
	for (node v : GC.nodes) {
		HLpos[v] = pos[GC.original(v)];
	}
	HL.restorePos(HLpos);
	FastHierarchyLayout L;
	L.call(HL, GA);
	GA.flipVertical();
}

void drawSVG(const Graph& G, const NodeArray<int>& lvl, const NodeArray<int>& pos,
		const string& name) {
	GraphAttributes GA {G, GraphAttributes::all};
	for (node n : G.nodes) {
		GA.label(n) = to_string(n->index());
	}
	layout(G, lvl, pos, GA);
	ogdf::GraphIO::write(GA, name);
}

void embedPLE(Graph& G, const vector<vector<node>>& emb, NodeArray<int>& lvl, NodeArray<int>& pos) {
	vector<adjEntry> adjs;
	AdjCompLess less(pos);
	for (int i = 0; i < emb.size(); ++i) {
		for (int j = 0; j < emb[i].size(); ++j) {
			node n = emb[i][j];
			lvl[n] = i;
			pos[n] = j;

			adjs.clear();
			adjs.reserve(n->degree());
			for (adjEntry adj : n->adjEntries) {
				adjs.push_back(adj);
			}
			sort(adjs.begin(), adjs.end(), less);
			// G.sort(n, adjs);
		}
	}
}

void checkPLE(const Graph& G, const vector<vector<node>>& emb, const NodeArray<int>& lvl,
		const NodeArray<int>& pos) {
	// check that the arrays are consistent
	int c = 0;
	for (int i = 0; i < emb.size(); ++i) {
		for (int j = 0; j < emb[i].size(); ++j) {
			OGDF_ASSERT(lvl[emb[i][j]] == i);
			OGDF_ASSERT(pos[emb[i][j]] == j);
			c++;
		}
	}
	OGDF_ASSERT(c == G.numberOfNodes());

	// check that all edges are proper
	for (edge e : G.edges) {
		OGDF_ASSERT(lvl[e->source()] + 1 == lvl[e->target()]);
	}

	// check that the embedding is bimodal
	for (node n : G.nodes) {
		if (n->degree() < 2) {
			continue;
		}
		bool up = false;
		adjEntry pred = n->lastAdj(), curr = n->firstAdj();
		while (curr != nullptr) {
			if (curr->isSource()) {
				up = true;
				OGDF_ASSERT(lvl[curr->twinNode()] == lvl[n] + 1);
			} else {
				OGDF_ASSERT(lvl[curr->twinNode()] == lvl[n] - 1);
				OGDF_ASSERT(!up);
			}
			if (curr != n->firstAdj()) {
				if (curr->isSource() && pred->isSource()) {
					OGDF_ASSERT(pos[curr->twinNode()] >= pos[pred->twinNode()]);
				} else if (!curr->isSource() && !pred->isSource()) {
					OGDF_ASSERT(pos[curr->twinNode()] <= pos[pred->twinNode()]);
				}
			}
			pred = curr;
			curr = curr->succ();
		}
	}

	// check that the embedding is planar
	for (int l = 0; l < emb.size(); l++) {
		int leftmost = -1;
		for (int i = 0; i < emb[l].size(); ++i) {
			if (emb[l][i]->outdeg() < 1) {
				continue;
			}
			for (adjEntry adj : emb[l][i]->adjEntries) {
				if (!adj->isSource()) {
					continue;
				}
				int p = pos[adj->twinNode()];
				OGDF_ASSERT(p >= leftmost);
				leftmost = max(p, leftmost);
			}
		}
	}
}

void randomProperMaximalLevelPlaneGraph(Graph& G, vector<vector<node>>& emb, int N, int K,
		bool radial) {
	// checks
	OGDF_ASSERT(N > 0);
	OGDF_ASSERT(K > 0);

	// init arrays
	G.clear();
	emb.clear();
	emb.resize(K);

	// create nodes, at least one on each level
	for (int i = 0; i < K; ++i) {
		emb[i].push_back(G.newNode());
	}
	for (int i = K; i < N; ++i) {
		emb[randomNumber(0, K - 1)].push_back(G.newNode());
	}

	// create all edges on each level
	for (int l = 0; l < K - 1; ++l) {
		int lp = 0;
		int up = 0;
		bool wrap_forw = false;
		bool wrap_back = false;
		while (true) {
			G.newEdge(emb[l][lp], emb[l + 1][up]);

			wrap_forw = wrap_forw || (lp + 1 == emb[l].size() && up == 0);
			wrap_back = wrap_back || (lp == 0 && up + 1 == emb[l + 1].size());

			// select the next node either on the upper or the lower level
			if (lp + 1 < emb[l].size() && up + 1 < emb[l + 1].size()) {
				int r = randomNumber(0, 1);
				lp += r;
				up += 1 - r;
			} else if (lp + 1 < emb[l].size()) {
				lp++;
			} else if (up + 1 < emb[l + 1].size()) {
				up++;
			} else {
				break;
			}
		}

		// add the wrap-around edge
		if (radial) {
			if (randomNumber(0, 1) == 1 || wrap_back) {
				if (!wrap_forw) {
					G.newEdge(emb[l].back(), emb[l + 1].front());
				}
			} else {
				G.newEdge(emb[l].front(), emb[l + 1].back());
			}
		}
	}
}

void pruneEdges(Graph& G, int max_edges, int min_deg) {
	vector<edge> edges;
	for (edge e : G.edges) {
		edges.push_back(e);
	}
	mt19937 mt(randomSeed());
	shuffle(edges.begin(), edges.end(), mt);
	for (edge e : edges) {
		if (e->source()->degree() > min_deg && e->target()->degree() > min_deg) {
			G.delEdge(e);
		}
		if (G.numberOfEdges() <= max_edges) {
			break;
		}
	}
}

bool AdjCompLess::operator()(adjEntry a, adjEntry b) {
	if (a->theNode() != b->theNode()) {
		return a->theNode()->index() < b->theNode()->index();
	}
	if (a->isSource()) {
		if (b->isSource()) {
			return pos[a->twinNode()] < pos[b->twinNode()];
		} else {
			return false;
		}
	} else {
		if (b->isSource()) {
			return true;
		} else {
			return pos[a->twinNode()] > pos[b->twinNode()];
		}
	}
}

void reduceLevelToCluster(const Graph& LG, const vector<vector<node>>& emb, Graph& G,
		ClusterGraph& CG) {
	NodeArray<pair<node, node>> map(LG);
	cluster p = CG.rootCluster();
	for (int l = emb.size() - 1; l >= 0; --l) {
		cluster c = CG.newCluster(p);
		for (int i = 0; i < emb[l].size(); ++i) {
			node n = emb[l][i];
			node u = G.newNode();
			node v = G.newNode();
			CG.reassignNode(u, c);
			CG.reassignNode(v, p);
			map[n] = {u, v};
			G.newEdge(u, v);
		}
		p = c;
	}
	for (edge e : LG.edges) {
		G.newEdge(map[e->source()].second, map[e->target()].first);
	}
}

}
