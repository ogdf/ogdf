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
#include <ogdf/basic/CombinatorialEmbedding.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/NodeSet.h>
#include <ogdf/basic/graph_generators/deterministic.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/layered/FastHierarchyLayout.h>
#include <ogdf/layered/SugiyamaLayout.h>

#include <list>
#include <vector>

using namespace ogdf;
using namespace std;

#ifdef OGDF_DEBUG
#	define safe_while(cond, lim)                         \
		for (int safe_while_cnt = 0; (cond) &&            \
				[&]() {                                   \
					OGDF_ASSERT(safe_while_cnt <= (lim)); \
					return true;                          \
				}();                                      \
				safe_while_cnt++)
#else
#	define safe_while(cond, lim) while (cond)
#endif

ogdf::Logger logger;

void drawSVG(const Graph& G, const NodeArray<int>& lvl, const NodeArray<int>& pos) {
	GraphAttributes GA {G, GraphAttributes::all};
	for (node n : G.nodes) {
		GA.label(n) = to_string(n->index());
	}

	Hierarchy H {G, lvl};
	HierarchyLevels HL {H};
	GraphCopy& GC = (GraphCopy&)H;
	NodeArray<int> HLpos {GC, -1};
	for (node v : GC.nodes) {
		HLpos[v] = pos[GC.original(v)];
	}
	HL.restorePos(HLpos);
	FastHierarchyLayout L;
	L.call(HL, GA);

	GA.flipVertical();
	ogdf::GraphIO::write(GA, "test.svg");
}

struct AdjCompLess {
	const NodeArray<int>& pos;

	explicit AdjCompLess(const NodeArray<int>& pos) : pos(pos) { }

	bool operator()(adjEntry a, adjEntry b) {
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
};

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
			G.sort(n, adjs);
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

void randomProperLevelPlaneGraph(Graph& G, vector<vector<node>>& emb, NodeArray<int>& lvl,
		NodeArray<int>& pos, int N, int M, int K, vector<int>& edge_lengths) {
	// checks
	OGDF_ASSERT(N > 0);
	OGDF_ASSERT(K > 0);
	OGDF_ASSERT(N >= K);

	// init arrays
	G.clear();
	lvl.init(G, -1);
	pos.init(G, -1);
	emb.clear();
	emb.resize(K);

	// aux array for selecting edge endpoints
	NodeArray<bool> is_virt(G, false);

	// track first and last position (inclusive) on next upward layer each node sees
	NodeArray<pair<int, int>> vis(G, {-1, -1});

	// function for updating visibilities after inserting edge emb[l][up] -> emb[l+1][vp]
	auto update_vis = [&vis, &emb](int l, int up, int vp) {
		// update visible range of nodes on level l to the left
		for (int i = up - 1; i >= 0; i--) {
			if (vis[emb[l][i]].second <= vp) {
				break;
			} else {
				// logger.lout() << "\tLevel " << l << " Pos " << i << " Visibility "
				// 			  << vis[emb[l][i]].second << " -> " << vp << std::endl;
				vis[emb[l][i]].second = vp;
			}
		}
		// update visible range of nodes on level l to the right
		for (int i = up + 1; i < emb[l].size(); i++) {
			if (vis[emb[l][i]].first >= vp) {
				break;
			} else {
				// logger.lout() << "\tLevel " << l << " Pos " << i << " Visibility "
				// 			  << vis[emb[l][i]].first << " -> " << vp << std::endl;
				vis[emb[l][i]].first = vp;
			}
		}
	};

	while (edge_lengths.size() < M) {
		edge_lengths.push_back(1);
	}

	// for each long edge
	vector<pair<int, int>> long_edges;
	long_edges.reserve(edge_lengths.size());
	int subdiv_cnt = 0;
	for (int el : edge_lengths) {
		OGDF_ASSERT(el > 0 && el <= K);
		subdiv_cnt += el - 1;
		// select and store its start layer
		int l = randomNumber(0, K - 1 - el);
		long_edges.emplace_back(l, el);
		// create subdivision vertices for all long edges
		for (int i = l + 1; i < l + el; i++) {
			emb[i].push_back(G.newNode());
			// note that this is not the exact subdivision vertex we will use
			// we just ensure that there are sufficiently many vertices available on level i
		}
	}

	// ensure each level has a (non-subdivision) node
	for (int i = 0; i < K; ++i) {
		emb[i].push_back(G.newNode());
	}
	// create further nodes as requested
	for (int i = K; i < N; ++i) {
		emb[randomNumber(0, K - 1)].push_back(G.newNode());
	}

	// correctly update all arrays from emb
	for (int i = 0; i < emb.size(); ++i) {
		int next_max;
		if (i + 1 < emb.size()) {
			next_max = emb[i + 1].size() - 1;
		} else {
			next_max = -1;
		}
		for (int j = 0; j < emb[i].size(); ++j) {
			node n = emb[i][j];
			lvl[n] = i;
			pos[n] = j;
			vis[n] = {0, next_max};
		}
	}

	for (int l = 0; l < emb.size(); ++l) {
		ostream& lg = logger.lout();
		lg << "Level " << l << ": ";
		for (int i = 0; i < emb[l].size(); ++i) {
			lg << emb[l][i] << " ";
		}
		lg << std::endl;
	}
	logger.lout() << std::endl;

	// sort long edges by descending end layer
	sort(long_edges.begin(), long_edges.end(), [](const pair<int, int>& a, const pair<int, int>& b) {
		return a.first + a.second > b.first + b.second;
	});

	for (const auto& p : long_edges) {
		logger.lout() << "Long edge " << p.first << " + " << p.second << " = "
					  << (p.first + p.second) << std::endl;
	}
	logger.lout() << std::endl;

	// create long edges by descending end layer
	for (const auto& p : long_edges) {
		logger.lout() << "Creating long edge " << p.first << " + " << p.second << " = "
					  << (p.first + p.second) << std::endl;
		int sl = p.first;
		int el = p.first + p.second;
		int sp = randomNumber(0, emb[sl].size() - 1);
		safe_while(is_virt[emb[sl][sp]], emb[sl].size()) { //
			sp = (sp + 1) % (emb[sl].size());
		}
		logger.lout() << "Start pos is " << sp << std::endl;
		Logger::Indent _(logger);

		node u = emb[sl][sp];
		int up = sp;
		for (int l = sl; l < el - 1; l++) {
			// logger.lout() << "Level " << l << " Pos " << up << " Vert " << u << " Visibility "
			// 			  << vis[u].first << " - " << vis[u].second << std::endl;
			int vp = randomNumber(vis[u].first, vis[u].second);
			int cnt = vis[u].second - vis[u].first;
			while (is_virt[emb[l + 1][vp]]) {
				if (cnt < 1) {
					break;
				}
				cnt--;
				vp++;
				if (vp > vis[u].second) {
					vp = vis[u].first;
				}
			}
			if (is_virt[emb[l + 1][vp]] && cnt < 1) {
				// bad news: we see no vertex that is non-virtual; thus, we now cheat:
				// use an already-used virtual vertex, delete a superfluous isolated vertex from
				// the current level, and split the doubly-used vertex later on
				for (int i = 0; i <= emb[l + 1].size(); ++i) {
					OGDF_ASSERT(i < emb[l + 1].size());
					if (!is_virt[emb[l + 1][i]] && emb[l + 1][i]->degree() == 0) {
						is_virt[emb[l + 1][i]] = true;
						break;
					}
				}
			}
			node v = emb[l + 1][vp];
			logger.lout() << "Edge " << u << " - " << v << " Pos " << up << " - " << vp << std::endl;
			G.newEdge(u, v);

			is_virt[v] = true;
			update_vis(l, up, vp);

			u = v;
			up = vp;
		}

		//last segment
		{
			int vp = randomNumber(vis[u].first, vis[u].second);
			int cnt = vis[u].second - vis[u].first;
			while (is_virt[emb[el][vp]]) {
				if (cnt < 1) {
					break;
				}
				cnt--;
				vp++;
				if (vp > vis[u].second) {
					vp = vis[u].first;
				}
			}
			if (is_virt[emb[el][vp]] && cnt < 1) {
				// bad news: we see no vertex that is non-virtual; thus, we now cheat:
				// use an already-used virtual vertex, delete a superfluous isolated vertex from
				// the current level, and split the doubly-used vertex later on
				// OGDF_ASSERT(vis[u].first + 1 == vis[u].second);
				for (int i = 0; i <= emb[el].size(); ++i) {
					OGDF_ASSERT(i < emb[el].size());
					if (!is_virt[emb[el][i]] && emb[el][i]->degree() == 0) {
						is_virt[emb[el][i]] = true;
						break;
					}
				}
			}
			node v = emb[el][vp];
			logger.lout() << "Edge " << u << " - " << v << " Pos " << up << " - " << vp << std::endl;
			G.newEdge(u, v);
			update_vis(el - 1, up, vp);
		}

		drawSVG(G, lvl, pos);
	}

	// // create remaining requested edges as short edges
	// for (int i = long_edges.size(); i < M; i++) {
	// 	node u = G.chooseNode([&](node n) { return !is_virt(n) && lvl[n] + 1 < emb.size(); });
	// 	int l = lvl[u];
	// 	int vp = randomNumber(vis[u].first, vis[u].second);
	// 	safe_while(is_virt[emb[l + 1][vp]], vis[u].second - vis[u].first) {
	// 		vp++;
	// 		if (vp > vis[u].second) {
	// 			vp = vis[u].first;
	// 		}
	// 	}
	// 	node v = emb[l + 1][vp];
	// 	logger.lout() << "Level " << l << " Edge " << u << " - " << v << " Pos " << pos[u] << " - " << vp << std::endl;
	// 	G.newEdge(u, v);
	// 	update_vis(l, pos[u], vp);
	// }

	// sort adjacency lists
	embedPLE(G, emb, lvl, pos);
}

int main(int argc, char* argv[]) {
	Graph G;
	NodeArray<int> lvl(G, -1);
	vector<vector<node>> emb;
	NodeArray<int> pos(G, -1);
	vector<int> vec {4, 3, 3, 2, 2, 2};
	randomProperLevelPlaneGraph(G, emb, lvl, pos, 10, 20, 5, vec);

	// draw svg
	drawSVG(G, lvl, pos);

	// check
	checkPLE(G, emb, lvl, pos);

	return 0;
}
