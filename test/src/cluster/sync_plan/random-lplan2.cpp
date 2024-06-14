#include <ogdf/basic/CombinatorialEmbedding.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/graph_generators/deterministic.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/layered/FastHierarchyLayout.h>
#include <ogdf/layered/SugiyamaLayout.h>

#include <list>
#include <vector>

using namespace ogdf;
using namespace std;


ogdf::Logger logger;

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

void checkPLE(const Graph& G, const NodeArray<int>& levels, const vector<vector<node>>& emb,
		const NodeArray<int>& pos) {
	int c = 0;
	for (int i = 0; i < emb.size(); ++i) {
		for (int j = 0; j < emb[i].size(); ++j) {
			OGDF_ASSERT(levels[emb[i][j]] == i);
			OGDF_ASSERT(pos[emb[i][j]] == j);
			c++;
		}
	}
	OGDF_ASSERT(c == G.numberOfNodes());

	for (node n : G.nodes) {
		if (n->degree() < 2) {
			continue;
		}
		bool up = false;
		adjEntry pred = n->lastAdj(), curr = n->firstAdj();
		while (curr != nullptr) {
			if (curr->isSource()) {
				up = true;
				OGDF_ASSERT(levels[curr->twinNode()] == levels[n] + 1);
			} else {
				OGDF_ASSERT(levels[curr->twinNode()] == levels[n] - 1);
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
	// TODO check planarity
}

struct VisInfo {
	int cnt = -1;
	vector<node> vis;
	adjEntry right;

	vector<node> back_sink;
	vector<adjEntry> back_adj;

	adjEntry adjForVis(node n) {
		adjEntry adj;
		if (n == vis.front()) {
			if (right == nullptr) {
				return nullptr; //TODO ???
			}
			adj = right->twin()->cyclicSucc();
			if (!adj->isSource()) {
				adj = nullptr;
			}
		} else {
			OGDF_ASSERT(n->indeg() == 0);
			adj = n->firstAdj();
		}
		return adj;
	}
};

struct LevelPlanVis {
	Graph& G;
	vector<vector<node>>& emb;
	NodeArray<int>& lvl;
	NodeArray<int>& pos;

	LevelPlanVis(Graph& g, vector<vector<node>>& emb, NodeArray<int>& lvl, NodeArray<int>& pos)
		: G(g), emb(emb), lvl(lvl), pos(pos), sink_vis(G), adj_vis(G) { }

	NodeArray<VisInfo> sink_vis;
	AdjEntryArray<VisInfo> adj_vis;

	VisInfo* getVis(node n, adjEntry right) {
		// Check Levels
		int l = lvl[n];
		int ll = l;
		while (ll < emb.size() && (l == ll || emb[ll].empty())) {
			ll++;
		}
		logger.lout() << "Node " << n << " Level " << l << " Pos " << pos[n] << " Next " << ll
					  << " Right " << right << std::endl;
		if (ll == emb.size()) {
			return nullptr;
		}

		// Find next upwards edge to the right if we are a sink
		VisInfo* info;
		if (right == nullptr) {
			OGDF_ASSERT(n->outdeg() == 0);
			info = &sink_vis[n];
			for (int i = pos[n] + 1; i < emb[l].size(); i++) {
				if (emb[l][i]->outdeg() == 0) {
					continue;
				}
				right = emb[l][i]->firstAdj();
				if (emb[l][i]->indeg() > 0) {
					for (int j = emb[l][i]->outdeg(); j > 0; j--) {
						right = right->cyclicPred();
					}
				}
				break;
			}
		} else {
			OGDF_ASSERT(right->theNode() == n);
			info = &adj_vis[right];
		}
		if (info->cnt >= 0) {
			OGDF_ASSERT(right == info->right);
			return info;
		}
		info->right = right;

		// rightmost node we can see
		node rightN = nullptr;
		if (right != nullptr) {
			OGDF_ASSERT(right->isSource());
			rightN = right->twinNode();
		} else {
			rightN = emb[ll].back();
			//OGDF_ASSERT(rightN->indeg() == 0);
		}
		OGDF_ASSERT(lvl[rightN] == ll);
		ostream& lg = logger.lout();
		lg << "Right " << right << " RightN " << rightN << " Vis " << rightN << " ";

		// Find all predecessors of rightN we can also see
		auto& vis = info->vis;
		vis.emplace_back(rightN);
		// no edge right of us || seeing rightN from the left (with || without) edges further up
		if (right == nullptr || right->twin()->cyclicSucc()->isSource()
				|| right->twin()->cyclicSucc() == rightN->firstAdj()) {
			while ((vis.back()->indeg() == 0 || vis.size() == 1) && pos[vis.back()] > 0) {
				int p = pos[vis.back()] - 1;
				node n = emb[ll][p];
				lg << n << " ";
				vis.emplace_back(n);
				OGDF_ASSERT(lvl[n] == ll);
				OGDF_ASSERT(pos[n] == p);
			}
		}
		lg << std::endl;
		info->cnt = vis.size();

		// recurse
		Logger::Indent _ {logger};
		for (int i = vis.size() - 2; i >= 0; --i) {
			VisInfo* vi = getVis(vis[i], info->adjForVis(vis[i]));
			if (!vi) {
				continue;
			}
			info->cnt += vi->cnt;
			if (n->outdeg() == 0) {
				vi->back_sink.push_back(n);
			} else {
				vi->back_adj.push_back(right);
			}
		}
		return info;
	}

	void getVisAt(VisInfo& info, int idx, vector<pair<VisInfo&, node>>& path) {
		OGDF_ASSERT(info.cnt > 0);
		OGDF_ASSERT(idx < info.cnt);
		for (node n : info.vis) {
			if (idx == 0) {
				logger.lout() << "idx " << idx << " node " << n << std::endl;
				path.emplace_back(info, n);
				return;
			}
			OGDF_ASSERT(n != info.vis.back());
			adjEntry a = info.adjForVis(n);
			VisInfo* nInfo = getVis(n, a);
			logger.lout() << "idx " << idx << " node " << n << " adj " << a << " cnt " << nInfo->cnt
						  << std::endl;
			idx--;
			if (idx < nInfo->cnt) {
				Logger::Indent _ {logger};
				path.emplace_back(info, n);
				return getVisAt(*nInfo, idx, path);
			}
			idx -= nInfo->cnt;
		}
		OGDF_ASSERT(false);
	}

	edge insertEdge(node predN, adjEntry adj, int idx) {
		checkPLE(G,lvl,emb,pos);

		VisInfo& info = *getVis(predN, adj);
		vector<pair<VisInfo&, node>> path;
		getVisAt(info, idx, path);

		adjEntry predAdj = nullptr;
		if (info.right) {
			predAdj = info.right->cyclicPred();
			if (predAdj->theNode() != predN) {
				OGDF_ASSERT(predN->outdeg() == 0);
				if (predN->indeg() > 0) {
					predAdj = predN->lastAdj();
				} else {
					predAdj = nullptr;
				}
			}
		}

		for (auto& pair : path) {
			if (pair.second == path.back().second) {
				break;
			}
			node pp = G.newNode();

			if (predAdj) {
				G.newEdge(predAdj, pp);
			} else {
				G.newEdge(predN, pp);
			}
			predN = pp;
			predAdj = pp->firstAdj();

			node p = pair.second;
			lvl[pp] = lvl[p];
			pos[pp] = pos[p];
			auto it = emb[lvl[p]].begin();
			advance(it, pos[p]);
			emb[lvl[p]].insert(it, pp);
			for (int i = pos[p] + 1; i < emb[lvl[p]].size(); i++) {
				pos[emb[lvl[p]][i]]++;
			}
		}

		OGDF_ASSERT(path.size() > 1);
		VisInfo& lInfo = path.at(path.size() - 1).first;
		node last = path.back().second;
		edge e;
		if (lInfo.right && last == lInfo.right->twinNode()) {
			adjEntry lastAdj = lInfo.right->twin();

			if (predAdj) {
				e = G.newEdge(predAdj, lastAdj);
			} else {
				e = G.newEdge(predN, lastAdj);
			}
		} else {
			if (predAdj) {
				e = G.newEdge(predAdj, last);
			} else {
				e = G.newEdge(predN, last);
			}
		}

		embedPLE(G, emb, lvl, pos);
		//checkPLE(G, lvl, emb, pos);

		return e;
	}
};

void randomProperLevelPlaneGraph(Graph& G, NodeArray<int>& levels, vector<vector<node>>& emb,
		NodeArray<int>& pos, int N, int M, int K) {
	G.clear();
	levels.init(G, -1);
	pos.init(G, -1);
	emb.clear();
	emb.resize(K);
	vector<node> nodes;
	nodes.reserve(N);

	for (int i = 0; i < N; ++i) {
		node n = G.newNode();
		int l = randomNumber(0, K - 1);
		levels[n] = l;
		pos[n] = emb[l].size();
		emb[l].push_back(n);
		nodes.push_back(n);
	}

	while (G.numberOfEdges() < M) {
		node u = nodes[randomNumber(0, N - 1)];
		adjEntry a = nullptr;
		if (u->outdeg() > 0) {
			a = u->firstAdj();
			for (int i = randomNumber(0, u->outdeg()); i > 0; i--) {
				a = a->cyclicPred();
			}
		}
		LevelPlanVis LPV(G, emb, levels, pos);
		VisInfo* vis = LPV.getVis(u, a);
		if (!vis || vis->cnt<3) {
			continue;
		}
		LPV.insertEdge(u, a, randomNumber(1, vis->cnt - 2));
	}
}

int main(int argc, char* argv[]) {
	Graph G;
	NodeArray<int> lvl(G, -1);
	vector<vector<node>> emb;
	NodeArray<int> pos(G, -1);
	randomProperLevelPlaneGraph(G, lvl, emb, pos, 10, 20, 5);

	////embedPLE(G,emb,lvl,pos);
	//checkPLE(G, lvl, emb, pos);


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

	return 0;
}

int main2(int argc, char* argv[]) {
	Graph G;
	NodeArray<int> levels(G, -1);
	vector<vector<node>> emb;
	NodeArray<int> pos(G, -1);

	Array<node> N;
	customGraph(G, 25,
			{{1, 3}, {2, 4}, {2, 5}, {2, 6}, {3, 6}, {3, 9}, {4, 10}, {4, 11}, {5, 12}, {7, 14},
					{8, 15}, {9, 15}, {11, 16}, {13, 17}, {13, 18}, {13, 19}, {14, 20}, {15, 20},
					{16, 21}, {17, 21}, {19, 22}, {20, 22}, {24, 17}},
			N);

	emb = {
			{N[0], N[1]},
			{N[2], N[3]},
			{N[4], N[5], N[6], N[7], N[8], N[9]},
			{N[10], N[11], N[12], N[23], N[24], N[13], N[14], N[15]},
			{N[16], N[17], N[18], N[19], N[20]},
			{N[21], N[22]},
	};

	embedPLE(G, emb, levels, pos);
	checkPLE(G, levels, emb, pos);
	LevelPlanVis LPV(G, emb, levels, pos);
	VisInfo* info = LPV.getVis(N[3], N[3]->lastAdj());
	logger.lout() << "Cnt " << info->cnt << std::endl;
	for (int i = 0; i < info->cnt; ++i) {
		logger.lout() << i << ":" << std::endl;
		Logger::Indent _ {logger};
		vector<pair<VisInfo&, node>> path;
		LPV.getVisAt(*info, i, path);
		node n = path.back().second;
		ostream& l = logger.lout();
		l << i << ": " << n << " P: ";
		for (auto& n : path) {
			l << n.second << " ";
		}
		l << std::endl;
	}

	LPV.insertEdge(N[3], N[3]->lastAdj(), 8);

	// SugiyamaLayout SL;
	// SL.setRanking(new FixedRankingModule(levels));


	// for (node n : G.nodes) {
	// 	if (n->outdeg() == 0) {
	// 		VisInfo* vi = LPV.getVis(n, nullptr);
	// 		logger.lout() << "Sink " << n << " Cnt " << (vi ? vi->cnt : -1) << std::endl;
	// 	} else {
	// 		logger.lout() << "Node " << n << std::endl;
	// 		Logger::Indent _ {logger};
	// 		for (adjEntry adj : n->adjEntries) {
	// 			if (!adj->isSource()) {
	// 				continue;
	// 			}
	// 			VisInfo* vi = LPV.getVis(n, adj);
	// 			logger.lout() << "Adj " << adj << " Cnt " << (vi ? vi->cnt : -1) << std::endl;
	// 		}
	// 	}
	// }
}
