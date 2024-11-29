/** \file
 * \brief Implementation of some randomized clustering generators
 *
 * \author Carsten Gutwenger, Markus Chimani, Jöran Schierbaum, Simon D. Fink
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


#include <ogdf/basic/Array.h>
#include <ogdf/basic/CombinatorialEmbedding.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/graph_generators/clustering.h>
#include <ogdf/basic/graph_generators/randomized.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/cluster/ClusterGraph.h>
#include <ogdf/cluster/sync_plan/PMatching.h>
#include <ogdf/cluster/sync_plan/SyncPlan.h>
#include <ogdf/cluster/sync_plan/utils/Bijection.h>

#include <chrono>
#include <cstdint>
#include <functional>
#include <ostream>
#include <random>
#include <unordered_set>
#include <vector>

using std::minstd_rand;
using std::uniform_int_distribution;
using std::uniform_real_distribution;
using std::unordered_set;
using tpc = std::chrono::high_resolution_clock;
using tp = std::chrono::time_point<std::chrono::high_resolution_clock>;

namespace ogdf {

static void constructCConnectedCluster(node v, ClusterGraph& C, minstd_rand& rng);
static void constructCluster(node v, ClusterGraph& C);
static void bfs(node v, SList<node>& newCluster, NodeArray<bool>& visited, ClusterGraph& C,
		minstd_rand& rng);

void randomClustering(ClusterGraph& C, int cNum) {
	auto& G = C.constGraph();
	int n = G.numberOfNodes();

	int count = 0;
	NodeArray<int> num(G);
	Array<node> numNode(0, n - 1, nullptr);
	for (node v : G.nodes) {
		num[v] = count;
		numNode[count] = v;
		count++;
	}

	minstd_rand rng(randomSeed());
	uniform_int_distribution<> dist(0, n - 1);

	for (int i = 0; i < cNum; i++) {
		constructCluster(numNode[dist(rng)], C);
	}

#ifdef OGDF_DEBUG
	C.consistencyCheck();
#endif
}

void randomCConnectedClustering(ClusterGraph& C, int cNum) {
	auto& G = C.constGraph();
	int n = G.numberOfNodes();

	int count = 0;
	NodeArray<int> num(G);
	Array<node> numNode(0, n - 1, nullptr);
	for (node v : G.nodes) {
		num[v] = count;
		numNode[count] = v;
		count++;
	}

	minstd_rand rng(randomSeed());
	uniform_int_distribution<> dist(0, n - 1);

	for (int i = 0; i < cNum; i++) {
		constructCConnectedCluster(numNode[dist(rng)], C, rng);
	}

	// By construction, clusters might have just one child.
	// remove these clusters
	SListPure<cluster> store;
	for (cluster c : C.clusters) {
		if ((c->cCount() + c->nCount()) == 1) {
			store.pushBack(c);
		}
	}
	while (!store.empty()) {
		cluster c = store.popFrontRet();
		if (c != C.rootCluster()) {
			C.delCluster(c);
		}
	}
	if ((C.rootCluster()->cCount() == 1) && (C.rootCluster()->nCount() == 0)) {
		cluster cl = *C.rootCluster()->cBegin();
		C.delCluster(cl);
	}

#ifdef OGDF_DEBUG
	C.consistencyCheck();
#endif
}

static void constructCConnectedCluster(node v, ClusterGraph& C, minstd_rand& rng) {
	SList<node> newCluster;
	newCluster.pushBack(v);
	NodeArray<bool> visited(C, false);
	visited[v] = true;
	bfs(v, newCluster, visited, C, rng);
	if (newCluster.size() > 1) {
		cluster cl = C.newCluster(C.clusterOf(v));
		while (!newCluster.empty()) {
			node w = newCluster.popFrontRet();
			C.reassignNode(w, cl);
		}
	}
}

// Construct new (child) cluster by randomly choosing nodes in v's cluster
static void constructCluster(node v, ClusterGraph& C) {
	if (C.clusterOf(v)->nCount() < 2) {
		return;
	}

	SList<node> newCluster;
	newCluster.pushBack(v);

	minstd_rand rng(randomSeed());
	uniform_int_distribution<> dist(0, 99);

	// store the cluster nodes for random selection
	// we could just randomly select by running up the list
	for (node u : C.clusterOf(v)->nodes) {
		if (u != v && dist(rng) > 65) {
			newCluster.pushBack(u);
		}
	}

	cluster cl = C.newCluster(C.clusterOf(v));
	while (!newCluster.empty()) {
		node w = newCluster.popFrontRet();
		C.reassignNode(w, cl);
	}
}

// Insert nodes in v's cluster to new cluster with a certain probability
static void bfs(node v, SList<node>& newCluster, NodeArray<bool>& visited, ClusterGraph& C,
		minstd_rand& rng) {
	uniform_int_distribution<> dist(0, 99);

	SListPure<node> bfsL;
	for (adjEntry adj : v->adjEntries) {
		edge e = adj->theEdge();
		node w = e->opposite(v);
		int probability = dist(rng);
		if (probability < 70 && !visited[w]) {
			visited[w] = true;
			if (C.clusterOf(v) == C.clusterOf(w)) {
				newCluster.pushBack(w);
				bfsL.pushBack(w);
			}
		} else {
			visited[w] = true;
		}
	}
	while (!bfsL.empty()) {
		bfs(bfsL.popFrontRet(), newCluster, visited, C, rng);
	}
}

void createClustersHelper(ClusterGraph& C, const node curr, const node pred, const cluster predC,
		List<cluster>& internal, List<cluster>& leaves) {
	cluster currC = predC ? C.createEmptyCluster(predC) : C.rootCluster();
	if (curr->degree() == 1 && pred != nullptr) {
		leaves.pushBack(currC);
	} else {
		for (adjEntry adj : curr->adjEntries) {
			node next = adj->twinNode();
			if (next == pred) {
				continue;
			}
			createClustersHelper(C, next, curr, currC, internal, leaves);
		}
		internal.pushBack(currC);
	}
}

void randomClustering(ClusterGraph& C, const node root, int moreInLeaves) {
	auto& G = C.constGraph();
	C.init(G);

	// Build cluster structure (and store which clusters are internal and which are leaves)
	List<cluster> internal;
	List<cluster> leaves;
	createClustersHelper(C, root, nullptr, nullptr, internal, leaves);

	// Assign nodes to clusters
	List<node> nodes;
	G.allNodes<List<node>>(nodes);

	// Step 1: Ensure two node per leaf-cluster
	nodes.permute();
	for (cluster c : leaves) {
		C.reassignNode(nodes.popFrontRet(), c);
		C.reassignNode(nodes.popFrontRet(), c);
	}

	// Step 2: Distribute the other nodes
	int n = G.numberOfNodes();
	int numI = internal.size();
	int numL = leaves.size();
	double chanceForInternal = (numI * n / double(numL * moreInLeaves + numI)) / double(n - 2 * numL);
	// a leaf-cluster should have (on average) moreInLeaves-times as many vertices as in internal-cluster.
	// #verticesInInternalCluster = n / (numL*moreInLeaves + numI)
	// #nodesToDistribute = n - 2*numL
	// => chance that a node goes into an internal cluster = numI * #verticesInInternalCluster / (n-2*numL)

	minstd_rand rng(randomSeed());
	uniform_real_distribution<> dist_0_1(0.0, 1.0);

	while (!nodes.empty()) {
		cluster cl;
		if (dist_0_1(rng) < chanceForInternal) {
			cl = *internal.get(uniform_int_distribution<>(0, internal.size() - 1)(rng));
		} else {
			cl = *leaves.get(uniform_int_distribution<>(0, leaves.size() - 1)(rng));
		}
		C.reassignNode(nodes.popFrontRet(), cl);
	}
}

std::ostream& operator<<(std::ostream& os, const RandomClusterConfig& config) {
	os << "max_nodes_in_cluster: " << config.max_nodes_in_cluster
	   << " prob_no_further_node: " << config.prob_no_further_node << " ("
	   << config.expected_nodes() << ")"
	   << " prob_no_further_cluster: " << config.prob_no_further_cluster << " ("
	   << 1.0 / config.prob_no_further_cluster << ")"
	   << " max_clusters: " << config.max_clusters << " min_root_nodes: " << config.min_root_nodes
	   << " timeout: " << config.timeout;
	return os;
}

class Clusterer {
	ClusterGraph& CG;
	const RandomClusterConfig& config;

	GraphCopy copy;
	NodeArray<cluster> clusters;
	NodeArray<bool> mark;
	int marked = 0;
	tp stop;

public:
	Clusterer(ClusterGraph& _CG, const RandomClusterConfig& _config)
		: CG(_CG), config(_config), copy(CG.constGraph()), clusters(copy, nullptr), mark(copy, false) {
		OGDF_ASSERT(CG.getGraph()->representsCombEmbedding());
		OGDF_ASSERT(CG.numberOfClusters() == 1);
		CG.adjAvailable(true);
		copy.setOriginalEmbedding();
		if (!isConnected(copy)) {
			makeConnected(copy);
		}
		triangulate(copy);
		OGDF_ASSERT(copy.representsCombEmbedding());

		if (config.timeout != 0) {
			stop = tpc::now() + std::chrono::seconds(config.timeout);
		} else {
			stop = tp::max();
		}
	}

	bool timedout() const { return tpc::now() > stop; }

	bool hasFreeNodes() const {
		return !timedout() && CG.rootCluster()->nCount() > config.min_root_nodes;
	}

	bool shouldAddCluster() const {
		return hasFreeNodes()
				&& (config.max_clusters == 0 || CG.numberOfClusters() < config.max_clusters)
				&& randomDouble(0, 1) > config.prob_no_further_cluster;
	}

	bool shouldAddNode(cluster c) const {
		return hasFreeNodes()
				&& (config.max_nodes_in_cluster == 0 || c->nCount() < config.max_nodes_in_cluster)
				&& randomDouble(0, 1) > config.prob_no_further_node;
	}

	bool isTriangulated() {
		for (node n : copy.nodes) {
			for (adjEntry adj : n->adjEntries) {
				adjEntry succ2 = adj->faceCycleSucc()->faceCycleSucc();
				if (succ2 != adj && succ2->faceCycleSucc() != adj) {
					return false;
				}
			}
		}
		return true;
	}

	void makeClusters() {
		while (shouldAddCluster()) {
			node n = copy.chooseNode([](node u) -> bool { return u->degree() > 0; });
			if (!n) {
				break;
			}
			cluster c = CG.createEmptyCluster();
			mergeIntoCluster(n, c);
			// l.lout() << "Created cluster " << c->index() << " with node " << n->index() << "."
			// 		 << std::endl;
			// Logger::Indent _(l);

			OGDF_ASSERT(!mark[n]);
			for (adjEntry adj : n->adjEntries) {
				if (!mark[adj->twinNode()]) {
					mark[adj->twinNode()] = true;
					marked++;
				}
			}

			while (shouldAddNode(c) && n->degree() > 0) {
				adjEntry rand_adj = nullptr;
				adjEntry adj = n->adjEntries.head();
				int pos = randomNumber(0, n->degree() - 1);
				if (pos < n->degree() / 2) {
					for (int i = 0; i < pos; ++i) {
						adj = adj->cyclicSucc();
					}
				} else {
					for (int i = 0; i < n->degree() - pos; ++i) {
						adj = adj->cyclicPred();
					}
				}
				for (int i = 0; i < n->degree(); ++i) {
					if (!canContract(adj)) {
						adj = adj->cyclicSucc();
					} else {
						rand_adj = adj;
						break;
					}
				}
				if (!rand_adj) {
					break;
				}

				OGDF_ASSERT(canContract(rand_adj));
				contractIntoCluster(rand_adj);

				OGDF_ASSERT(copy.representsCombEmbedding());
				node c1, c2;
				OGDF_ASSERT(isTriconnected(copy, c1, c2));
				OGDF_ASSERT(isTriangulated());
			}

			OGDF_ASSERT(!mark[n]);
			for (adjEntry adj : n->adjEntries) {
				if (mark[adj->twinNode()]) {
					mark[adj->twinNode()] = false;
					marked--;
				}
			}
			OGDF_ASSERT(marked == 0);
			makeClusterAdjs(n);
		}
		OGDF_ASSERT(CG.representsCombEmbedding());
	}

	bool canContract(adjEntry adj) {
		if (config.cconnected && !copy.original(adj->theEdge())) {
			return false;
		}
		node neigh1 = nullptr, neigh2 = nullptr;
		for (adjEntry a : adj->twinNode()->adjEntries) {
			if (!mark[a->twinNode()]) {
				continue;
			}
			if (a->twinNode() == neigh1) {
				continue;
			}
			if (a->twinNode() == neigh2) {
				continue;
			}

			if (neigh1 == nullptr) {
				neigh1 = a->twinNode();
			} else if (neigh2 == nullptr) {
				neigh2 = a->twinNode();
			} else {
				return false;
			}
		}
		return true;
	}

	void contractIntoCluster(adjEntry adj) {
		edge to_contract = adj->theEdge();
		node u = adj->theNode();
		node v = adj->twinNode();
		OGDF_ASSERT(!to_contract->isSelfLoop());

		// l.lout() << "Contracting node " << v << " (°" << v->degree() << ") into the cluster of "
		// 		 << u << " (°" << u->degree() << ")" << std::endl;
		mergeIntoCluster(v, clusters[u]);

		OGDF_ASSERT(mark[v] != false);
		mark[v] = false;
		marked--;

		for (adjEntry a : v->adjEntries) {
			if (a->twinNode() != u && !mark[a->twinNode()]) {
				mark[a->twinNode()] = true;
				marked++;
			}
		}
		if (!adj->isSource()) {
			copy.reverseEdge(adj->theEdge());
		}
		node r = copy.contract(adj->theEdge());
		OGDF_ASSERT(r == u);
#ifdef OGDF_DEBUG
		CG.consistencyCheck();
#endif
	}

	void mergeIntoCluster(node n, cluster c) {
		if (clusters[n] != nullptr) {
			if (clusters[n] == c) {
				return;
			}
			OGDF_ASSERT(clusters[n] != nullptr);
			CG.moveCluster(clusters[n], c);
		} else {
			CG.reassignNode(copy.original(n), c);
		}
		clusters[n] = c;
	}

	void makeClusterAdjs(node n) {
		// std::ostream& ls = l.lout();
		// ls << "Cluster border crossing adjacencies: ";
		cluster c = clusters[n];
		c->adjEntries.clear();
		for (adjEntry adj : n->adjEntries) {
			edge Ge = copy.original(adj->theEdge());
			if (!Ge) {
				continue;
			}
			if (c->isDescendant(CG.clusterOf(Ge->source()), true)) {
				OGDF_ASSERT(!c->isDescendant(CG.clusterOf(Ge->target()), true));
				c->adjEntries.pushBack(Ge->adjSource());
				// ls << Ge->source()->index() << " --" << Ge->index() << "-> "
				//    << Ge->target()->index() << ", ";
			} else {
				OGDF_ASSERT(c->isDescendant(CG.clusterOf(Ge->target()), true));
				c->adjEntries.pushBack(Ge->adjTarget());
				// ls << Ge->target()->index() << " <-" << Ge->index() << "-- "
				//    << Ge->source()->index() << ", ";
			}
		}
		// ls << std::endl;
		CG.adjAvailable(true);
		OGDF_HEAVY_ASSERT(CG.representsCombEmbedding());
	}
};

bool randomPlanarClustering(ClusterGraph& CG, const RandomClusterConfig& config) {
	Clusterer c(CG, config);
	c.makeClusters();
	return !c.timedout();
}

void randomClusterPlanarGraph(Graph& G, ClusterGraph& CG, int clusters, int node_per_cluster,
		int edges_per_cluster) {
	OGDF_ASSERT(&CG.constGraph() == &G);
	// connected is ok as cut vertices will turn into disconnected clusters
	randomPlanarConnectedGraph(G, node_per_cluster, edges_per_cluster);

	for (int i = 0; i < clusters; ++i) {
		node n = G.chooseNode([](node x) { return x->degree() >= 4; });
		if (n == nullptr) {
			break;
		}
		cluster c = CG.newCluster(CG.clusterOf(n));

		Graph H;
		randomPlanarConnectedGraph(H, node_per_cluster, edges_per_cluster);
		node u = H.chooseNode([n](node x) { return x->degree() == n->degree(); });
		if (u == nullptr) {
			u = H.chooseNode([n](node x) { return x->degree() >= n->degree(); });
			if (u != nullptr) {
				while (u->degree() > n->degree()) {
					H.delEdge(u->adjEntries.head()->theEdge());
				}
			} else {
				CG.delCluster(c);
				continue;
			}
		}
		NodeArray<node> nodeMap(H, nullptr);
		EdgeArray<edge> edgeMap(H, nullptr);
		G.insert(H, nodeMap, edgeMap);

		for (node x : H.nodes) {
			CG.reassignNode(nodeMap[x], c);
		}

		node v = nodeMap[u];
		OGDF_ASSERT(n->degree() == v->degree());
		using namespace sync_plan;
		PipeBij bij;
		getPipeBijection(n, v, bij);
		sync_plan::join(G, n, v, bij);
	}
}

void randomSyncPlanInstance(sync_plan::SyncPlan& pq, int pipe_count, int min_deg) {
	for (int i = 0; i < pipe_count; ++i) {
		node u = pq.G->chooseNode(
				[&](node n) { return !pq.matchings.isMatchedPVertex(n) && n->degree() >= min_deg; });
		if (u == nullptr) {
			return;
		}
		node v = pq.G->chooseNode([&](node n) {
			return !pq.matchings.isMatchedPVertex(n) && n->degree() == u->degree() && n != u;
		});
		if (v == nullptr) {
			return;
		}
		pq.matchings.matchNodes(u, v);
	}
}

void addEdges(Graph* g, std::vector<edge>& added, int cnt) {
	CombinatorialEmbedding E(*g);
	for (int i = 0; i < cnt; ++i) {
		face f = E.chooseFace([](face ff) { return ff->size() > 3; });
		if (f == nullptr) {
			return;
		}
		adjEntry a = f->firstAdj();
		for (int j = 0; j < randomNumber(0, f->size() - 1); ++j) {
			a = a->faceCycleSucc();
		}
		adjEntry b = a;
		for (int j = 0; j < randomNumber(2, f->size() - 2); ++j) {
			b = b->faceCycleSucc();
		}
		OGDF_ASSERT(b != a);
		OGDF_ASSERT(b != a->faceCycleSucc());
		OGDF_ASSERT(b != a->faceCyclePred());
		if (a->theNode() == b->theNode()) {
			continue;
		}
		if (g->searchEdge(a->theNode(), b->theNode())) {
			continue;
		}
		added.push_back(E.splitFace(a, b));
	}
}

void randomSEFEInstanceBySharedGraph(Graph* sefe, EdgeArray<uint8_t, true>& edge_types, int edges1,
		int edges2) {
	OGDF_ASSERT(!sefe->empty());
	OGDF_ASSERT(isConnected(*sefe));
	OGDF_ASSERT(sefe->representsCombEmbedding());
	OGDF_ASSERT(edge_types.graphOf() == sefe);
	for (edge e : sefe->edges) {
		edge_types[e] = 3;
	}

	std::vector<edge> added1;
	addEdges(sefe, added1, edges1);
	Graph::HiddenEdgeSet h1(*sefe);
	for (edge e : added1) {
		edge_types[e] = 1;
		h1.hide(e);
	}

	std::vector<edge> added2;
	addEdges(sefe, added2, edges2);
	for (edge e : added2) {
		edge_types[e] = 2;
	}

	h1.restore();
}

void randomSEFEInstanceByUnionGraph(const Graph* sefe, EdgeArray<uint8_t, true>& edge_types,
		double frac_shared, double frac_g1) {
	OGDF_ASSERT(edge_types.graphOf() == sefe);
	for (edge e : sefe->edges) {
		double r = randomDouble(0, 1);
		if (r < frac_shared) {
			edge_types[e] = 3;
		} else if (r < frac_shared + frac_g1) {
			edge_types[e] = 1;
		} else {
			edge_types[e] = 2;
		}
	}
}
}
