//*********************************************************
//  Regression test for planarity tests and embeddings
//
//  Tested classes:
//    - BoothLueker
//    - BoyerMyrvold
//    - NonPlanarCore
//
//  Author: Carsten Gutwenger, Tilo Wiedera, Mirko Wagner
//*********************************************************

#include <random>

#include <bandit/bandit.h>
#include <resources.h>

#include <ogdf/basic/graph_generators.h>
#include <ogdf/planarity/BoothLueker.h>
#include <ogdf/planarity/BoyerMyrvold.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/planarity/NonPlanarCore.h>
#include <ogdf/planarity/PlanarizationLayout.h>
#include <ogdf/planarity/SubgraphPlanarizer.h>
#include <ogdf/planarity/FixedEmbeddingInserter.h>
#include <ogdf/graphalg/MaxFlowSTPlanarItaiShiloach.h>

using namespace ogdf;
using namespace bandit;
using std::minstd_rand;
using std::uniform_int_distribution;
using ReturnType = CrossingMinimizationModule::ReturnType;

const int MIN_N = 500;
const int MAX_N = 2500;
const int STEP_N = 50;
const int STEPS = (MAX_N - MIN_N) / STEP_N;


void removeEdges(Graph &G, const List<edge> &delEdges){
	ListConstIterator<edge> it;
	for(it = delEdges.begin(); it.valid(); ++it){
		G.delEdge(*it);
	}
}

void randomizeAdjLists(Graph G, minstd_rand &rng){
	for(node v : G.nodes){
		List<adjEntry> L;
		v->allAdjEntries(L);
		L.permute(rng);
		G.sort(v, L);
	}
}

void addRandomEdges(Graph &G, int m, minstd_rand &rng){
	const int n = G.numberOfNodes();

	Array<node> nodes(n);
	int i = 0;
	for(node v : G.nodes){
		nodes[i++] = v;
	}

	uniform_int_distribution<> dist_0(0, n - 1), dist_1(1, n - 1);
	for(i = 0; i < m; ++i){
		int i = dist_0(rng);
		int j = dist_1(rng);
		G.newEdge(nodes[i], nodes[(i + j) % n]);
	}
}

void addRandomLoops(Graph &G, int m, minstd_rand &rng){
	const int n = G.numberOfNodes();

	Array<node> nodes(n);
	int i = 0;
	for(node v : G.nodes){
		nodes[i++] = v;
	}

	uniform_int_distribution<> dist(0, n - 1);
	for(i = 0; i < m; ++i){
		node v = nodes[dist(rng)];
		G.newEdge(v, v);
	}
}

void addRandomMultiEdges(Graph &G, int add, minstd_rand &rng){
	const int m = G.numberOfEdges();

	Array<edge> edges(m);
	int i = 0;
	for(edge e : G.edges){
		edges[i++] = e;
	}

	uniform_int_distribution<> dist(0, m - 1);
	for(i = 0; i < add; ++i){
		edge e = edges[dist(rng)];
		G.newEdge(e->source(), e->target());
	}
}


void describeModule(const std::string &name, PlanarityModule &pm){
describe(name, [&](){
	minstd_rand rng(42);
	srand(4711);

	Graph G;

	it("recognizes planar biconnected graphs", [&](){
		int64_t time = 0, T;

		for(int n = MIN_N; n <= MAX_N; n += STEP_N){
			int mValues[] = {3 * n / 2, 2 * n, 5 * n / 2};
			for(auto m : mValues){
				planarBiconnectedGraph(G, n, m);
				addRandomLoops(G, 10, rng);
				randomizeAdjLists(G, rng);

				System::usedRealTime(T);
				AssertThat(pm.isPlanar(G), Equals(true));
				time += System::usedRealTime(T);
			}
		}
		cout << endl << "      average time was " << time / STEPS / 3 << "ms" << endl;
	});

	it("recognizes planar connected graphs", [&](){
		int64_t time = 0, T;

		for(int n = MIN_N; n <= MAX_N; n += STEP_N){
			int nb = n / 10;
			int mValues[] = {3 * nb / 2, 2 * nb, 5 * nb / 2};
			for(auto m : mValues){
				planarCNBGraph(G, nb, m, 10);
				addRandomLoops(G, 10, rng);
				randomizeAdjLists(G, rng);

				System::usedRealTime(T);
				AssertThat(pm.isPlanar(G), Equals(true));
				time += System::usedRealTime(T);
			}
		}
		cout << endl << "      average time was " << time / STEPS / 3 << "ms" << endl;
	});

	it("creates a planar embedding on biconnected graphs", [&](){
		int64_t time = 0, T;

		for(int n = MIN_N; n <= MAX_N; n += STEP_N){
			int nb = n / 10;
			int mValues[] = {3 * nb / 2 + 1, 2 * nb + 1, 5 * nb / 2 + 1};
			for(auto m : mValues){
				planarBiconnectedGraph(G, n, m);
				addRandomLoops(G, 10, rng);
				addRandomMultiEdges(G, 50, rng);
				randomizeAdjLists(G, rng);

				System::usedRealTime(T);
				pm.planarEmbed(G);
				time += System::usedRealTime(T);
				AssertThat(G.representsCombEmbedding(), Equals(true));
			}
		}
		cout << endl << "      average time was " << time / STEPS / 3 << "ms" << endl;
	});

	it("works on connected graphs", [&](){
		int64_t time = 0, T;

		for(int n = MIN_N; n <= MAX_N; n += STEP_N){
			int nb = n / 10;
			int mValues[] = {3 * nb / 2 + 1, 2 * nb + 1, 5 * nb / 2 + 1};
			for(auto m : mValues){
				int addEdges = randomNumber(0, 4);
				planarCNBGraph(G, nb, m, 10);
				addRandomEdges(G, addEdges, rng);
				addRandomMultiEdges(G, 50, rng);
				randomizeAdjLists(G, rng);

				System::usedRealTime(T);
				time += System::usedRealTime(T);
				AssertThat(pm.isPlanar(G), Equals(isPlanar(G)));
			}
		}
		cout << endl << "      average time was " << time / STEPS / 3 << "ms" << endl;
	});

	it("detects non-planarity in small complete graphs", [&](){
		for(int i = 5; i < 50; i++){
			completeGraph(G, i);
			AssertThat(pm.isPlanar(G), Equals(i < 5));
		}
	});
});
}

void testNonPlanarCore()
{
	for_each_graph_it("returns a simple core", {"north/g.41.26.gml", "north/g.73.8.gml"},
	                  [&](Graph &graph, const string &filename){
		                  makeBiconnected(graph);
		                  NonPlanarCore npc(graph);
		                  const Graph &core = npc.core();
		                  AssertThat(isSimpleUndirected(core), IsTrue());

		                  for(edge e : core.edges){
			                  AssertThat(npc.cost(e), IsGreaterThan(0));
			                  if(!npc.isVirtual(e)){
				                  AssertThat(npc.realEdge(e), !IsNull());
			                  }
		                  }
	                  });

	it("works on a minimal previously failing instance (2 x K5)", [](){
		Graph graph;
		EdgeArray<int> weight(graph);

		node s = graph.newNode();
		node t = graph.newNode();
		graph.newEdge(t, s);

		node v = graph.newNode();
		graph.newEdge(s, v);
		graph.newEdge(v, t);

		for(int k = 0; k < 2; k++){
			List<node> nodes;
			nodes.pushBack(s);
			nodes.pushBack(t);

			for(int i = 0; i < 3; i++){
				nodes.pushBack(graph.newNode());
			}

			for(node v : nodes){
				for(node w : nodes){
					if(v->index() < w->index() && (v != s || w != t)){
						graph.newEdge(v, w);
					}
				}
			}
		}

		NonPlanarCore npc(graph);
		const Graph &core = npc.core();

		for(edge e : core.edges) {
			if(npc.isVirtual(e)) {
				for(auto eCut : npc.mincut(e)) {
					if(eCut.e->source() == npc.original(e->source()) || eCut.e->target() == npc.original(e->target())) {
						AssertThat(eCut.dir, IsTrue());
					} else {
						AssertThat(eCut.dir, IsFalse());
					}
				}
			}
		}
		AssertThat(isLoopFree(core), IsTrue());
		AssertThat(isSimpleUndirected(core), IsTrue());
		AssertThat(core.numberOfNodes(), Equals(graph.numberOfNodes() - 1));
		AssertThat(core.numberOfEdges(), Equals(graph.numberOfEdges() - 2));
	});

	it("recognizes weight with the default MaxFlowModule Dijkstra", [&](){
		Graph graph;
		completeGraph(graph, 5);
		EdgeArray<int> weight(graph, 1);
		edge e = graph.chooseEdge();
		edge f = graph.newEdge(e->target(), e->source());
		weight[graph.split(e)] = 32;
		weight[graph.split(f)] = 64;
		weight[graph.newEdge(e->target(), f->target())] = 4;
		weight[e] = 8;
		weight[f] = 16;
		NonPlanarCore npc(graph, weight);
		const Graph &core = npc.core();
		for(edge e : core.edges) {
			if(npc.isVirtual(e)){
				AssertThat(npc.cost(e), Equals(28));
			}
		}
	});

	it("recognizes weight with ItaiShiloach as the MaxFlowModule", [&](){
		Graph graph;
		completeGraph(graph, 5);
		EdgeArray<int> weight(graph, 1);
		edge e = graph.chooseEdge();
		edge f = graph.newEdge(e->target(), e->source());
		weight[graph.split(e)] = 32;
		weight[graph.split(f)] = 64;
		weight[graph.newEdge(e->target(), f->target())] = 4;
		weight[e] = 8;
		weight[f] = 16;
		MaxFlowSTPlanarItaiShiloach<int> maxflow;
		NonPlanarCore npc(graph, weight, maxflow);
		const Graph &core = npc.core();
		for(edge e : core.edges) {
			if(npc.isVirtual(e)){
				AssertThat(npc.cost(e), Equals(28));
			}
		}
	});

	for_each_graph_it("retransforms", {"north/g.41.26.gml", "north/g.73.8.gml"}, [&](Graph &graph, const string &file) {
		makeBiconnected(graph);
		NonPlanarCore C(graph);
		const Graph &core = C.core();
		AssertThat(isPlanar(core), IsFalse());
		AssertThat(core.numberOfNodes(), !Equals(0));
		SubgraphPlanarizer SP;
		PlanRep planarCore(core);

		GraphCopy endGraph(graph);
		int crossingNumber = 0;
		ReturnType ret = SP.call(planarCore, 0, crossingNumber, &C.cost());
		AssertThat(ret == ReturnType::retTimeoutFeasible || ret == ReturnType::retFeasible ||
				   ret == ReturnType::retOptimal, IsTrue());
		AssertThat(planarEmbed(planarCore), IsTrue());
		planarCore.removePseudoCrossings();

		C.retransform(planarCore, endGraph);

		AssertThat(isPlanar(endGraph), IsTrue());
		AssertThat(endGraph.genus(), Equals(0));

		// now the embedding of the endGraph is tested to assert that the embedding of planarCore was
		// used to embed the endGraph
		for(node v : planarCore.nodes){
			if(planarCore.isDummy(v)){
				continue;
			}
			node endNode = endGraph.copy(C.original(planarCore.original(v)));
			List<adjEntry> adjEntries;
			endNode->allAdjEntries(adjEntries);
			int stComponentCounter(0);
			Array<int> componentList;
			componentList.grow(adjEntries.size(), -1);
			for(adjEntry pcAdj : v->adjEntries){
				edge coreEdge = planarCore.original(pcAdj->theEdge());
				node stNode = (pcAdj == pcAdj->theEdge()->adjSource() ? C.sNode(coreEdge) : C.tNode(coreEdge));
				EdgeArray<edge> &mapE = *C.mapE(coreEdge);
				for(adjEntry stAdj : stNode->adjEntries){
					List<edge> chain = endGraph.chain(mapE[stAdj->theEdge()]);
					adjEntry endAdj = nullptr;
					for(edge e : chain){
						if(e->source() == endNode){
							endAdj = e->adjSource();
						}
						if(e->target() == endNode){
							endAdj = e->adjTarget();
						}

					}
					auto searchIt = adjEntries.search(endAdj);
					AssertThat(searchIt.valid(), IsTrue());
					int position = adjEntries.pos(searchIt);
					componentList[position] = stComponentCounter;
				}
				stComponentCounter++;
			}
			int before(*componentList.rbegin());
			for(int i : componentList){
				if(i != before){
					AssertThat((before + 1) % stComponentCounter, Equals(i))
					before = i;
				}
			}
		}
	});

	it("contracts chains", [&](){
		Graph graph;
		GraphAttributes GA(graph);
		GA.initAttributes(GraphAttributes::nodeType | GraphAttributes::edgeType | GraphAttributes::nodeLabel |
		                  GraphAttributes::nodeStyle | GraphAttributes::edgeLabel | GraphAttributes::edgeStyle |
		                  GraphAttributes::edgeArrow);
		for(int i = 0; i < 13; i++){
			node curr = graph.newNode();
			GA.label(curr) = to_string(curr->index());
			GA.fillColor(curr) = Color::Turquoise;
		}

		List<node> v;
		graph.allNodes(v);

		graph.newEdge(*v.get(0), *v.get(1));
		graph.newEdge(*v.get(1), *v.get(2));
		graph.newEdge(*v.get(2), *v.get(4));
		graph.newEdge(*v.get(1), *v.get(3));
		graph.newEdge(*v.get(4), *v.get(3));
		graph.newEdge(*v.get(3), *v.get(5));
		graph.newEdge(*v.get(5), *v.get(6));
		graph.newEdge(*v.get(5), *v.get(2));
		graph.newEdge(*v.get(4), *v.get(6));
		edge e67 = graph.newEdge(*v.get(6), *v.get(7));
		edge e78 = graph.newEdge(*v.get(7), *v.get(8));
		graph.newEdge(*v.get(0), *v.get(11));
		graph.newEdge(*v.get(0), *v.get(10));
		graph.newEdge(*v.get(11), *v.get(12));
		graph.newEdge(*v.get(10), *v.get(12));
		graph.newEdge(*v.get(10), *v.get(9));
		graph.newEdge(*v.get(9), *v.get(8));
		graph.newEdge(*v.get(5), *v.get(4));
		graph.newEdge(*v.get(12), *v.get(8));
		graph.newEdge(*v.get(11), *v.get(9));

		EdgeArray<int> weight(graph, 1);
		weight[e67] = 2;
		weight[e78] = 3;
		NonPlanarCore C(graph, weight);
		const Graph &core = C.core();
		node v6, v8;
		for(node w : core.nodes){
			if(C.original(w) == *v.get(6)){
				v6 = w;
			}
			if(C.original(w) == *v.get(8)){
				v8 = w;
			}
		}
		edge virt = nullptr;
		for(edge e : core.edges){
			if((e->source() == v6 && e->target() == v8)
			   || (e->source() == v8 && e->target() == v6)){
				virt = e;
			}
		}
		AssertThat(virt, !Equals((void *) nullptr));
		AssertThat(C.isVirtual(virt), IsTrue());
		AssertThat(C.cost(virt), Equals(2));
	});

	it("eliminates multiedges", [](){
		Graph graph;
		completeGraph(graph, 5);
		edge e = graph.chooseEdge();
		graph.newEdge(e->source(), e->target());
		e = graph.chooseEdge();
		graph.newEdge(e->target(), e->source());
		NonPlanarCore npc(graph);
		const Graph &core = npc.core();
		AssertThat(isSimpleUndirected(core), IsTrue());
		AssertThat(core.numberOfNodes(), Equals(graph.numberOfNodes()));
		AssertThat(core.numberOfEdges(), Equals(10));
	});

	it("returns a list of original edges of a core edge", [](){
		Graph graph;
		completeGraph(graph, 5);
		edge e = graph.chooseEdge();
		edge f = graph.split(e);
		NonPlanarCore npc(graph);
		for(edge eCore : npc.core().edges) {
			List<edge> list = npc.original(eCore);
			if(npc.isVirtual(eCore)){
				AssertThat(list.size(), Equals(2));
				if(list.front() == e) {
					AssertThat(list.back(), Equals(f));
				} else {
					AssertThat(list.front(), Equals(f));
					AssertThat(list.back(), Equals(e));
				}
			} else {
				AssertThat(list.size(), Equals(1));
				AssertThat(list.front(), Equals(npc.realEdge(eCore)));
			}
		}
	});
}

go_bandit([](){
	describe("Planarity tests", [](){
		BoothLueker bl;
		describeModule("Booth-Lueker", bl);
		BoyerMyrvold bm;
		describeModule("Boyer-Myrvold", bm);

		it("transforms based on the right graph, when it's a GraphCopySimple", [&](){
			Graph G;
			randomRegularGraph(G, 10, 6);
			GraphCopySimple gcs(G);
			BoyerMyrvold bm;
			SList<KuratowskiWrapper> kur_subs;
			SList<KuratowskiSubdivision> lksGcs;
			SList<KuratowskiSubdivision> lksG;

			bm.planarEmbed(gcs, kur_subs, BoyerMyrvoldPlanar::doFindUnlimited);
			bm.transform(kur_subs,lksGcs,gcs);
			bm.transform(kur_subs,lksG,G);
		});
	});

	describe("NonPlanarCore", [](){
		testNonPlanarCore();
	});
});
