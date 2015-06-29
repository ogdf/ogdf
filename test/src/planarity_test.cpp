//*********************************************************
//  Regression test for planarity tests and embeddings
//
//  Tested classes:
//    - BoothLueker
//    - BoyerMyrvold
//
//  Author: Carsten Gutwenger, Tilo Wiedera
//*********************************************************

#include <random>

#include <bandit/bandit.h>
#include <resources.h>

#include <ogdf/basic/graph_generators.h>
#include <ogdf/basic/graph_generators.h>
#include <ogdf/planarity/BoothLueker.h>
#include <ogdf/planarity/BoyerMyrvold.h>
#include <ogdf/planarity/FastPlanarSubgraph.h>
#include <ogdf/planarity/MaximalPlanarSubgraphSimple.h>
#include <ogdf/planarity/BoyerMyrvoldSubgraph.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/planarity/NonPlanarCore.h>

using namespace ogdf;
using namespace bandit;
using std::minstd_rand;
using std::uniform_int_distribution;

const int MIN_N = 500;
const int MAX_N = 2500;
const int STEP_N = 50;
const int STEPS = (MAX_N - MIN_N)/STEP_N;


void removeEdges(Graph &G, const List<edge> &delEdges)
{
	ListConstIterator<edge> it;
	for(it = delEdges.begin(); it.valid(); ++it)
		G.delEdge(*it);
}

void randomizeAdjLists(Graph G, minstd_rand &rng)
{
	for(node v : G.nodes) {
		List<adjEntry> L;
		G.adjEntries(v,L);
		L.permute(rng);
		G.sort(v,L);
	}
}

void addRandomEdges(Graph &G, int m, minstd_rand &rng)
{
	const int n = G.numberOfNodes();

	Array<node> nodes(n);
	int i = 0;
	for (node v : G.nodes)
		nodes[i++] = v;

	uniform_int_distribution<> dist_0(0,n-1), dist_1(1,n-1);
	for(i = 0; i < m; ++i) {
		int i = dist_0(rng);
		int j = dist_1(rng);
		G.newEdge(nodes[i],nodes[(i+j)%n]);
	}
}

void addRandomLoops(Graph &G, int m, minstd_rand &rng)
{
	const int n = G.numberOfNodes();

	Array<node> nodes(n);
	int i = 0;
	for(node v : G.nodes)
		nodes[i++] = v;

	uniform_int_distribution<> dist(0,n-1);
	for(i = 0; i < m; ++i) {
		node v = nodes[dist(rng)];
		G.newEdge(v,v);
	}
}

void addRandomMultiEdges(Graph &G, int add, minstd_rand &rng)
{
	const int m = G.numberOfEdges();

	Array<edge> edges(m);
	int i = 0;
	for(edge e : G.edges)
		edges[i++] = e;

	uniform_int_distribution<> dist(0,m-1);
	for(i = 0; i < add; ++i) {
		edge e = edges[dist(rng)];
		G.newEdge(e->source(), e->target());
	}
}


void describeModule(const std::string &name, PlanarityModule &pm)
{
describe(name.c_str(), [&](){
	minstd_rand rng(42);
	srand(4711);

	Graph G;
	FastPlanarSubgraph fps;
	fps.runs(16);

	it("recognizes planar biconnected graphs", [&](){
		int64_t time = 0, T;

		for(int n = MIN_N; n <= MAX_N; n += STEP_N)
		{
			int mValues[] = { 3*n/2, 2*n, 5*n/2 };
			for(int i = 0; i < 3; i++) {
				int m = mValues[i];

				planarBiconnectedGraph(G, n, m);
				addRandomLoops(G, 10, rng);
				randomizeAdjLists(G, rng);

				System::usedRealTime(T);
				AssertThat(pm.isPlanar(G), Equals(true));
				time += System::usedRealTime(T);
			}
		}
		cout << endl << "      average time was " << time/STEPS/3 << "ms" << endl;
	});

	it("recognizes planar connected graphs", [&](){
		int64_t time = 0, T;

		for(int n = MIN_N; n <= MAX_N; n += STEP_N) {
			int nb = n/10;
			int mValues[] = { 3*nb/2, 2*nb, 5*nb/2 };
			for(int i = 0; i < 3; i++) {
				int m = mValues[i];

				planarCNBGraph(G, nb, m, 10);
				addRandomLoops(G, 10, rng);
				randomizeAdjLists(G, rng);

				System::usedRealTime(T);
				AssertThat(pm.isPlanar(G), Equals(true));
				time += System::usedRealTime(T);
			}
		}
		cout << endl << "      average time was " << time/STEPS/3 << "ms" << endl;
	});

	it("creates a planar embedding on biconnected graphs", [&](){
		int64_t time = 0, T;

		for(int n = MIN_N; n <= MAX_N; n += STEP_N) {
			int nb = n/10;
			int mValues[] = { 3*nb/2+1, 2*nb+1, 5*nb/2+1 };
			for(int i = 0; i < 3; i++) {
				int m = mValues[i];

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
		cout << endl << "      average time was " << time/STEPS/3 << "ms" << endl;
	});

	it("works on connected graphs", [&](){
		int64_t time = 0, T;

		for(int n = MIN_N; n <= MAX_N; n += STEP_N) {
			int nb = n/10;
			int mValues[] = { 3*nb/2+1, 2*nb+1, 5*nb/2+1 };
			for(int i = 0; i < 3; i++) {
				int m = mValues[i];

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
		cout << endl << "      average time was " << time/STEPS/3 << "ms" << endl;
	});

	it("works on biconnected graphs with removed fast planar subgraph", [&](){
		int64_t time = 0, T;

		for(int n = MIN_N; n <= MAX_N; n += STEP_N) {
			int nb = n/10;
			int mValues[] = { 3*nb/2+1, 2*nb+1, 5*nb/2+1 };
			for(int i = 0; i < 3; i++) {
				int m = mValues[i];

				planarBiconnectedGraph(G, n, m);
				addRandomMultiEdges(G, 5, rng);
				addRandomEdges(G, 5, rng);
				randomizeAdjLists(G, rng);

				List<edge> delEdges;
				fps.call(G, delEdges);
				removeEdges(G, delEdges);

				System::usedRealTime(T);
				time += System::usedRealTime(T);
				AssertThat(pm.isPlanar(G), Equals(true));
			}
		}
		cout << endl << "      average time was " << time/STEPS/3 << "ms" << endl;
	});

	it("detects non-planarity in small complete graphs", [&](){
		for(int i = 5; i < 50; i++) {
			completeGraph(G, i);
			AssertThat(pm.isPlanar(G), Equals(i < 5));
		}
	});
});
}

void testSubgraph(Graph &graph, PlanarSubgraphModule &sm, PlanarityModule &tester, bool assertMaximality, bool weightEdges)
{
	EdgeArray<int> costs;
	edge mustHaveEdge = nullptr;

	if(weightEdges) {
		costs.init(graph);

		for(edge e : graph.edges) {
			costs[e] = 1;
		}

		mustHaveEdge = graph.chooseEdge();
		costs[mustHaveEdge] = 1 + graph.numberOfEdges();
	}


	List<edge> removedEdges;
	if(weightEdges) {
		sm.call(graph, costs, removedEdges);
	} else {
		sm.call(graph, removedEdges);
	}

	cout << std::endl << "      removed " << removedEdges.size() << " edges" << std::endl;
	bool connected = isConnected(graph);

	Graph::HiddenEdgeSetHandle handle = graph.newHiddenEdgeSet();
	List<edge> hiddenEdges;
	for(edge e : removedEdges) {
		graph.hideEdge(handle, e);
		hiddenEdges.pushBack(e);
		AssertThat(e, Is().Not().EqualTo(mustHaveEdge));
	}

	AssertThat(isConnected(graph), Equals(connected));
	AssertThat(tester.isPlanar(graph), Equals(true));

	if(assertMaximality) {
		ListConstIterator<edge> it = removedEdges.begin();
		for(edge e : hiddenEdges) {
			graph.restoreEdge(handle, e);
			AssertThat(tester.isPlanar(graph), Equals(false));
			graph.hideEdge(handle, e);
			it = it.succ();
		}
	}
}

void testNonPlanarCore()
{
	for_each_graph_it("returns a simple core", {"north/g.41.26.gml", "north/g.73.8.gml"},
		[&](Graph &graph, const string &filename){
			makeBiconnected(graph);
			NonPlanarCore npc(graph);
			const Graph &core = npc.core();
			AssertThat(isSimpleUndirected(core), IsTrue());

			for(edge e : core.edges) {
				AssertThat(npc.cost(e), IsGreaterThan(0));

				if(!npc.isVirtual(e)) {
					AssertThat(npc.realEdge(e), Is().Not().EqualTo((void*) nullptr));
				}
			}
		});

	it("works on a minimal previously failing instance (2 x K5)", []() {
		Graph graph;

		node s = graph.newNode();
		node t = graph.newNode();
		graph.newEdge(s, t);

		node v = graph.newNode();
		graph.newEdge(s, v);
		graph.newEdge(v, t);

		for (int k = 0; k < 2; k++) {
			List<node> nodes;
			nodes.pushBack(s);
			nodes.pushBack(t);

			for (int i = 0; i < 3; i++) {
				nodes.pushBack(graph.newNode());
			}

			for (node v : nodes) {
				for (node w : nodes) {
					if (v->index() < w->index() && (v != s || w != t)) {
						graph.newEdge(v, w);
					}
				}
			}
		}

		NonPlanarCore npc(graph);
		const Graph &core = npc.core();
		AssertThat(isSimpleUndirected(core), IsTrue());
		AssertThat(core.numberOfNodes(), Equals(graph.numberOfNodes() - 1));
		AssertThat(core.numberOfEdges(), Equals(graph.numberOfEdges() - 2));

		for (edge e : core.edges) {
			node v = npc.original(e->source());
			node w = npc.original(e->target());
			if((v == s && w == t) || (v == t && w == s)) {
				AssertThat(npc.cost(e), Equals(2));
			} else {
				AssertThat(npc.cost(e), Equals(1));
			}
		}
	});

	it("eliminates self-loops", [](){
		Graph graph;
		completeGraph(graph, 5);

		for(node v : graph.nodes) {
			graph.newEdge(v, v);
		}

		NonPlanarCore npc(graph);
		const Graph &core = npc.core();
		AssertThat(isSimpleUndirected(core), IsTrue());
		AssertThat(core.numberOfNodes(), Equals(graph.numberOfNodes()));
		AssertThat(core.numberOfEdges(), Equals(10));
	});
}

go_bandit([](){
	describe("Planarity tests", [](){
		BoothLueker bl;
		describeModule("Booth-Lueker", bl);
		BoyerMyrvold bm;
		describeModule("Boyer-Myrvold", bm);
	});

	describe("Planar Subgraphs", [](){
		describe("Boyer-Myrvold-Subgraph", [](){
			bool maximizeSubgraph = true;
			BoothLueker bl;
			BoyerMyrvoldSubgraph bms(maximizeSubgraph, 10, 1);
			minstd_rand rng(42);

			for(int n = 5; n < 30; n++) {
				it(string("works on a K" + to_string(n)).c_str(), [&](){
					Graph graph;
					completeGraph(graph, n);
					randomizeAdjLists(graph, rng);
					testSubgraph(graph, bms, bl, maximizeSubgraph, false);
				});
			}

			for(int n = 5; n < 30; n++) {
				int m = n*(n-1)/3;
				it(string("works on a random graph with " + to_string(n) + " nodes and " + to_string(m) + " edges").c_str(), [&](){
					Graph graph;
					randomGraph(graph, n, m);
					randomizeAdjLists(graph, rng);
					testSubgraph(graph, bms, bl, maximizeSubgraph, false);
				});
			}

			for(int n = 5; n < 20; n++) {
				it(string("works on a K" + to_string(n) + " with weighted edges").c_str(), [&](){
					Graph graph;
					completeGraph(graph, n);
					randomizeAdjLists(graph, rng);
					BoyerMyrvoldSubgraph myBms(false, 1, 0);
					testSubgraph(graph, myBms, bl, false, true);
				});
			}

			vector<string> instaces = {
				"north/g.61.11.gml",
				"rome/grafo3703.45.lgr.gml.pun",
				"rome/grafo5745.50.lgr.gml.pun"
			};

			for(int i = 0; i < 2; i++) {
				string tags = "";
				bool weighted = i;
				if(weighted) { tags += " weighted,"; }

				for_each_graph_it(("works on a" + tags + " formerly failing instance").c_str(), instaces,
						[&](Graph &graph, const string &filename){
							int randomness = weighted ? 0 : 1;
							BoyerMyrvoldSubgraph myBms(maximizeSubgraph, 1, randomness);
							testSubgraph(graph, myBms, bl, maximizeSubgraph, weighted);
						});
			}
		});
	});

	describe("NonPlanarCore", [](){
		testNonPlanarCore();
	});
});
