#include <bandit/bandit.h>

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/graph_generators.h>
#include <ogdf/basic/simple_graph_alg.h>

using namespace ogdf;
using namespace bandit;

/**
 * Assert that isAcylic()/isAcyclicUndirected() returns the correct value and
 * that the list of collected backedges is filled correctly. For cyclic graphs
 * assert that removing all backedges makes the graph acyclic but maintains
 * connectivity.
 *
 * @param G is the graph to be tested.
 * @param directed sets whether isAcyclic() or isAcyclicUndirected() is tested.
 * @param expected is the expected result of the function call.
 */
void isAcyclicAssert(Graph &G, bool directed, bool expected)
{
	List<edge> backedges;
	bool result = directed ?
		          isAcyclic(G, backedges) : isAcyclicUndirected(G, backedges);

	if (expected) {
		AssertThat(result, IsTrue());
		AssertThat(backedges.empty(), IsTrue());
	} else {
		AssertThat(result, IsFalse());
		AssertThat(backedges.size(), IsGreaterThan(0));
		AssertThat(backedges.size(), IsLessThan(G.numberOfEdges() + 1));

		bool connected = isConnected(G);

		for (edge e : backedges) {
			G.delEdge(e);
		}

		result = directed ?
		         isAcyclic(G, backedges) : isAcyclicUndirected(G, backedges);
		AssertThat(result, IsTrue());
		AssertThat(backedges.empty(), IsTrue());
		AssertThat(isConnected(G), Equals(connected));
	}
}

/**
 * Perform tests for isAcylic() or isAcyclicUndirected().
 *
 * @param directed sets whether isAcyclic() or isAcyclicUndirected() is tested.
 */
void describeIsAcyclic(bool directed)
{
	Graph G;

	before_each([&](){
		G.clear();
	});

	it("works on an empty graph", [&](){
		isAcyclicAssert(G, directed, true);
	});

	it("works on a graph with a single node", [&](){
		G.newNode();
		isAcyclicAssert(G, directed, true);
	});

	it("works on a graph with a self-loop", [&](){
		node v1 = G.newNode();
		G.newEdge(v1, v1);

		isAcyclicAssert(G, directed, false);
	});

	it("works on a graph with parallel edges", [&](){
		node v1 = G.newNode();
		node v2 = G.newNode();
		G.newEdge(v1, v2);
		G.newEdge(v2, v1);

		isAcyclicAssert(G, directed, false);
	});

	it("works on an acylic graph", [&](){
		node v1 = G.newNode();
		node v2 = G.newNode();
		node v3 = G.newNode();
		G.newEdge(v1, v2);
		G.newEdge(v1, v3);

		isAcyclicAssert(G, directed, true);
	});

	it("works on a cyclic graph", [&](){
		node v1 = G.newNode();
		node v2 = G.newNode();
		node v3 = G.newNode();
		G.newEdge(v1, v2);
		G.newEdge(v2, v3);
		G.newEdge(v3, v1);

		isAcyclicAssert(G, directed, false);
	});

	it("works on a disconnected acyclic graph", [&](){
		G.newNode();
		node v2 = G.newNode();
		node v3 = G.newNode();
		node v4 = G.newNode();
		G.newEdge(v2, v3);
		G.newEdge(v2, v4);

		isAcyclicAssert(G, directed, true);
	});

	it("works on a disconnected cyclic graph", [&](){
		G.newNode();
		node v2 = G.newNode();
		node v3 = G.newNode();
		node v4 = G.newNode();
		G.newEdge(v2, v3);
		G.newEdge(v3, v4);
		G.newEdge(v4, v2);

		isAcyclicAssert(G, directed, false);
	});

	it("works on an acyclic graph requiring multiple dfs starts if directed", [&](){
		node v1 = G.newNode();
		node v2 = G.newNode();
		node v3 = G.newNode();
		node v4 = G.newNode();
		G.newEdge(v1, v2);
		G.newEdge(v2, v3);
		G.newEdge(v4, v2);

		isAcyclicAssert(G, directed, true);
	});

	it("works on a cyclic graph requiring multiple dfs starts if directed", [&](){
		node v1 = G.newNode();
		node v2 = G.newNode();
		node v3 = G.newNode();
		node v4 = G.newNode();
		G.newEdge(v1, v2);
		G.newEdge(v2, v3);
		G.newEdge(v3, v1);
		G.newEdge(v4, v3);

		isAcyclicAssert(G, directed, false);
	});

	it("works on a directed acyclic but undirected cyclic graph", [&](){
		node v1 = G.newNode();
		node v2 = G.newNode();
		node v3 = G.newNode();
		G.newEdge(v1, v2);
		G.newEdge(v1, v3);
		G.newEdge(v2, v3);

		isAcyclicAssert(G, directed, directed);
	});

	it("works on an extremely large acyclic graph", [&](){
		randomTree(G, 125000, 1, 0);
		isAcyclicAssert(G, directed, true);
	});

	it("works on an extremely large cyclic graph", [&](){
		randomBiconnectedGraph(G, 125000, 250000);
		isAcyclicAssert(G, directed, false);
	});
}

go_bandit([]() {
	describe("Simple Graph Algorithms", [](){
		describe("isBiconnected", [](){
			Graph G;
			node cutVertex;

			node v1 = G.newNode();
			node v2 = G.newNode();
			node v3 = G.newNode();
			G.newEdge(v1, v2);

			it("works on a disconnected graph", [&](){
				AssertThat(isBiconnected(G, cutVertex), IsFalse());
				AssertThat(cutVertex, Equals(nullptr));
			});

			G.newEdge(v1, v3);

			it("works on a connected but not biconnected graph", [&](){
				AssertThat(isBiconnected(G, cutVertex), IsFalse());
				AssertThat(cutVertex, Equals(v1));
			});

			G.newEdge(v2, v3);

			it("works on a simple biconnected graph", [&](){
				AssertThat(isBiconnected(G, cutVertex), IsTrue());
				AssertThat(cutVertex, Equals(nullptr));
			});

			it("works on an extremely large graph", [&](){
				randomBiconnectedGraph(G, 250000, 500000);
				AssertThat(isBiconnected(G), IsTrue());
			});
		});

		describe("makeBiconnected", [](){
			Graph G;
			List<edge> added;

			before_each([&](){
				G.clear();
				added.clear();
			});

			it("works on a disconnected graph", [&](){
				node v1 = G.newNode();
				node v2 = G.newNode();
				G.newNode(); // node without edge
				G.newEdge(v1, v2);

				makeBiconnected(G, added);
				AssertThat(isBiconnected(G), IsTrue());
				AssertThat(added.size(), Equals(2));
			});

			it("works on a connected but not biconnected graph", [&](){
				node v1 = G.newNode();
				node v2 = G.newNode();
				node v3 = G.newNode();
				G.newEdge(v1, v2);
				G.newEdge(v1, v3);

				makeBiconnected(G, added);
				AssertThat(isBiconnected(G), IsTrue());
				AssertThat(added.size(), Equals(1));
			});

			it("works on a simple biconnected graph", [&](){
				randomBiconnectedGraph(G, 10, 20);
				AssertThat(isBiconnected(G), IsTrue());

				makeBiconnected(G, added);
				AssertThat(isBiconnected(G), IsTrue());
				AssertThat(added.empty(), IsTrue());
			});

			it("works on an extremely large graph", [&](){
				emptyGraph(G, 250000);
				AssertThat(isBiconnected(G), IsFalse());

				// A graph with n nodes needs at least n edges to be biconnected
				makeBiconnected(G, added);
				AssertThat(isBiconnected(G), IsTrue());
				AssertThat(added.size(), IsGreaterThan(250000));
			});
		});

		describe("isAcyclic", [](){
			describeIsAcyclic(true);
		});

		describe("isAcyclicUndirected", [](){
			describeIsAcyclic(false);
		});
	});
});
