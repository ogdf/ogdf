#include <bandit/bandit.h>
#include <resources.h>

#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/graph_generators.h>
#include <ogdf/planarity/EmbedderMaxFace.h>
#include <ogdf/planarity/EmbedderMaxFaceLayers.h>
#include <ogdf/planarity/EmbedderMinDepth.h>
#include <ogdf/planarity/EmbedderMinDepthMaxFace.h>
#include <ogdf/planarity/EmbedderMinDepthMaxFaceLayers.h>
#include <ogdf/planarity/EmbedderMinDepthPiTa.h>
#include <ogdf/planarity/EmbedderOptimalFlexDraw.h>
#include <ogdf/planarity/SimpleEmbedder.h>

using namespace bandit;
using namespace ogdf;

constexpr int numberOfNodes = 42;

void validateCopy(const Graph &graph, const GraphCopy &copy) {
	AssertThat(graph.numberOfNodes(), Equals(copy.numberOfNodes()));
	AssertThat(graph.numberOfEdges(), Equals(copy.numberOfEdges()));

	for(node v : copy.nodes) {
		AssertThat(copy.isDummy(v), IsFalse());
	}

	for(edge e : copy.edges) {
		AssertThat(copy.isDummy(e), IsFalse());
		edge f = copy.original(e);
		AssertThat(f->source(), Equals(copy.original(e->source())));
		AssertThat(f->target(), Equals(copy.original(e->target())));
	}
}

void shuffleEmbedding(Graph &graph) {
	for(node v : graph.nodes) {
		for(adjEntry adj : v->adjEntries) {
			graph.swapAdjEdges(adj, randomNumber(0, 1) ? v->firstAdj() : v->lastAdj());
		}
	}
}

void testEmbedder(EmbedderModule &embedder, const Graph &graph, bool repeat = true) {
	GraphCopy copy(graph);
	if(repeat) {
		shuffleEmbedding(copy);
	}
	// initialize adjExternal with a corrupt value
	adjEntry adjExternal = graph.firstNode()->firstAdj();

	embedder(copy, adjExternal);

	validateCopy(graph, copy);
	AssertThat(adjExternal, !IsNull());
#ifdef OGDF_DEBUG
	AssertThat(adjExternal->graphOf(), Equals(&copy));
#endif
	AssertThat(copy.representsCombEmbedding(), IsTrue());

	// test planarly embedded input
	if(repeat) {
		testEmbedder(embedder, copy, false);
	}
}

void describeEmbedder(const string &title, EmbedderModule &embedder) {
	describe(title, [&]() {
		Graph graph;

		before_each([&]() {
			graph.clear();
		});

		it("works on a tree", [&]() {
			randomTree(graph, numberOfNodes);
			testEmbedder(embedder, graph);
		});

		it("works on a fully triangulated graph", [&]() {
			planarConnectedGraph(graph, numberOfNodes, 3*numberOfNodes - 6);
			testEmbedder(embedder, graph);
		});

		for(int n = numberOfNodes; n < numberOfNodes + 10; n++) {
			it("works on a random planar graph containing " + to_string(n) + " nodes", [&]() {
				planarTriconnectedGraph(graph, n, randomNumber(int(1.5*n), 3*n-6));
				testEmbedder(embedder, graph);
			});
		}
#ifdef OGDF_DEBUG
		it("fails on a K5", [&]() {
			adjEntry adjExternal;
			completeGraph(graph, 5);
			AssertThrows(PreconditionViolatedException, embedder(graph, adjExternal));
		});

		it("fails on a random non-planar graph", [&]() {
			adjEntry adjExternal;
			randomGraph(graph, numberOfNodes, 3*numberOfNodes-5);
			AssertThrows(PreconditionViolatedException, embedder(graph, adjExternal));
		});
#endif
	});
}

template<typename EmbedderType>
void describeEmbedder(const string &title) {
	EmbedderType embedder;
	describeEmbedder(title, embedder);
}

template<>
void describeEmbedder<EmbedderMinDepthPiTa>(const string &title) {
	EmbedderMinDepthPiTa embedder;
	bool extendedDD = embedder.useExtendedDepthDefinition();
	describeEmbedder(title + " [extendedDD=" + to_string(extendedDD) + "]", embedder);
	embedder.useExtendedDepthDefinition(!extendedDD);
	describeEmbedder(title + " [extendedDD=" + to_string(!extendedDD) + "]", embedder);
}

template<>
void describeEmbedder<EmbedderOptimalFlexDraw>(const string &title) {
	describe(title, [](){
		EmbedderOptimalFlexDraw embedder;
		describeEmbedder("Non-Weighted Version", embedder);

		describe("Weighted Edges", [&](){
			it("works on a random graph", [&](){
				Graph graph;
				planarConnectedGraph(graph, numberOfNodes, 2*numberOfNodes);
				EdgeArray<int> costs(graph);

				for(edge e : graph.edges) {
					costs[e] = randomNumber(1, numberOfNodes);
				}

				testEmbedder(embedder, graph);
			});
		});
	});
}

go_bandit([]() {
	describe("Embedders", []() {
		describeEmbedder<EmbedderMaxFace>("EmbedderMaxFace");
		describeEmbedder<EmbedderMaxFaceLayers>("EmbedderMaxFaceLayers");
		describeEmbedder<EmbedderMinDepth>("EmbedderMinDepth");
		describeEmbedder<EmbedderMinDepthMaxFace>("EmbedderMinDepthMaxFace");
		describeEmbedder<EmbedderMinDepthMaxFaceLayers>("EmbedderMinDepthMaxFaceLayers");
#if 0
		describeEmbedder<EmbedderMinDepthPiTa>("EmbedderMinDepthPiTa");
		describeEmbedder<EmbedderOptimalFlexDraw>("EmbedderOptimalFlexDraw");
#endif
		describeEmbedder<SimpleEmbedder>("SimpleEmbedder");
	});
});
