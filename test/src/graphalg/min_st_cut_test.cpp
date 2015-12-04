#include <bandit/bandit.h>

#include <ogdf/graphalg/MaxFlowEdmondsKarp.h>
#include <ogdf/graphalg/MinSTCut.h>
#include <ogdf/basic/graph_generators.h>

using namespace ogdf;
using namespace bandit;
using std::string;

template<typename T>
void describeSuite(const string &name) {
describe(string("MinSTCut<" + name + ">"), [](){
	it("works on a simple example", [](){
		Graph graph;
		node s = graph.newNode();
		node t = graph.newNode();
		node v1 = graph.newNode();
		node v2 = graph.newNode();

		graph.newEdge(s, v1);
		graph.newEdge(v2, t);

		EdgeArray<T> weights(graph, 4);
		weights[graph.newEdge(s, v2)] = 1;
		weights[graph.newEdge(v1, t)] = 2;

		EdgeArray<T> flow;
		MaxFlowEdmondsKarp<T> mfek(graph);
		mfek.computeFlow(weights, s, t, flow);
		MinSTCut<T> minSTCut;
		minSTCut.call(weights, flow, s, t);

		AssertThat(minSTCut.isInFrontCut(s), Equals(true));
		AssertThat(minSTCut.isInFrontCut(v1), Equals(true));
		AssertThat(minSTCut.isInBackCut(t), Equals(true));
		AssertThat(minSTCut.isInBackCut(v2), Equals(true));
	});

	it("works on a more complex example", [](){
		Graph graph;
		for(int i = 0; i < 8; i++) {
			graph.newNode();
		}
		List<node> nodes;
		graph.allNodes(nodes);
		EdgeArray<T> weights(graph);
		weights[graph.newEdge(*nodes.get(0), *nodes.get(1))] = 16;
		weights[graph.newEdge(*nodes.get(0), *nodes.get(2))] = 13;
		weights[graph.newEdge(*nodes.get(1), *nodes.get(2))] = 10;
		weights[graph.newEdge(*nodes.get(1), *nodes.get(3))] = 12;
		weights[graph.newEdge(*nodes.get(2), *nodes.get(1))] =  4;
		weights[graph.newEdge(*nodes.get(2), *nodes.get(4))] = 14;
		weights[graph.newEdge(*nodes.get(3), *nodes.get(2))] =  9;
		weights[graph.newEdge(*nodes.get(3), *nodes.get(5))] = 20;
		weights[graph.newEdge(*nodes.get(4), *nodes.get(3))] =  7;
		weights[graph.newEdge(*nodes.get(4), *nodes.get(5))] =  4;

		weights[graph.newEdge(*nodes.get(5), *nodes.get(6))] =  100;
		weights[graph.newEdge(*nodes.get(7), *nodes.get(5))] =  100;

		EdgeArray<T> flow;
		MaxFlowEdmondsKarp<T> mfek(graph);
		mfek.computeFlow(weights, *nodes.get(0), *nodes.get(5), flow);
		MinSTCut<T> minSTCut;
		minSTCut.call(weights, flow, *nodes.get(0), *nodes.get(5));

		AssertThat(minSTCut.isInFrontCut(*nodes.get(0)), Equals(true));
		AssertThat(minSTCut.isInFrontCut(*nodes.get(1)), Equals(true));
		AssertThat(minSTCut.isInFrontCut(*nodes.get(2)), Equals(true));
		AssertThat(minSTCut.isInFrontCut(*nodes.get(4)), Equals(true));

		AssertThat(minSTCut.isInBackCut(*nodes.get(3)), Equals(true));
		AssertThat(minSTCut.isInBackCut(*nodes.get(5)), Equals(true));

		AssertThat(minSTCut.isInFrontCut(*nodes.get(6)), Equals(false));
		AssertThat(minSTCut.isInBackCut(*nodes.get(6)), Equals(false));
		AssertThat(minSTCut.isOfType(*nodes.get(6), MinSTCut<T>::NO_CUT), Equals(true));
		AssertThat(minSTCut.isInBackCut(*nodes.get(7)), Equals(true));
		AssertThat(minSTCut.isInFrontCut(*nodes.get(7)), Equals(false));
	});
});
}

go_bandit([](){
	describe("MinSTCut class", []() {
		describeSuite<int>("int");
		describeSuite<double>("double");
	});
});
