//*********************************************************
//  Bandit tests for strong component algorithms
//
//  Author: Tilo Wiedera
//*********************************************************

#include <bandit/bandit.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/basic/graph_generators.h>

using namespace bandit;
using namespace ogdf;

/** Checks wheter a path from source to target
 *  does exist by traversing the graph.
 */
bool pathExists(const Graph &graph, const node source, const node target)
{
	OGDF_ASSERT(source != target);
	OGDF_ASSERT(source->graphOf() == &graph);
	OGDF_ASSERT(target->graphOf() == &graph);

	List<node> queue;
	NodeArray<bool> visited(graph, false);
	visited[source] = true;
	queue.pushBack(source);

	bool result = false;
	while(!queue.empty() && !result) {
		node v = queue.popFrontRet();
		for(adjEntry adj : v->adjEntries) {
			node w = adj->theEdge()->target();
			if(!visited[w]) {
				result |= w == target;
				visited[w] = true;
				queue.pushBack(w);
			}
		}
	}

	return result;
}

go_bandit([](){
	describe("strong components", [](){
		for(int n = 0; n < 75; n++) {
			it(string("works on a random graph of size " + to_string(n)), [&](){
				Graph graph;
				randomDiGraph(graph, n, randomDouble(0, 1));

				NodeArray<int> components(graph);
				int nComponents = strongComponents(graph, components);
				List<node> nodes;
				graph.allNodes(nodes);
				for(node v = graph.firstNode(); v; v = v->succ()) {
					for(node w = v->succ(); w; w = w->succ()) {
						AssertThat(components[v], IsGreaterThan(-1) && IsLessThan(nComponents));
						AssertThat(components[w], IsGreaterThan(-1) && IsLessThan(nComponents));
						if(components[v] == components[w]) {
							AssertThat(pathExists(graph, v, w), IsTrue());
							AssertThat(pathExists(graph, w, v), IsTrue());
						} else {
							AssertThat(pathExists(graph, w, v) && pathExists(graph, v, w), IsFalse());
						}
					}
				}
			});
		}

		it("works on a predefined graph with overlapping circles", [](){
			Graph graph;
			emptyGraph(graph, 8);
			List<node> nodes;
			graph.allNodes(nodes);
			graph.newEdge(*nodes.get(2), *nodes.get(5));
			graph.newEdge(*nodes.get(3), *nodes.get(6));
			graph.newEdge(*nodes.get(4), *nodes.get(7));
			graph.newEdge(*nodes.get(5), *nodes.get(4));
			graph.newEdge(*nodes.get(6), *nodes.get(5));
			graph.newEdge(*nodes.get(6), *nodes.get(1));
			graph.newEdge(*nodes.get(7), *nodes.get(2));
			graph.newEdge(*nodes.get(7), *nodes.get(3));
			graph.newEdge(*nodes.get(7), *nodes.get(6));

			NodeArray<int> components(graph);
			int nComponents = strongComponents(graph, components);
			AssertThat(nComponents, Equals(3));
			AssertThat(components[*nodes.get(0)], !Equals(components[*nodes.get(1)]));
			AssertThat(components[*nodes.get(0)], !Equals(components[*nodes.get(2)]));
			AssertThat(components[*nodes.get(1)], !Equals(components[*nodes.get(2)]));
			AssertThat(components[*nodes.get(2)],  Equals(components[*nodes.get(7)]));
			AssertThat(components[*nodes.get(3)],  Equals(components[*nodes.get(7)]));
			AssertThat(components[*nodes.get(4)],  Equals(components[*nodes.get(7)]));
			AssertThat(components[*nodes.get(5)],  Equals(components[*nodes.get(7)]));
			AssertThat(components[*nodes.get(6)],  Equals(components[*nodes.get(7)]));


		});
	});
});
