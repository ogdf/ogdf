/** \file
 * \brief Tests for the basic graph class
 *
 * \author Tilo Wiedera
 *
 * \par License:
 * This file is part of the Open myGraph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.txt in the root directory of the OGDF installation for details.
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
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/

#include <bandit/bandit.h>

#include <ogdf/basic/Graph.h>
#include <resources.h>

using namespace ogdf;
using namespace bandit;

/**
 * Returns an arbitrary edge where both nodes have at least \c minDegree incident edges.
 *
 * @param graph graph to investigate
 * @param minDegree minimal number of incident edges
 * @return the chosen edge or \c nullptr if none could be found
 */
edge chooseEdge(const Graph &graph, int minDegree) {
	for(edge e : graph.edges) {

		if(e->source()->degree() >= minDegree && e->target()->degree() >= minDegree) {
			return e;
		}
	}

	return nullptr;
}

/**
 * Returns an arbitrary node with at least \c minDegree incident edges.
 *
 * @param graph graph to investigate
 * @param minDegree minimal number of incident edges
 * @return the chosen node or \c nullptr if none could be found
 */
node chooseNode(const Graph &graph, int minDegree) {
	for(node v : graph.nodes) {

		if(v->degree() >= minDegree) {
			return v;
		}
	}

	return nullptr;
}

/**
 * Returns an arbitrary node which does not equal \c v.
 *
 * @param graph graph to investigate
 * @param v the node to exclude from selection
 * @return the chosen node or \c nullptr if none could be found
 */
node chooseNode(const Graph &graph, node v) {
	for(node w : graph.nodes) {

		if(w != v) {
			return w;
		}
	}

	return nullptr;
}

go_bandit([](){
describe("Graph Class", [](){
	std::vector<std::string> files = {"rome/grafo3703.45.lgr.gml.pun", "rome/grafo5745.50.lgr.gml.pun", "north/g.41.26.gml", "north/g.61.11.gml", "north/g.73.8.gml"};

	it("is initialized correctly", [](){
		Graph graph;

		AssertThat(graph.empty(), IsTrue());
		AssertThat(graph.numberOfNodes(), Equals(0));
		AssertThat(graph.numberOfEdges(), Equals(0));
		AssertThat(graph.maxNodeIndex(), IsLessThan(0));
		AssertThat(graph.maxEdgeIndex(), IsLessThan(0));
		AssertThat(graph.maxAdjEntryIndex(), IsLessThan(0));
		AssertThat(graph.nodeArrayTableSize(), IsGreaterThan(0));
		AssertThat(graph.edgeArrayTableSize(), IsGreaterThan(0));
		AssertThat(graph.adjEntryArrayTableSize(), IsGreaterThan(0));
		AssertThat(graph.firstNode(), Equals((void*) nullptr));
		AssertThat(graph.lastNode(), Equals((void*) nullptr));
		AssertThat(graph.firstEdge(), Equals((void*) nullptr));
		AssertThat(graph.lastEdge(), Equals((void*) nullptr));
	});

	for_each_graph_it("finds a existing edges", files, [](Graph &graph, const string file){
		edge e = graph.chooseEdge();
		AssertThat(graph.searchEdge(e->source(), e->target()), Equals(e));
	});

	for_each_graph_it("finds a reverse edge", files, [](Graph &graph, const string file){
		edge e = graph.chooseEdge();
		AssertThat(graph.searchEdge(e->target(), e->source()), Equals(e));
	});

	for_each_graph_it("does not find non-existent edges", files, [](Graph &graph, const string file){
		edge e = graph.chooseEdge();
		graph.delEdge(e);
		AssertThat(graph.searchEdge(e->source(), e->target()), Equals((void*) nullptr));
	});

	for_each_graph_it("can be assigned", files, [](Graph &graph, const string file){
		int m = graph.numberOfEdges();

		int *degreeCounter = new int[m];

		for(int i = 0; i < m; i++) {
			degreeCounter[i] = 0;
		}

		for(node v : graph.nodes) {
			degreeCounter[v->degree()]++;
		}

		Graph copy = graph;

		AssertThat(copy.numberOfNodes(), Equals(graph.numberOfNodes()));
		AssertThat(copy.numberOfEdges(), Equals(m));

		for(node v : copy.nodes) {
			degreeCounter[v->degree()]--;
		}

		for(node v : graph.nodes) {
			AssertThat(degreeCounter[v->degree()], Equals(0));
		}

		delete[] degreeCounter;
	});

	it("adds nodes", [](){
		Graph graph;
		const int numberOfNodes = 100;

		for(int i = 0; i < numberOfNodes; i++) {
			graph.newNode();
		}

		AssertThat(graph.empty(), IsFalse());
		AssertThat(graph.numberOfNodes(), Equals(numberOfNodes));
		AssertThat(graph.numberOfEdges(), Equals(0));
		AssertThat(graph.maxNodeIndex(), IsGreaterThan(numberOfNodes - 2));
		AssertThat(graph.firstNode(), Is().Not().EqualTo((void*) nullptr));
		AssertThat(graph.lastNode(), Is().Not().EqualTo((void*) nullptr));

		int maxIndex = graph.maxNodeIndex();
		bool *visited = new bool[maxIndex + 1];

		for(int i = 0; i <= maxIndex; i++) {
			visited[i] = false;
		}

		int count = 0;

		for(node v : graph.nodes) {
			int index = v->index();
			AssertThat(index, IsGreaterThan(-1));
			AssertThat(index, IsLessThan(maxIndex + 1));
			AssertThat(visited[index], IsFalse());
			visited[index] = true;
			count++;
		}

		AssertThat(count, Equals(numberOfNodes));

		delete[] visited;
	});

	it("adds edges", [](){
		Graph graph;

		for(int i = 0; i < 100; i++) {
			graph.newNode();
		}

		int count = 0;

		for(node v : graph.nodes) {
			for(node w : graph.nodes) {
				if((v->index() + w->index()) % 3 == 0) {
					graph.newEdge(v, w);
					count++;
				}
			}
		}

		AssertThat(graph.numberOfEdges(), Equals(count));
		AssertThat(graph.maxEdgeIndex(), IsGreaterThan(count - 2));
		AssertThat(graph.maxAdjEntryIndex(), IsGreaterThan(count - 2));
		AssertThat(graph.firstEdge(), Is().Not().EqualTo((void*) nullptr));
		AssertThat(graph.lastEdge(), Is().Not().EqualTo((void*) nullptr));

		int maxIndex = graph.maxEdgeIndex();
		bool *visited = new bool[maxIndex + 1];

		for(int i = 0; i <= maxIndex; i++) {
			visited[i] = false;
		}

		int iterCount = 0;

		for(edge e : graph.edges) {
			int index = e->index();
			AssertThat(index, IsGreaterThan(-1));
			AssertThat(index, IsLessThan(maxIndex + 1));
			AssertThat(visited[index], IsFalse());
			visited[index] = true;
			iterCount++;
		}

		AssertThat(iterCount, Equals(count));

		delete[] visited;
	});

	for_each_graph_it("removes a node", files, [](Graph &graph, const string file){
		int n = graph.numberOfNodes();
		int m = graph.numberOfEdges();

		node v = graph.chooseNode();
		int deg = v->degree();

		graph.delNode(v);

		AssertThat(graph.numberOfNodes(), Equals(n - 1));
		AssertThat(graph.numberOfEdges(), Equals(m - deg));
	});

	for_each_graph_it("removes an edge", files, [](Graph &graph, const string file){
		int n = graph.numberOfNodes();
		int m = graph.numberOfEdges();

		edge e = graph.chooseEdge();
		node s = e->source();
		node t = e->target();

		graph.delEdge(e);

		AssertThat(graph.searchEdge(s, t), Equals((void*) nullptr));
		AssertThat(graph.numberOfNodes(), Equals(n));
		AssertThat(graph.numberOfEdges(), Equals(m - 1));
	});

	for_each_graph_it("can be cleared", files, [](Graph &graph, const string file){
		graph.clear();

		AssertThat(graph.empty(), IsTrue());
		AssertThat(graph.numberOfNodes(), Equals(0));
		AssertThat(graph.numberOfEdges(), Equals(0));
	});

	for_each_graph_it("hides an edge and restores it", files, [](Graph &graph, const string file){
		int n = graph.numberOfNodes();
		int m = graph.numberOfEdges();

		edge e = graph.chooseEdge();
		Graph::HiddenEdgeSetHandle handle = graph.newHiddenEdgeSet();
		graph.hideEdge(handle, e);

		AssertThat(graph.numberOfNodes(), Equals(n));
		AssertThat(graph.numberOfEdges(), Equals(m - 1));
		AssertThat(graph.searchEdge(e->source(), e->target()), Equals((void*) nullptr));

		graph.restoreEdge(handle, e);

		AssertThat(graph.numberOfEdges(), Equals(m));
		AssertThat(graph.searchEdge(e->source(), e->target()), Equals(e));
	});

	for_each_graph_it("restores all hidden edges", files, [](Graph &graph, const string file){
		int m = graph.numberOfEdges();

		Graph::HiddenEdgeSetHandle handle = graph.newHiddenEdgeSet();

		for(int i = 0; i < m / 2; i++) {
			graph.hideEdge(handle, graph.chooseEdge());
		}

		AssertThat(graph.numberOfEdges(), Equals(m - m / 2));
		graph.restoreEdges(handle);
		AssertThat(graph.numberOfEdges(), Equals(m));
	});

	for_each_graph_it("hides all edges across 10 sets", files, [](Graph &graph, const string file){
		int m = graph.numberOfEdges();
		int maxIndex = graph.maxNodeIndex();

		int *inDeg = new int[maxIndex + 1];
		int *outDeg = new int[maxIndex + 1];

		for(node v : graph.nodes) {
			inDeg[v->index()] = v->indeg();
			outDeg[v->index()] = v->outdeg();
		}

		List<Graph::HiddenEdgeSetHandle> handles;

		for(int i = 0; i < 10; i++) {
			handles.pushFront(graph.newHiddenEdgeSet());

			for(int k = 0; k < m / 10; k++) {
				graph.hideEdge(handles.front(), graph.chooseEdge());
			}
		}

		handles.permute();

		while(graph.numberOfEdges() > 0) {
			graph.hideEdge(handles.front(), graph.chooseEdge());
		}

		for(node v : graph.nodes) {
			AssertThat(v->indeg(), Equals(0));
			AssertThat(v->outdeg(), Equals(0));
		}

		for(Graph::HiddenEdgeSetHandle handle : handles) {
			graph.restoreEdges(handle);
		}

		AssertThat(graph.numberOfEdges(), Equals(m));

		for(node v : graph.nodes) {
			AssertThat(v->indeg(), Equals(inDeg[v->index()]));
			AssertThat(v->outdeg(), Equals(outDeg[v->index()]));
		}

		delete[] inDeg;
		delete[] outDeg;
	});

	for_each_graph_it("reverses an edge", files, [](Graph &graph, const string file){
		edge e = chooseEdge(graph, 5);
		node s = e->source();
		node t = e->target();

		int inT = t->indeg();
		int outT = t->outdeg();

		int inS = s->indeg();
		int outS = s->outdeg();

		graph.reverseEdge(e);

		AssertThat(e->source(), Equals(t));
		AssertThat(e->target(), Equals(s));
		AssertThat(e->source()->degree(), Equals(inT + outT));
		AssertThat(e->target()->degree(), Equals(inS + outS));
		AssertThat(e->source()->indeg(), Equals(inT - 1));
		AssertThat(e->source()->outdeg(), Equals(outT + 1));
	});

	for_each_graph_it("reverses all edges", files, [](Graph &graph, const string file){
		int maxIndex = graph.maxEdgeIndex();
		node *sources = new node[maxIndex + 1];
		node *targets = new node[maxIndex + 1];

		for(int i = 0; i <= maxIndex; i++) {
			sources[i] = targets[i] = nullptr;
		}

		for(edge e : graph.edges) {
			sources[e->index()] = e->source();
			targets[e->index()] = e->target();
		}

		graph.reverseAllEdges();

		for(edge e : graph.edges) {
			AssertThat(e->source(), Equals(targets[e->index()]));
			AssertThat(e->target(), Equals(sources[e->index()]));
		}

		delete[] sources;
		delete[] targets;
	});

	for_each_graph_it("moves an adjacency entry", files, [](Graph &graph, const string file){
		adjEntry adj = chooseEdge(graph, 5)->adjSource();
		adjEntry adjSucc = adj->cyclicSucc();

		graph.moveAdj(adj, Direction::after, adjSucc);

		AssertThat(adjSucc->cyclicSucc(), Equals(adj));
		AssertThat(adj->cyclicSucc(), Is().Not().EqualTo(adjSucc));

		graph.moveAdj(adj, Direction::before, adjSucc);

		AssertThat(adj->cyclicSucc(), Equals(adjSucc));
		AssertThat(adjSucc->cyclicSucc(), Is().Not().EqualTo(adj));
	});

	for_each_graph_it("swaps the target of an edge", files, [](Graph &graph, const string file){
		edge e = graph.chooseEdge();
		node s = e->source();
		node t = e->target();

		node v = chooseNode(graph, t);

		graph.moveTarget(e, v);

		AssertThat(e->source(), Equals(s));
		AssertThat(e->target(), Equals(v));
	});

	for_each_graph_it("swaps the source of an edge", files, [](Graph &graph, const string file){
		edge e = graph.chooseEdge();
		node s = e->source();
		node t = e->target();

		node v = chooseNode(graph, s);

		graph.moveSource(e, v);

		AssertThat(e->source(), Equals(v));
		AssertThat(e->target(), Equals(t));
	});

	for_each_graph_it("splits an edge", files, [](Graph &graph, const string file){
		int n = graph.numberOfNodes();
		int m = graph.numberOfEdges();

		edge e = graph.chooseEdge();
		node v = e->target();

		edge f = graph.split(e);

		AssertThat(f->source(), Equals(e->target()));
		AssertThat(f->target(), Equals(v));
		AssertThat(f->source()->degree(), Equals(2));
		AssertThat(graph.numberOfNodes(), Equals(n + 1));
		AssertThat(graph.numberOfEdges(), Equals(m + 1));
	});

	for_each_graph_it("un-splits an edge by dummy-node", files, [](Graph &graph, const string file){
		int n = graph.numberOfNodes();
		int m = graph.numberOfEdges();

		edge e = graph.chooseEdge();
		node s = e->source();
		node t = e->target();

		graph.split(e);

		node v = e->target();

		graph.unsplit(v);

		AssertThat(graph.numberOfNodes(), Equals(n));
		AssertThat(graph.numberOfEdges(), Equals(m));
		AssertThat(e->source(), Equals(s));
		AssertThat(e->target(), Equals(t));
		AssertThat(graph.searchEdge(s, t), Equals(e));
	});

	for_each_graph_it("un-splits an edge by dummy-edge", files, [](Graph &graph, const string file){
		int n = graph.numberOfNodes();
		int m = graph.numberOfEdges();

		edge e = graph.chooseEdge();
		node s = e->source();
		node t = e->target();

		edge f = graph.split(e);
		graph.unsplit(e, f);

		AssertThat(graph.numberOfNodes(), Equals(n));
		AssertThat(graph.numberOfEdges(), Equals(m));
		AssertThat(e->source(), Equals(s));
		AssertThat(e->target(), Equals(t));
		AssertThat(graph.searchEdge(s, t), Equals(e));
	});

	for_each_graph_it("splits nodes", files, [](Graph &graph, const string file){
		node vLeft = chooseNode(graph, 6);

		int degree = vLeft->degree();
		List<adjEntry> entries;
		graph.adjEntries(vLeft, entries);
		adjEntry adjFirstRight = *entries.get(degree / 2);
		node vRight = graph.splitNode(vLeft->firstAdj(), adjFirstRight);
		int count = 0;

		for(adjEntry adj = vLeft->firstAdj()->succ(); adj != nullptr; adj = adj->succ()) {
			AssertThat(adj, Equals(*entries.get(count++)));
		}

		for(adjEntry adj = vRight->firstAdj()->succ(); adj != nullptr; adj = adj->succ()) {
			AssertThat(adj, Equals(*entries.get(count++)));
		}

		AssertThat(count, Equals(degree));
		AssertThat(vLeft->degree() + vRight->degree(), Equals(degree + 2));
	});

	for_each_graph_it("contracts an edge", files, [](Graph &graph, const string file){
		edge e = chooseEdge(graph, 5);
		node s = e->source();
		node t = e->target();

		// create the list of expected adjacency order
		List<node> nodes;
		List<edge> edges;
		graph.adjEdges(s, edges);

		for(edge f : edges) {
			nodes.pushBack(f->opposite(s));
		}

		// to prevent ambiguity, delete would-be multi-edges
		List<edge> deleteMe;
		graph.adjEdges(t, edges);
		ListIterator<node> it = nodes.search(t);

		while(edges.front() != e) {
			edges.moveToBack(edges.begin());
		}

		edges.del(edges.begin());

		for(edge f : edges) {
			if(nodes.search(f->opposite(t)).valid()) {
				deleteMe.pushBack(f);
			} else {
				nodes.insertBefore(f->opposite(t), it);
			}
		}

		nodes.del(it);

		for(edge f : deleteMe) {
			graph.delEdge(f);
		}

		node v = graph.contract(e);
		edge f = graph.searchEdge(v, nodes.front());

		AssertThat(v == t || v == s, IsTrue());
		AssertThat(v->degree(), Equals(nodes.size()));
		AssertThat(f, Is().Not().EqualTo((void*) nullptr));

		adjEntry adj = f->source() == v ? f->adjSource() : f->adjTarget();
		for(node w : nodes) {
			AssertThat(adj->twinNode(), Equals(w));
			adj = adj->cyclicSucc();
		}

	});

	for_each_graph_it("collapses half of all nodes", files, [](Graph &graph, const string file){
		int m = graph.numberOfEdges();

		List<node> nodes;
		int maxIndex = graph.maxNodeIndex();
		bool *adjacent = new bool[maxIndex + 1];

		for(int i = 0; i <= maxIndex; i++) {
			adjacent[i] = false;
		}

		for(node v : graph.nodes) {
			if(v->index() % 2) {
				nodes.pushBack(v);
			}
		}

		int minRemoved = 0;
		for(edge e : graph.edges) {
			int target = e->target()->index();
			int source = e->source()->index();

			if(source % 2 && target % 2 == 0) {
				adjacent[target] = true;
			}

			if(source % 2 == 0 && target % 2) {
				adjacent[source] = true;
			}

			minRemoved += source % 2 && target % 2;
		}

		node v = nodes.front();
		graph.collapse(nodes);

		AssertThat(nodes.empty(), IsTrue());
		AssertThat(graph.numberOfEdges(), IsLessThan(1 + m - minRemoved));

		for(adjEntry adj = v->firstAdj(); adj != nullptr; adj = adj->succ()) {
			adjacent[adj->twinNode()->index()] = false;
		}

		for(int i = 0; i <= maxIndex; i++) {
			AssertThat(adjacent[i], IsFalse());
		}

		delete[] adjacent;
	});

	for_each_graph_it("sorts adjacency lists", files, [](Graph &graph, const string file){
		node v = chooseNode(graph, 6);

		List<adjEntry> entries;
		graph.adjEntries(v, entries);

		entries.permute();

		graph.sort(v, entries);

		AssertThat(v->firstAdj(), Equals(entries.front()));
		AssertThat(v->lastAdj(), Equals(entries.back()));

		adjEntry adjBefore = nullptr;
		for(adjEntry adj : entries) {
			if(adjBefore != nullptr) {
				AssertThat(adjBefore->succ(), Equals(adj));
				AssertThat(adj->pred(), Equals(adjBefore));
			}

			adjBefore = adj;
		}
	});

	for_each_graph_it("reverses the order of all edges adjacent to a given node", files, [](Graph &graph, const string file){
		node v = chooseNode(graph, 6);
		List<edge> edges;
		graph.adjEdges(v, edges);

		graph.reverseAdjEdges(v);
		edges.reverse();

		adjEntry adj = v->firstAdj();
		for(edge e : edges) {
			AssertThat(adj, Is().Not().EqualTo((void*) nullptr));
			AssertThat(adj->theEdge(), Equals(e));

			adj = adj->succ();
		}
	});

	for_each_graph_it("swaps adjacency entries", files, [](Graph &graph, const string file){
		edge e = chooseEdge(graph, 5);
		adjEntry adj = e->adjSource()->cyclicSucc()->cyclicSucc();

		graph.swapAdjEdges(e->adjSource(), adj);

		AssertThat(adj->cyclicSucc()->cyclicSucc(), Equals(e->adjSource()));
		AssertThat(e->adjSource()->cyclicSucc()->cyclicSucc(), Is().Not().EqualTo(adj));
	});

	for_each_graph_it("does not return a negative genus", files, [](Graph &graph, const string file){
		AssertThat(graph.genus(), IsGreaterThan(-1));
	});

	for_each_graph_it("detects a combinatorial embedding", files, [](Graph &graph, const string file){
		AssertThat(graph.representsCombEmbedding(), Equals(graph.genus() == 0));
	});
});
});
