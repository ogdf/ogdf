/** \file
 * \brief Tests for Dijkstra's algorithm
 *
 * \author Jöran Schierbaum
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

#include <tuple>
#include <ogdf/basic/graph_generators/deterministic.h>
#include <ogdf/graphalg/Dijkstra.h>

#include <testing.h>
#include <graphs.h>

//! Storage container for an edge with a weight, indexed by integers
//! We use this instead of std::tuple to make use of aggregate initialization
template<typename T>
struct DijkstraTestEdge
{
	DijkstraTestEdge(int from, int to, T weight) : m_from(from), m_to(to), m_weight(weight) { }
	int m_from;
	int m_to;
	T m_weight;
};

//! Storage container for a node with a distance, indexed by integers
template<typename T>
struct DijkstraTestNode
{
	DijkstraTestNode(int idx, int predecessor, T distance) : m_idx(idx), m_pred(predecessor), m_distance(distance) { }
	int m_idx;
	int m_pred;
	T m_distance;
};

//! An instance of a graph Dijkstra can be called on, with expected resulting shortest path tree.
template<typename T>
struct DijkstraTestInstance
{

	//! Constructs a new test instance.
	//! @note In order to be valid, an instance needs to be supplied with a size using setSize(int).
	DijkstraTestInstance() : m_size(-1), m_directed(false), m_terminationTarget(-1), m_terminationDistance(-1) { }

	//! Configures the size of the instance.
	/**
	 * @param n number of nodes for the instance
	 * @return itself for method chaining
	 */
	DijkstraTestInstance& setSize(int n) {
		m_size = n;
		return *this;
	}

	//! Configure this instance to regard the edges as directed.
	/**
	 * @return itself for method chaining
	 */
	DijkstraTestInstance& setDirected() {
		m_directed = true;
		return *this;
	}

	//! Sets the edges with their weights
	/**
	 * @param edgeList a List of edges, each of which consists of its source node, target node and weight.
	 * @return itself for method chaining
	 */
	DijkstraTestInstance& edges(List<DijkstraTestEdge<T>> edgeList) {
		m_edges = edgeList;
		return *this;
	}

	//! Sets the expected node data after a normal Dijkstra run
	/**
	 * The predecessor relation and distance information resulting from a Dijkstra call on the configured
	 * instance are verified to equal this data in order for the test to succeed.
	 *
	 * @param expectedNodeData An array of node data, where each entry of this
	 *        array is a DijkstraTestNode. Each node holds a node identifier, its predecessor and the
	 *        distance to the instance's source node.
	 *        This array should have exactly one entry for every node of the graph.
	 * @return itself for method chaining
	 */
	DijkstraTestInstance& expectedPredecessorRelation(Array<DijkstraTestNode<T>> expectedNodeData) {
		m_nodes = expectedNodeData;
		return *this;
	}

	//! Sets the expected node data after a Dijkstra run bound by a target node
	/**
	 * @param node the id of a node that Dijkstra is bound by. The algorithm will terminate early as soon
	 *        as this node has been fully handled.
	 * @copydetails expectedPredecessorRelation(Array<DijkstraTestNode<T>>)
	 */
	DijkstraTestInstance& terminationTarget(int node, Array<DijkstraTestNode<T>> expectedNodeData) {
		m_terminationTarget = node;
		m_terminationTarget_nodes = expectedNodeData;
		return *this;
	}

	//! Sets the expected node data after a Dijkstra run bound by maximum distance
	/**
	 * @param distance the distance that Dijkstra is bound by. The algorithm will terminate early as soon
	 *        as this distance has been exhausted.
	 * @copydetails expectedPredecessorRelation(Array<DijkstraTestNode<T>>)
	 */
	DijkstraTestInstance& terminationDistance(T distance, Array<DijkstraTestNode<T>> expectedNodeData) {
		m_terminationDistance = distance;
		m_terminationDistance_nodes = expectedNodeData;
		return *this;
	}

	//! Call the standard Dijkstra algorithm and check the resulting predecessor relation
	//! @note needs to be called from within \c bandit::it()
	void testBasic() const {
		performTest({-1, -1, std::numeric_limits<T>::max()}, m_nodes);
	}

	//! Call the early terminated Dijkstra with parameters that are equal to the last expected node.
	//! This should produce the same result as the basic test.
	void testEarlyTerminationSafe() const {
		DijkstraTestNode<T> c{0, 0, 0};
		for (auto n : m_nodes) {
			if (n.m_distance > c.m_distance) {
				c = n;
			}
		}
		performTest(c, m_nodes);
	}

	//! Call Dijkstra with the configured target node and check if it adapts accordingly.
	void testEarlyTerminationTarget() const {
		// Only test if a scenario has been configured
		if (m_terminationTarget != -1) {
			performTest({m_terminationTarget, -1, std::numeric_limits<T>::max()}, m_terminationTarget_nodes);
		}
	}

	//! Call Dijkstra with the configured maximum distance and check if it adapts accordingly.
	void testEarlyTerminationDistance() const {
		// Only test if a scenario has been configured
		if (m_terminationDistance != -1) {
			performTest({-1, -1, m_terminationDistance}, m_terminationDistance_nodes);
		}
	}


protected:
	//! Call Dijkstra with the given targets and check the resulting predecessor relation against the parameter.
	/**
	 * @param target Configures early termination. If its id is not \c -1, use as target node id. Its distance
	 *        is used as the maximum distance for early termination.
	 * @param expectedNodeData expected results of Dijkstra
	 */
	void performTest(DijkstraTestNode<T> target, const Array<DijkstraTestNode<T>>& expectedNodeData) const {
		OGDF_ASSERT(m_size > 0);
		Graph G;
		Array<node> nodes;
		EdgeArray<T> weights(G);

		customGraph(G, m_size, {}, nodes);
		for (auto e : m_edges) {
			weights[G.newEdge(nodes[e.m_from], nodes[e.m_to])] = e.m_weight;
		}

		Dijkstra<T, PairingHeap> dij;
		NodeArray<edge> predecessor(G);
		NodeArray<T> distance(G);

		node targ = (target.m_idx == -1 ? nullptr : nodes[target.m_idx]);
		dij.call(G, weights, nodes[0], predecessor, distance, m_directed, false, targ, target.m_distance);

		for (int i = 0; i < m_size; ++i) {
			node current = nodes[expectedNodeData[i].m_idx];
			if (expectedNodeData[i].m_pred == -1) {
				AssertThat(predecessor[current], IsNull());
			}
			else {
				AssertThat(predecessor[current], !IsNull());
				AssertThat(predecessor[current]->opposite(current), Equals(nodes[expectedNodeData[i].m_pred]));
			}
			AssertThat(distance[current], Equals(expectedNodeData[i].m_distance));
		}
	}

	int m_size;
	bool m_directed;
	List<DijkstraTestEdge<T>> m_edges;
	Array<DijkstraTestNode<T>> m_nodes;
	int m_terminationTarget;
	Array<DijkstraTestNode<T>> m_terminationTarget_nodes;
	T m_terminationDistance;
	Array<DijkstraTestNode<T>> m_terminationDistance_nodes;
};


template<typename T>
void forAllInstances(std::function<void(const DijkstraTestInstance<T>&)> doTest) {

	auto testInstance = [&](const string& desc, std::function<void(DijkstraTestInstance<T>&)> populateInstance) {
		it("works on a " + desc, [&] {
			DijkstraTestInstance<T> instance;
			populateInstance(instance);
			doTest(instance);
		});
	};

	testInstance("three node graph", [](DijkstraTestInstance<T>& instance) {
		instance.setSize(3).edges({
			{0, 1,  5},
			{0, 2, 10},
			{1, 2,  4}
		}).expectedPredecessorRelation({
			{0, -1, 0},
			{1, 0, 5},
			{2, 1, 9}
		}).terminationTarget(1, {
			{0, -1, 0},
			{1, 0, 5},
			{2, 0, 10}
		});
	});

	testInstance("K5 with uniform weights", [](DijkstraTestInstance<T>& instance) {
		instance.setSize(5).edges({
			{0, 1, 1},
			{0, 2, 1},
			{0, 3, 1},
			{0, 4, 1},
			{1, 2, 1},
			{1, 3, 1},
			{1, 4, 1},
			{2, 3, 1},
			{2, 4, 1},
			{3, 4, 1}
		}).expectedPredecessorRelation({
			{0, -1, 0},
			{1, 0, 1},
			{2, 0, 1},
			{3, 0, 1},
			{4, 0, 1}
		}).terminationTarget(2, {
			{0, -1, 0},
			{1, 0, 1},
			{2, 0, 1},
			{3, 0, 1},
			{4, 0, 1}
		});
	});

	testInstance("disconnected graph", [](DijkstraTestInstance<T>& instance) {
		instance.setSize(5).edges({
			{0, 1, 5},
			{0, 2, 3},
			{1, 2, 1},
			{3, 4, 2}
		}).expectedPredecessorRelation({
			{0, -1, 0},
			{1, 2, 4},
			{2, 0, 3},
			{3, -1, std::numeric_limits<T>::max()},
			{4, -1, std::numeric_limits<T>::max()}
		}).terminationDistance(3, {
			{0, -1, 0},
			{1, -1, std::numeric_limits<T>::max()},
			{2, 0, 3},
			{3, -1, std::numeric_limits<T>::max()},
			{4, -1, std::numeric_limits<T>::max()}
		});
	});

	testInstance("graph without edges", [](DijkstraTestInstance<T>& instance) {
		instance.setSize(7).edges({
		}).expectedPredecessorRelation({
			{0, -1, 0},
			{1, -1, std::numeric_limits<T>::max()},
			{2, -1, std::numeric_limits<T>::max()},
			{3, -1, std::numeric_limits<T>::max()},
			{4, -1, std::numeric_limits<T>::max()},
			{5, -1, std::numeric_limits<T>::max()},
			{6, -1, std::numeric_limits<T>::max()}
		});
	});

	testInstance("graph with similar path lengths", [](DijkstraTestInstance<T>& instance) {
		instance.setSize(4).edges({
			{0, 1, 20},
			{0, 2, 20},
			{1, 3, 30},
			{2, 3, 29}
		}).expectedPredecessorRelation({
			{0, -1, 0},
			{1, 0, 20},
			{2, 0, 20},
			{3, 2, 49}
		});
	});

	testInstance("graph with similar path lengths on two symmetric sides", [](DijkstraTestInstance<T>& instance) {
		instance.setSize(8).edges({
			{0, 1, 1},
			{1, 2, 2},
			{1, 3, 5},
			{1, 4, 10},
			{2, 3, 2},
			{3, 4, 5},
			{0, 5, 1},
			{5, 6, 2},
			{5, 7, 5},
			{5, 4, 10},
			{6, 7, 2},
			{7, 4, 4}
		}).expectedPredecessorRelation({
			{0, -1, 0},
			{1, 0, 1},
			{2, 1, 3},
			{3, 2, 5},
			{4, 7, 9},
			{5, 0, 1},
			{6, 5, 3},
			{7, 6, 5}
		}).terminationTarget(2, {
			{0, -1, 0},
			{1, 0, 1},
			{2, 1, 3},
			{3, 1, 6},
			{4, 5, 11},
			{5, 0, 1},
			{6, 5, 3},
			{7, 5, 6}
		}).terminationDistance(5, {
			{0, -1, 0},
			{1, 0, 1},
			{2, 1, 3},
			{3, 2, 5},
			{4, -1, std::numeric_limits<T>::max()},
			{5, 0, 1},
			{6, 5, 3},
			{7, 6, 5}
		});
	});

	testInstance("directed graph", [](DijkstraTestInstance<T>& instance) {
		instance.setSize(5).setDirected().edges({
			{0, 1, 4},
			{0, 2, 2},
			{1, 2, 1},
			{2, 3, 5},
			{3, 4, 1},
			{4, 0, 1}
		}).expectedPredecessorRelation({
			{0, -1, 0},
			{1, 0, 4},
			{2, 0, 2},
			{3, 2, 7},
			{4, 3, 8}
		}).terminationTarget(2, {
			{0, -1, 0},
			{1, 0, 4},
			{2, 0, 2},
			{3, -1, std::numeric_limits<T>::max()},
			{4, -1, std::numeric_limits<T>::max()}
		}).terminationDistance(3, {
			{0, -1, 0},
			{1, -1, std::numeric_limits<T>::max()},
			{2, 0, 2},
			{3, -1, std::numeric_limits<T>::max()},
			{4, -1, std::numeric_limits<T>::max()}
		});
	});
}

template<typename T>
void compareDijkstraAlgorithms(const Graph &G, bool testMultipleSource)
{
	EdgeArray<T> weights(G, 1);

	Dijkstra<T, PairingHeap> dij;
	NodeArray<edge> predecessorBasic(G);
	NodeArray<edge> predecessorEarlyTerminated(G);
	NodeArray<edge> predecessorNotEarlyTerminated(G);
	NodeArray<T> distanceBasic(G);
	NodeArray<T> distanceEarlyTerminated(G);
	NodeArray<T> distanceNotEarlyTerminated(G);

	List<node> sources;
	if (testMultipleSource) {
		sources.pushBack(G.chooseNode());
		sources.pushBack(G.chooseNode([&](node s) { return !sources.search(s).valid(); }));
		sources.pushBack(G.chooseNode([&](node s) { return !sources.search(s).valid(); }));
		dij.callUnbound(G, weights, sources, predecessorBasic, distanceBasic, false);
	}
	else {
		dij.callUnbound(G, weights, G.firstNode(), predecessorBasic, distanceBasic, false);
	}

	T maxDistance = 0;
	node lastNode = nullptr;
	for (auto it = distanceBasic.begin(); it != distanceBasic.end(); ++it) {
		if (it.value() > maxDistance) {
			maxDistance = it.value();
			lastNode = it.key();
		}
	}
	if (testMultipleSource) {
		dij.callBound(G, weights, sources, predecessorEarlyTerminated, distanceEarlyTerminated, false, false, lastNode, maxDistance);
		dij.callBound(G, weights, sources, predecessorNotEarlyTerminated, distanceNotEarlyTerminated, false, false, nullptr, std::numeric_limits<T>::max());
	}
	else {
		dij.callBound(G, weights, G.firstNode(), predecessorEarlyTerminated, distanceEarlyTerminated, false, false, lastNode, maxDistance);
		dij.callBound(G, weights, G.firstNode(), predecessorNotEarlyTerminated, distanceNotEarlyTerminated, false, false, nullptr, std::numeric_limits<T>::max());
	}

	// We cannot compare predecessorBasic with the other predecessor relations, as the early
	// termination algorithm selects nodes in a different order.
	AssertThat(predecessorEarlyTerminated, EqualsContainer(predecessorNotEarlyTerminated));

	AssertThat(distanceBasic, EqualsContainer(distanceEarlyTerminated));
	AssertThat(distanceBasic, EqualsContainer(distanceNotEarlyTerminated));
}

template<typename T>
void performTestsSingleSource()
{
	describe("Finding a shortest path tree on simple instances", [] {
		forAllInstances<T>([](const DijkstraTestInstance<T>& instance) {
			instance.testBasic();
		});
	});
	describe("Finding the same shortest path tree with early termination parameters set to last node", [] {
		forAllInstances<T>([](const DijkstraTestInstance<T>& instance) {
			instance.testEarlyTerminationSafe();
		});
	});
	describe("Terminating early on a given target node", [] {
		forAllInstances<T>([](const DijkstraTestInstance<T>& instance) {
			instance.testEarlyTerminationTarget();
		});
	});
	describe("Terminating early on a given maximal distance", [] {
		forAllInstances<T>([](const DijkstraTestInstance<T>& instance) {
			instance.testEarlyTerminationDistance();
		});
	});
	describe("Finding the same shortest paths using different setups", [] {
		forEachGraphItWorks({}, [&] (const Graph& G) {
			if (G.numberOfNodes() == 0) return;
			compareDijkstraAlgorithms<T>(G, false);
		});
	});
}

template<typename T>
void performTestsMultipleSource()
{
	describe("Finding the same shortest paths using different setups", [] {
		forEachGraphItWorks({}, [&] (const Graph& G) {
			if (G.numberOfNodes() < 3) return;
			compareDijkstraAlgorithms<T>(G, true);
		});
	});
}
go_bandit([] {
	describe("Dijkstra<int> single source", [] {
		performTestsSingleSource<int>();
	});
	describe("Dijkstra<double> single source", [] {
		performTestsSingleSource<double>();
	});
	describe("Dijkstra<int> multiple sources", [] {
		performTestsMultipleSource<int>();
	});
	describe("Dijkstra<double> multiple sources", [] {
		performTestsMultipleSource<double>();
	});
});
