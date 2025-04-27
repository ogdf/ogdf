/** \file
 * \brief Tests for the GraphObserver class
 *
 * \author Max Ilsen
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

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/graph_generators.h>

#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <testing.h>

//! GraphObserver that counts how often its methods are called.
class CountingGraphObserver : public GraphObserver {
public:
	int m_nodesAdded = 0;
	int m_nodesDeleted = 0;
	int m_edgesAdded = 0;
	int m_edgesDeleted = 0;
	int m_cleared = 0;
	int m_reregistered = 0;
	const Graph* m_graph;

	explicit CountingGraphObserver(const Graph* G) : GraphObserver(), m_graph(nullptr) {
		reregister(G);
	}

	void nodeDeleted(node v) override {
		consistencyCheck();
		m_nodesDeleted++;
	};

	void nodeAdded(node v) override {
		consistencyCheck();
		m_nodesAdded++;
	};

	void edgeDeleted(edge e) override {
		consistencyCheck();
		m_edgesDeleted++;
	};

	void edgeAdded(edge e) override {
		consistencyCheck();
		adjEntry adjSrc = e->adjSource();
		adjEntry adjTgt = e->adjTarget();
		AssertThat(adjSrc, !IsNull());
		AssertThat(adjSrc->theNode(), !IsNull());
		AssertThat(adjSrc->twin(), !IsNull());
		AssertThat(adjTgt, !IsNull());
		AssertThat(adjTgt->theNode(), !IsNull());
		AssertThat(adjTgt->twin(), !IsNull());
		m_edgesAdded++;
	};

	void cleared() override {
		consistencyCheck();
		m_cleared++;
	};

	void registrationChanged(const Graph* old) override {
		AssertThat(m_graph, Equals(old));
		m_graph = getGraph();
#ifdef OGDF_DEBUG
		if (old) {
			old->consistencyCheck();
		}
#endif
		consistencyCheck();
		m_reregistered++;
	}

	void consistencyCheck() {
		AssertThat(m_graph, Equals(getGraph()));
#ifdef OGDF_DEBUG
		if (getGraph()) {
			getGraph()->consistencyCheck();
		}
#endif
	}
};

enum class GOEvent { NODE_DEL, NODE_ADD, EDGE_DEL, EDGE_ADD, CLEAR, REREGISTER };

class LoggingGraphObserver : public CountingGraphObserver {
public:
	explicit LoggingGraphObserver(const Graph* g) : CountingGraphObserver(g) { }

	std::vector<std::pair<GOEvent, int>> m_events;

	void nodeDeleted(node v) override {
		CountingGraphObserver::nodeDeleted(v);
		m_events.emplace_back(GOEvent::NODE_DEL, v->index());
	}

	void nodeAdded(node v) override {
		CountingGraphObserver::nodeAdded(v);
		m_events.emplace_back(GOEvent::NODE_ADD, v->index());
	}

	void edgeDeleted(edge e) override {
		CountingGraphObserver::edgeDeleted(e);
		m_events.emplace_back(GOEvent::EDGE_DEL, e->index());
	}

	void edgeAdded(edge e) override {
		CountingGraphObserver::edgeAdded(e);
		m_events.emplace_back(GOEvent::EDGE_ADD, e->index());
	}

	void cleared() override {
		CountingGraphObserver::cleared();
		m_events.emplace_back(GOEvent::CLEAR, 0);
	}

	void registrationChanged(const Graph* old) override {
		CountingGraphObserver::registrationChanged(old);
		m_events.emplace_back(GOEvent::REREGISTER, 0);
	}
};

go_bandit([]() {
	describe("GraphObserver", []() {
		it("observes Graph::newNode and Graph::delNode", []() {
			Graph G1;
			LoggingGraphObserver GO(&G1);
			AssertThat(GO.m_nodesAdded, Equals(0));
			AssertThat(GO.m_nodesDeleted, Equals(0));
			node v = G1.newNode();
			AssertThat(GO.m_nodesAdded, Equals(1));
			node w = G1.newNode();
			AssertThat(GO.m_nodesAdded, Equals(2));
			G1.delNode(v);
			AssertThat(GO.m_nodesDeleted, Equals(1));
			G1.delNode(w);
			AssertThat(GO.m_nodesDeleted, Equals(2));
			AssertThat(GO.m_events,
					Equals(std::vector<std::pair<GOEvent, int>> {
							{GOEvent::NODE_ADD, 0},
							{GOEvent::NODE_ADD, 1},
							{GOEvent::NODE_DEL, 0},
							{GOEvent::NODE_DEL, 1},
					}));
		});

		it("observes edge deletion caused by Graph::delNode", []() {
			Graph G1;
			LoggingGraphObserver GO(&G1);
			G1.newEdge(G1.newNode(), G1.newNode());
			AssertThat(GO.m_nodesAdded, Equals(2));
			AssertThat(GO.m_edgesAdded, Equals(1));
			G1.delNode(G1.nodes.head());
			AssertThat(GO.m_nodesDeleted, Equals(1));
			AssertThat(GO.m_edgesDeleted, Equals(1));
			AssertThat(GO.m_events,
					Equals(std::vector<std::pair<GOEvent, int>> {
							{GOEvent::NODE_ADD, 0},
							{GOEvent::NODE_ADD, 1},
							{GOEvent::EDGE_ADD, 0},
							{GOEvent::EDGE_DEL, 0},
							{GOEvent::NODE_DEL, 0},
					}));
		});

		it("observes Graph::newEdge and Graph::delEdge", []() {
			Graph G1;
			LoggingGraphObserver GO(&G1);
			AssertThat(GO.m_edgesAdded, Equals(0));
			AssertThat(GO.m_edgesDeleted, Equals(0));
			edge e = G1.newEdge(G1.newNode(), G1.newNode());
			AssertThat(GO.m_edgesAdded, Equals(1));
			edge f = G1.newEdge(G1.newNode(), G1.newNode());
			AssertThat(GO.m_edgesAdded, Equals(2));
			G1.delEdge(e);
			AssertThat(GO.m_edgesDeleted, Equals(1));
			G1.delEdge(f);
			AssertThat(GO.m_edgesDeleted, Equals(2));
			AssertThat(GO.m_events,
					Equals(std::vector<std::pair<GOEvent, int>> {
							{GOEvent::NODE_ADD, 0},
							{GOEvent::NODE_ADD, 1},
							{GOEvent::EDGE_ADD, 0},
							{GOEvent::NODE_ADD, 2},
							{GOEvent::NODE_ADD, 3},
							{GOEvent::EDGE_ADD, 1},
							{GOEvent::EDGE_DEL, 0},
							{GOEvent::EDGE_DEL, 1},
					}));
		});

		it("observes Graph::clear", []() {
			Graph G1;
			LoggingGraphObserver GO(&G1);
			AssertThat(GO.m_cleared, Equals(0));
			G1.clear();
			AssertThat(GO.m_cleared, Equals(1));
			randomGraph(G1, 5, 7); // calls clear
			AssertThat(GO.m_cleared, Equals(2));
			G1.clear();
			AssertThat(GO.m_cleared, Equals(3));
			AssertThat(GO.m_nodesDeleted, Equals(0));
			AssertThat(GO.m_edgesDeleted, Equals(0));

			AssertThat(GO.m_events.size(), Equals(size_t(1 + 1 + 5 + 7 + 1)));
			AssertThat(GO.m_events.at(0).first, Equals(GOEvent::CLEAR));
			AssertThat(GO.m_events.at(1).first, Equals(GOEvent::CLEAR));
			for (int i = 0; i < 5 + 7; i++) {
				if (GO.m_events.at(i + 2).first != GOEvent::NODE_ADD) {
					AssertThat(GO.m_events.at(i + 2).first, Equals(GOEvent::EDGE_ADD));
				}
			}
			AssertThat(GO.m_events.at(2 + 5 + 7).first, Equals(GOEvent::CLEAR));
		});

		it("observes Graph::operator=", []() {
			Graph G1;
			randomGraph(G1, 5, 7);
			Graph G2;
			LoggingGraphObserver GO(&G2);
			G2 = G1;
			AssertThat(GO.m_cleared, Equals(1));
			AssertThat(GO.m_nodesAdded, Equals(5));
			AssertThat(GO.m_edgesAdded, Equals(7));
			AssertThat(GO.m_edgesDeleted, Equals(0));
			AssertThat(GO.m_edgesDeleted, Equals(0));

			AssertThat(GO.m_events.size(), Equals(size_t(1 + 5 + 7)));
			AssertThat(GO.m_events.at(0).first, Equals(GOEvent::CLEAR));
			for (int i = 0; i < 5 + 7; i++) {
				if (GO.m_events.at(i + 1).first != GOEvent::NODE_ADD) {
					AssertThat(GO.m_events.at(i + 1).first, Equals(GOEvent::EDGE_ADD));
				}
			}
		});

		it("observes Graph::insert(Graph)", []() {
			Graph G1;
			randomGraph(G1, 5, 7);
			Graph G2;
			LoggingGraphObserver GO(&G2);
			G2.insert(G1);
			AssertThat(GO.m_cleared, Equals(0));
			AssertThat(GO.m_nodesAdded, Equals(5));
			AssertThat(GO.m_edgesAdded, Equals(7));
			AssertThat(GO.m_edgesDeleted, Equals(0));
			AssertThat(GO.m_edgesDeleted, Equals(0));

			AssertThat(GO.m_events.size(), Equals(size_t(5 + 7)));
			for (int i = 0; i < 5 + 7; i++) {
				if (GO.m_events.at(i).first != GOEvent::NODE_ADD) {
					AssertThat(GO.m_events.at(i).first, Equals(GOEvent::EDGE_ADD));
				}
			}
		});

		it("observes Graph deconstruction", []() {
			std::unique_ptr<Graph> G1 {new Graph()};
			LoggingGraphObserver GO(G1.get());
			AssertThat(GO.getGraph(), Equals(G1.get()));
			AssertThat(GO.m_reregistered, Equals(1));
			AssertThat(GO.m_graph, Equals(G1.get()));
			G1->newNode();
			AssertThat(GO.m_nodesAdded, Equals(1));
			G1.reset();
			AssertThat(GO.getGraph(), Equals(nullptr));
			AssertThat(GO.m_graph, Equals(nullptr));
			AssertThat(GO.m_reregistered, Equals(2));
			AssertThat(GO.m_cleared, Equals(0));
			AssertThat(GO.m_nodesAdded, Equals(1));
			AssertThat(GO.m_edgesAdded, Equals(0));
			AssertThat(GO.m_edgesDeleted, Equals(0));
			AssertThat(GO.m_edgesDeleted, Equals(0));

			{
				Graph G2;
				GO.reregister(&G2);
				AssertThat(GO.getGraph(), Equals(&G2));
				AssertThat(GO.m_graph, Equals(&G2));
				AssertThat(GO.m_reregistered, Equals(3));
				G2.newNode();
				AssertThat(GO.m_nodesAdded, Equals(2));
			}
			AssertThat(GO.getGraph(), Equals(nullptr));
			AssertThat(GO.m_graph, Equals(nullptr));
			AssertThat(GO.m_reregistered, Equals(4));
			AssertThat(GO.m_cleared, Equals(0));
			AssertThat(GO.m_nodesAdded, Equals(2));
			AssertThat(GO.m_edgesAdded, Equals(0));
			AssertThat(GO.m_edgesDeleted, Equals(0));
			AssertThat(GO.m_edgesDeleted, Equals(0));
		});

		it("observes Graph::split and Graph::unsplit", []() {
			Graph G1;
			LoggingGraphObserver GO(&G1);
			edge e = G1.newEdge(G1.newNode(), G1.newNode());
			AssertThat(e->index(), Equals(0));
			edge f = G1.split(e);
			AssertThat(f->index(), Equals(1));
			node c = e->commonNode(f);
			AssertThat(c->index(), Equals(2));
			G1.unsplit(c);
			AssertThat(GO.m_events,
					Equals(std::vector<std::pair<GOEvent, int>> {
							{GOEvent::NODE_ADD, 0},
							{GOEvent::NODE_ADD, 1},
							{GOEvent::EDGE_ADD, 0},
							{GOEvent::NODE_ADD, 2},
							{GOEvent::EDGE_ADD, 1},
							{GOEvent::EDGE_DEL, 1},
							{GOEvent::NODE_DEL, 2},
					}));
		});

		it("observes Graph::contract and Graph::collapse", []() {
			Graph G1;
			node c = G1.newNode();
			List<edge> E;
			for (int i = 0; i < 5; ++i) {
				node n = G1.newNode();
				E.pushBack(G1.newEdge(c, n));
			}
			node last = G1.nodes.tail();
			edge x = G1.newEdge(last, G1.newNode());
			LoggingGraphObserver GO1(&G1);
			for (edge e : E) {
				G1.contract(e);
			}
			AssertThat(x->source(), Equals(c));
			AssertThat(GO1.m_events,
					Equals(std::vector<std::pair<GOEvent, int>> {
							{GOEvent::EDGE_DEL, 0},
							{GOEvent::NODE_DEL, 1},
							{GOEvent::EDGE_DEL, 1},
							{GOEvent::NODE_DEL, 2},
							{GOEvent::EDGE_DEL, 2},
							{GOEvent::NODE_DEL, 3},
							{GOEvent::EDGE_DEL, 3},
							{GOEvent::NODE_DEL, 4},
							{GOEvent::EDGE_DEL, 4},
							{GOEvent::NODE_DEL, 5},
					}));

			Graph G2;
			List<node> N;
			c = G2.newNode();
			N.pushBack(c);
			for (int i = 0; i < 5; ++i) {
				node n = G2.newNode();
				G2.newEdge(c, n);
				N.pushBack(n);
			}
			last = G2.nodes.tail();
			x = G2.newEdge(last, G2.newNode());
			LoggingGraphObserver GO2(&G2);
			G2.collapse(N);
			AssertThat(x->source(), Equals(c));
			AssertThat(GO1.m_events, Equals(GO2.m_events));
		});

		it("cleanly unregisters from its Graph", []() {
			Graph G1;
			{
				LoggingGraphObserver GO(&G1);
				G1.newNode();
				AssertThat(GO.m_events,
						Equals(std::vector<std::pair<GOEvent, int>> {{GOEvent::NODE_ADD, 0}}));
				GO.reregister(nullptr);
				G1.newNode();
				AssertThat(GO.m_events,
						Equals(std::vector<std::pair<GOEvent, int>> {
								{GOEvent::NODE_ADD, 0}, {GOEvent::REREGISTER, 0}}));
			}
			{
				LoggingGraphObserver GO(&G1);
				G1.newNode();
				AssertThat(GO.m_events,
						Equals(std::vector<std::pair<GOEvent, int>> {{GOEvent::NODE_ADD, 2}}));
			}
			// a callback now would cause a memory error because the observers
			// are not in scope anmore
			G1.newNode();
		});
	});
});
