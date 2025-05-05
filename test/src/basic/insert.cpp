/** \file
 * \brief Tests for ogdf::Graph::insert(...).
 *
 * \author Simon D. Fink
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

#include <ogdf/basic/CombinatorialEmbedding.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/GraphSets.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/graph_generators.h>
#include <ogdf/basic/simple_graph_alg.h>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <graphs.h>

#include <testing.h>

void randomTriangulate(Graph& G) {
	ogdf::makeConnected(G);
	ogdf::planarEmbed(G);

	CombinatorialEmbedding emb(G);
	std::vector<face> faces;
	for (face f : emb.faces) {
		if (f->size() > 3) {
			faces.push_back(f);
		}
	}
	while (!faces.empty()) {
		face f = faces.back();
		faces.pop_back();
		int a = randomNumber(0, f->size() - 1);
		int b = a;
		while (b == a) {
			b = randomNumber(0, f->size() - 1);
		}
		if (a > b) {
			std::swap(a, b);
		}
		b = (b - a) + 1;
		adjEntry adjA = f->firstAdj();
		while (a >= 1) {
			adjA = adjA->faceCycleSucc();
			a--;
		}
		adjEntry adjB = adjA;
		while (b >= 1) {
			adjB = adjB->faceCycleSucc();
			b--;
		}
		edge e = emb.splitFace(adjA, adjB);

		f = emb.leftFace(e->adjSource());
		if (f->size() > 3) {
			faces.push_back(f);
		}
		f = emb.leftFace(e->adjTarget());
		if (f->size() > 3) {
			faces.push_back(f);
		}
	}
	OGDF_ASSERT(emb.maximalFace()->size() <= 3);
}

class InsertLoggingGraph : public GraphCopySimple {
public:
	bool m_copyEmbedding = false;
	bool m_copyIDs = false;
	bool m_notifyObservers = false;
	bool m_edgeFilter = false;
	NodeArray<node>* m_nodeMap = nullptr;
	EdgeArray<edge>* m_edgeMap = nullptr;
	int* m_newNodes = nullptr;
	int* m_newEdges = nullptr;
	void* m_usrData = nullptr;
	std::set<node> m_nodeSet;
	std::set<edge> m_edgeSet;

protected:
	void* preInsert(bool copyEmbedding, bool copyIDs, bool notifyObservers, bool edgeFilter,
			NodeArray<node>& nodeMap, EdgeArray<edge>& edgeMap, int* newNodes,
			int* newEdges) override {
		m_copyEmbedding = copyEmbedding;
		m_copyIDs = copyIDs;
		m_notifyObservers = notifyObservers;
		m_edgeFilter = edgeFilter;
		m_nodeMap = &nodeMap;
		m_edgeMap = &edgeMap;
		m_newNodes = newNodes;
		m_newEdges = newEdges;
		m_usrData = GraphCopySimple::preInsert(copyEmbedding, copyIDs, notifyObservers, edgeFilter,
				nodeMap, edgeMap, newNodes, newEdges);
		return m_usrData;
	}

	void nodeInserted(void* userData, node original, node copy) override {
		AssertThat(userData, Equals(m_usrData));
		AssertThat((*m_nodeMap)[original], Equals(copy));
		GraphCopySimple::nodeInserted(userData, original, copy);
		m_nodeSet.insert(original);
	}

	void edgeInserted(void* userData, edge original, edge copy) override {
		AssertThat(userData, Equals(m_usrData));
		AssertThat((*m_edgeMap)[original], Equals(copy));
		GraphCopySimple::edgeInserted(userData, original, copy);
		m_edgeSet.insert(original);
	}

	void postInsert(void* userData, int newNodes, int newEdges) override {
		AssertThat(userData, Equals(m_usrData));
		GraphCopySimple::postInsert(userData, newNodes, newEdges);
		AssertThat(*m_newNodes, Equals(newNodes));
		AssertThat(*m_newEdges, Equals(newEdges));
		if (m_notifyObservers) {
			AssertThat(m_nodeSet.size(), Equals(size_t(newNodes)));
			AssertThat(m_edgeSet.size(), Equals(size_t(newEdges)));
		} else {
			AssertThat(m_nodeSet.size(), Equals(size_t(0)));
			AssertThat(m_edgeSet.size(), Equals(size_t(0)));
		}
	}
};

bool hasEvenEndpoints(edge e) {
	return e->source()->index() % 2 == 0 && e->target()->index() % 2 == 0;
}

std::unique_ptr<InsertLoggingGraph> ILG;
std::pair<int, int> res;
NodeArray<node> nodeMap;
EdgeArray<edge> edgeMap;

std::function<void()> reset(Graph& to) {
	return [&to] {
		ILG.reset(new InsertLoggingGraph());
		ILG->setOriginalGraph(to);
		nodeMap.init(to, nullptr);
		edgeMap.init(to, nullptr);
		res = {-1, -1};
	};
}

void checkDefaultFlags() {
	AssertThat(ILG->m_copyEmbedding, IsTrue());
	AssertThat(ILG->m_copyIDs, IsFalse());
	AssertThat(ILG->m_notifyObservers, IsTrue());
}

void checkRes(int n, int m) {
	AssertThat(res.first, Equals(n));
	AssertThat(res.second, Equals(m));
	if (ILG->m_notifyObservers) {
		AssertThat(ILG->m_nodeSet.size(), Equals(size_t(n)));
		AssertThat(ILG->m_edgeSet.size(), Equals(size_t(m)));
	}
}

template<typename NF, typename EF, typename NL, typename EL>
void checkSameGraph(const NF& nodeFilter, const EF& edgeFilter, const NL& nodes, const EL& edges,
		bool notifyObservers = true, bool copyEmbedding = true) {
	AssertThat(size_t(ILG->numberOfNodes()), Equals(size_t(nodes.size())));
	AssertThat(size_t(ILG->numberOfEdges()), Equals(size_t(edges.size())));
	if (notifyObservers) {
		AssertThat(ILG->m_nodeSet, Equals(std::set<node>(nodes.begin(), nodes.end())));
		AssertThat(ILG->m_edgeSet, Equals(std::set<edge>(edges.begin(), edges.end())));
	}

	for (node n : nodes) {
		if (!nodeFilter(n)) {
			continue;
		}

		node cnode = nodeMap[n];
		AssertThat(cnode, !IsNull());
		if (notifyObservers) {
			AssertThat(ILG->copy(n), Equals(cnode));
		}

		adjEntry cadj_it = cnode->firstAdj();
		std::vector<adjEntry> exp_emb, act_emb;
		for (adjEntry adj : n->adjEntries) {
			if (!edgeFilter(adj->theEdge())) {
				continue;
			}

			edge cedge = edgeMap[adj->theEdge()];
			AssertThat(edgeMap[adj->theEdge()], !IsNull());
			adjEntry cadj = adj->isSource() ? cedge->adjSource() : cedge->adjTarget();
			AssertThat(cadj->theNode(), Equals(nodeMap[n]));
			AssertThat(cadj->twinNode(), Equals(nodeMap[adj->twinNode()]));
			if (copyEmbedding) {
				exp_emb.emplace_back(cadj);
				act_emb.emplace_back(cadj_it);
				cadj_it = cadj_it->cyclicSucc();
			}
			if (notifyObservers) {
				AssertThat(cadj, Equals(ILG->copy(adj)));
				AssertThat(cedge, Equals(ILG->copy(adj->theEdge())));
				AssertThat(adj, Equals(ILG->original(cadj)));
				AssertThat(adj->theEdge(), Equals(ILG->original(cedge)));
				AssertThat(adj->theNode(), Equals(ILG->original(cadj->theNode())));
			}
		}
		AssertThat(act_emb, Equals(exp_emb));
	}
}

template<bool copyEmbedding = true, bool copyIDs = false, bool notifyObservers = true>
std::function<void()> testFiltered(Graph& orig, const std::function<bool(node)>& nodeFilter,
		const std::function<bool(edge)>& edgeFilter) {
	return [&]() {
		std::vector<node> nodes;
		std::vector<edge> edges;
		for (node n : orig.nodes) {
			if (nodeFilter(n)) {
				nodes.push_back(n);
			}
		}
		for (edge e : orig.edges) {
			if (edgeFilter(e)) {
				edges.push_back(e);
			}
		}

		std::default_random_engine rng(randomSeed());
		std::shuffle(nodes.begin(), nodes.end(), rng);
		std::shuffle(edges.begin(), edges.end(), rng);

		before_each(reset(orig));
		after_each([&] {
			if (ILG->m_edgeFilter) {
				AssertThat(ILG->m_copyEmbedding, Equals(true));
			} else {
				AssertThat(ILG->m_copyEmbedding, Equals(copyEmbedding));
			}
			AssertThat(ILG->m_copyIDs, Equals(copyIDs));
			AssertThat(ILG->m_notifyObservers, Equals(notifyObservers));
			checkRes(nodes.size(), edges.size());
			checkSameGraph(nodeFilter, edgeFilter, nodes, edges, notifyObservers, copyEmbedding);
		});

		if (copyEmbedding == true && copyIDs == false && notifyObservers == true) {
			// only available with default flags
			it("works via insert(nodeList, edgeList)", [&] {
				res = ILG->insert<std::vector<node>, std::vector<edge>>(nodes, edges, nodeMap,
						edgeMap);
				AssertThat(ILG->m_edgeFilter, IsFalse());
			});
		}
		it("works via insert(nodeBegin, nodesEnd, edgeFilter)", [&] {
			res = ILG->insert<std::vector<node>::iterator, std::function<bool(edge)>, copyIDs,
					notifyObservers>(nodes.begin(), nodes.end(), edgeFilter, nodeMap, edgeMap);
			AssertThat(ILG->m_edgeFilter, IsTrue());
		});
		it("works via insert(nodeBegin, nodesEnd, edgesBegin, edgesEnd)", [&] {
			res = ILG->insert<std::vector<node>::iterator, std::vector<edge>::iterator,
					copyEmbedding, copyIDs, notifyObservers>(nodes.begin(), nodes.end(),
					edges.begin(), edges.end(), nodeMap, edgeMap);
			AssertThat(ILG->m_edgeFilter, IsFalse());
		});
	};
}

go_bandit([]() {
	Graph orig;
	{
		Graph H;
		randomPlanarCNBGraph(H, 10, 20, 10);
		orig.insert(H);
		H.clear();
		randomPlanarCNBGraph(H, 10, 20, 5);
		orig.insert(H);
		H.clear();
		randomPlanarBiconnectedGraph(H, 50, 100, true);
		orig.insert(H);
		edge e = orig.chooseEdge();
		orig.newEdge(e->source(), e->target());
		AssertThat(orig.numberOfNodes(), IsGreaterThan(100));
		AssertThat(orig.numberOfEdges(), IsGreaterThan(200));
		AssertThat(connectedComponents(orig), Equals(3));
		EdgeArray<int> out(orig);
		AssertThat(biconnectedComponents(orig, out), Equals(16));
		AssertThat(planarEmbed(orig), IsTrue());
		orig.reverseEdge(orig.newEdge(e->adjSource(), e->adjTarget()->cyclicSucc()));
		AssertThat(isPlanar(orig), IsTrue());
	}

	describe("inserting the whole graph", [&]() {
		before_each(reset(orig));
		after_each([&] {
			checkDefaultFlags();
			checkRes(orig.numberOfNodes(), orig.numberOfEdges());
			checkSameGraph(filter_any_node, filter_any_edge, orig.nodes, orig.edges);
		});

		it("works via insert(G)", [&] {
			res = ILG->insert(orig);
			for (node n : orig.nodes) {
				nodeMap[n] = ILG->copy(n);
			}
			for (edge e : orig.edges) {
				edgeMap[e] = ILG->copy(e);
			}
			AssertThat(ILG->m_edgeFilter, IsTrue());
		});
		it("works via insert(G, nodeMap, edgeMap)", [&] {
			res = ILG->insert(orig, nodeMap, edgeMap);
			AssertThat(ILG->m_edgeFilter, IsTrue());
		});
		it("works via insert(CCsInfo, ccNr)", [&] {
			Graph::CCsInfo ccs(orig);
			res = {0, 0};
			std::set<node> nodes;
			std::set<edge> edges;
			for (int i = 0; i < ccs.numberOfCCs(); ++i) {
				auto r = ILG->insert(ccs, i, nodeMap, edgeMap);
				res.first += r.first;
				res.second += r.second;

				AssertThat(r.first, Equals(ccs.numberOfNodes(i)));
				AssertThat(r.second, Equals(ccs.numberOfEdges(i)));
				AssertThat(ILG->m_nodeSet.size(), Equals(size_t(ccs.numberOfNodes(i))));
				AssertThat(ILG->m_edgeSet.size(), Equals(size_t(ccs.numberOfEdges(i))));

				nodes.insert(ILG->m_nodeSet.begin(), ILG->m_nodeSet.end());
				edges.insert(ILG->m_edgeSet.begin(), ILG->m_edgeSet.end());
				ILG->m_nodeSet.clear();
				ILG->m_edgeSet.clear();
			}
			std::swap(nodes, ILG->m_nodeSet);
			std::swap(edges, ILG->m_edgeSet);
			AssertThat(ILG->m_edgeFilter, IsTrue());
		});
		it("works via insert(G.nodes, G.edges)", [&] {
			res = ILG->insert(orig.nodes, orig.edges, nodeMap, edgeMap);
			AssertThat(ILG->m_edgeFilter, IsFalse());
		});
		it("works via insert(G.nodes.begin, G.nodes.end, filter_any_edge)", [&] {
			res = ILG->insert(orig.nodes.begin(), orig.nodes.end(), filter_any_edge, nodeMap,
					edgeMap);
			AssertThat(ILG->m_edgeFilter, IsTrue());
		});
		it("works via insert(G.nodes.begin, G.nodes.end, G.edges.begin, G.edges.end)", [&] {
			res = ILG->insert(orig.nodes.begin(), orig.nodes.end(), orig.edges.begin(),
					orig.edges.end(), nodeMap, edgeMap);
			AssertThat(ILG->m_edgeFilter, IsFalse());
		});
	});

	describe("inserting the empty graph",
			testFiltered(
					orig, [](node n) { return false; }, [](edge e) { return false; }));
	describe("inserting all nodes and edges", testFiltered(orig, filter_any_node, filter_any_edge));
	describe("inserting some nodes and no edges",
			testFiltered(
					orig, [](node n) { return n->index() % 2 == 0; }, [](edge e) { return false; }));
	describe("inserting some nodes and all their edges",
			testFiltered(
					orig, [](node n) { return n->index() % 2 == 0; }, hasEvenEndpoints));
	describe("inserting some nodes and some edges (even)",
			testFiltered(
					orig, [](node n) { return n->index() % 2 == 0; },
					[](edge e) { return e->index() % 2 == 0 && hasEvenEndpoints(e); }));

	describe("inserting some nodes and some edges (triangulation)", [&]() {
		GraphCopySimple origCon(orig);
		NodeSet origN(origCon);
		EdgeSet origE(origCon);
		for (node n : origCon.nodes) {
			origN.insert(n);
		}
		for (edge e : origCon.edges) {
			origE.insert(e);
		}
		{
			Graph H;
			ogdf::wheelGraph(H, 10);
			origCon.insert(H);
		}
		randomTriangulate(origCon);

		describe("with template flags", [&] {
			describe("<copyEmbedding=true, copyIDs=true, notifyObservers=true>",
					testFiltered<true, true, true>(origCon, origN, origE));
			describe("<copyEmbedding=true, copyIDs=true, notifyObservers=false>",
					testFiltered<true, true, false>(origCon, origN, origE));
			describe("<copyEmbedding=true, copyIDs=false, notifyObservers=true>",
					testFiltered<true, false, true>(origCon, origN, origE)); // default
			describe("<copyEmbedding=true, copyIDs=false, notifyObservers=false>",
					testFiltered<true, false, false>(origCon, origN, origE));
			describe("<copyEmbedding=false, copyIDs=true, notifyObservers=true>",
					testFiltered<false, true, true>(origCon, origN, origE));
			describe("<copyEmbedding=false, copyIDs=true, notifyObservers=false>",
					testFiltered<false, true, false>(origCon, origN, origE));
			describe("<copyEmbedding=false, copyIDs=false, notifyObservers=true>",
					testFiltered<false, false, true>(origCon, origN, origE));
			describe("<copyEmbedding=false, copyIDs=false, notifyObservers=false>",
					testFiltered<false, false, false>(origCon, origN, origE));
		});

		describe("with an EdgeSet as filter", testFiltered(origCon, origN, origE));

		describe("with other insert variants", [&]() {
			before_each(reset(origCon));
			after_each([&] {
				checkDefaultFlags();
				checkRes(origN.size(), origE.size());
				checkSameGraph(origN, origE, origN, origE);
			});

			it("works via insert(nodeList, EdgeSet)", [&] {
				res = ILG->insert(origN.elements(), origE, nodeMap, edgeMap);
				AssertThat(ILG->m_edgeFilter, IsTrue());
			});
			it("works via insert(nodeList, EdgeSet.elements List)", [&] {
				res = ILG->insert(origN.elements(), origE.elements(), nodeMap, edgeMap);
				AssertThat(ILG->m_edgeFilter, IsFalse());
			});
			it("works via insert(nodeBegin, nodesEnd, EdgeSet.begin + end)", [&] {
				res = ILG->insert(origN.begin(), origN.end(), origE.begin(), origE.end(), nodeMap,
						edgeMap);
				AssertThat(ILG->m_edgeFilter, IsFalse());
			});
		});
	});

	describe("inserting preserves direction", [] {
		it("works forward", [] {
			Graph G;
			customGraph(G, 2, {{0, 1}});
			reset(G)();
			res = ILG->insert(G, nodeMap, edgeMap);

			checkDefaultFlags();
			checkRes(2, 1);
			checkSameGraph(filter_any_node, filter_any_edge, G.nodes, G.edges);
		});
		it("works backward", [] {
			Graph G;
			customGraph(G, 2, {{1, 0}});
			reset(G)();
			res = ILG->insert(G, nodeMap, edgeMap);

			checkDefaultFlags();
			checkRes(2, 1);
			checkSameGraph(filter_any_node, filter_any_edge, G.nodes, G.edges);
		});
	});
	describe("inserting works for all kinds of graphs", [] {
		forEachGraphItWorks({}, [](Graph& G) {
			reset(G)();
			res = ILG->insert(G, nodeMap, edgeMap);

			checkDefaultFlags();
			checkRes(G.numberOfNodes(), G.numberOfEdges());
			checkSameGraph(filter_any_node, filter_any_edge, G.nodes, G.edges);
		});
	});
});
