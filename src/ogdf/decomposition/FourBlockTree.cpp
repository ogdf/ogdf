/** \file
 * \brief Implementation of FourBlockTree using FourBlockTreeBuilder.
 *
 * \author Gregor Diatzko
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
#include <ogdf/basic/basic.h>
#include <ogdf/decomposition/FourBlockTree.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <exception>
#include <initializer_list>
#include <memory>
#include <utility>
#include <vector>

#ifdef OGDF_DEBUG
#	include <ogdf/basic/extended_graph_alg.h>
#	include <ogdf/basic/simple_graph_alg.h>
#endif // OGDF_DEBUG

using namespace ogdf;

template<typename _C, typename _L, typename _T = typename _C::value_type>
std::vector<_T> countingSort(const _C& input, _L key) {
	std::vector<size_t> counts;
	std::vector<std::vector<_T>> buckets;
	for (const auto& x : input) {
		const size_t k = key(x);
		if (k >= counts.size()) {
			counts.resize(k + 1, 0);
		}
		++counts[k];
	}
	for (size_t i = 1; i < counts.size(); ++i) {
		counts[i] += counts[i - 1];
	}
	std::vector<_T> res(input.size());
	for (auto it = input.rbegin(); it != input.rend(); ++it) {
		res[--counts[key(*it)]] = *it;
	}
	return res;
}

/**
 * A class that constructs the 4-block tree of a given graph.
 */
class FourBlockTreeBuilder {
	using triangle_t = std::array<edge, 3>;
	using returnSide_t = unsigned char;
	static constexpr returnSide_t LEFT = 0b10, RIGHT = 0b01; // when root is drawn at the top

	Graph& m_g;
	NodeArray<node>& m_originalNodes;

	/**
	 * A half-edge in g such that the external face of g lies to its left and
	 * m_externalFace == m_externalFace->theNode()->firstAdj().
	 */
	adjEntry m_externalFace;

	/**
	 * A node in the external face.
	 *
	 * This will be used as the root for DFS.
	 */
	node m_root;

	/**
	 * Index of every adjEntry in its adjencency list.
	 */
	AdjEntryArray<size_t> m_indices;

	/**
	 * Populate m_indices.
	 */
	void populateIndices();

	/**
	 * All separating triangles.
	 */
	std::vector<triangle_t> m_sepTriangles;

	/**
	 * Populate m_sepTriangles.
	 */
	void populateSepTriangles();

	/**
	 * The length of the root-v-path along tree edges for each node v.
	 */
	NodeArray<size_t> m_depth;

	/**
	 * The half-edge along which each node was found during the first DFS.
	 *
	 * m_parentEdge[v]->theNode() = v
	 */
	NodeArray<adjEntry> m_parentEdge; // m_parentEdge[v]->theNode() = v

	/**
	 * Whether each edge is a tree edge or a back edge.
	 */
	EdgeArray<bool> m_isTreeEdge;

	/**
	 * The lowpoint of each edge.
	 */
	EdgeArray<node> m_lowpoint;

	/**
	 * The return side of each edge.
	 */
	EdgeArray<returnSide_t> m_returnSide;

	/**
	 * The angular distance of each edge.
	 */
	EdgeArray<size_t> m_angularDistance;

	/**
	 * Run the first DFS.
	 *
	 * Populate m_depth, m_parentEdge, m_isTreeEdge, m_lowpoint, m_returnSide,
	 * m_angularDistance.
	 */
	void firstDfs();

	/**
	 * Order the elements of m_sepTriangles inside out.
	 */
	void orderTriangles();

	/**
	 * Split m_g along the elements of m_sepTriangles and build the resulting
	 * 4-block tree.
	 */
	std::unique_ptr<FourBlockTree> buildTree();

public:
	/**
	 * Prepares all necessary data structures.
	 *
	 * @param g The plane triangulated graph whose 4-block tree shall be constructed.
	 *          This graph will be used destructively.
	 *          Edge directions in g are not respected.
	 *          The order of edges at each node is used as the combinatorial
	 *          embedding.
	 * @param originalNodes The nodes in the original graph corresponding
	 *                      to those in g. This will be used to populate
	 *                      FourBlockTree::originalNodes.
	 * @param externalFace A half-edge in g such that the external face of g
	 *                     lies to its left.
	 *
	 * \pre externalFace == externalFace->theNode()->firstAdj()
	 */
	FourBlockTreeBuilder(Graph& g, NodeArray<node>& originalNodes, adjEntry externalFace);

	/**
	 * Run the algorithm.
	 */
	std::unique_ptr<FourBlockTree> call();
};

void FourBlockTreeBuilder::populateIndices() {
	for (const auto v : m_g.nodes) {
		size_t i = 0;
		for (const auto a : v->adjEntries) {
			m_indices[a] = i++;
		}
	}
}

void FourBlockTreeBuilder::populateSepTriangles() {
	/* find all triangles [Chiba, Nishizeki] */

	/* sort nodes by degree using CountingSort */
	const auto byDescDegree =
			countingSort(m_g.nodes, [&](node v) { return m_g.numberOfNodes() - v->degree(); });

	NodeArray<bool> marked(m_g, false);
	NodeArray<bool> deleted(m_g, false);
	NodeArray<adjEntry> a_vw(m_g, nullptr); // a_vw[w]->theEdge() = {v,w} for all w adj to v

	/* for all nodes ordered by non-ascending degree */
	for (const auto v : byDescDegree) {
		/* mark all neighbors of v */
		for (const auto a_vw_ : v->adjEntries) {
			const auto w = a_vw_->twinNode();
			marked[w] = true;
			a_vw[w] = a_vw_;
		}

		for (const auto a_vu : v->adjEntries) // for each marked node u
		{
			const auto u = a_vu->twinNode();
			if (deleted[u]) {
				continue;
			}
			for (const auto a_uw : u->adjEntries) // for each node w adj to u
			{
				const auto w = a_uw->twinNode();
				if (deleted[w]) {
					continue;
				}
				if (marked[w]) {
					/* {u,v,w} is a triangle */

					const auto ix_vu = m_indices[a_vu];
					const auto ix_vw = m_indices[a_vw[w]];
					if ((v->degree() + ix_vu - ix_vw + 1) % v->degree() > 2) { // <=> |ix_vu - ix_vw| > 1
						/* {u,v,w} is separating */
						m_sepTriangles.push_back(
								{a_vu->theEdge(), a_vw[w]->theEdge(), a_uw->theEdge()});
					}
				}
			}

			marked[u] = false;
		}

		deleted[v] = true;
	}
}

void FourBlockTreeBuilder::firstDfs() {
	NodeArray<bool> visited(m_g, false);
	EdgeArray<bool> traversed(m_g, false);
	std::vector<adjEntry> stack = {m_externalFace};
	visited[m_root] = traversed[m_externalFace->theEdge()] = true;
	m_parentEdge[m_root] = m_externalFace;

	/* if stack.back() = nullptr, we are done with the most recently visited node on level stack.size() */
	/* otherwise we are always in the subtree rooted at stack.back()->theNode(), traversing stack.back()->theEdge() */
	/* stack[m_depth[v]]->theNode() = v forall v on currently active path from root */

	while (!stack.empty()) {
		const auto adj = stack.back();
		if (adj) {
			const auto w = adj->twinNode();
			const auto e = adj->theEdge();
			if (!visited[w]) {
				/* e is tree edge */
				m_depth[w] = stack.size(); // = m_depth[v] + 1
				m_parentEdge[w] = adj->twin();
				visited[w] = true;
				traversed[e] = true;
				m_isTreeEdge[e] = true;
				if (!adj->isSource()) {
					m_g.reverseEdge(e);
				}
				stack.push_back(w->firstAdj());
			} else if (!traversed[e]) {
				/* e is back edge */
				traversed[e] = true;
				m_isTreeEdge[e] = false;
				m_lowpoint[e] = w;
				if (!adj->isSource()) {
					m_g.reverseEdge(e);
				}

				/* indices at w */
				const auto ix_e = m_indices[adj->twin()];
				const auto ix_p = m_indices[m_parentEdge[w]];
				const auto ix_d = m_indices[stack[m_depth[w]]];

				const size_t ang_e_d = (w->degree() + ix_d - ix_e) % w->degree();
				const size_t ang_e_p = (w->degree() + ix_p - ix_e) % w->degree();
				if (ang_e_d < ang_e_p) {
					/* in clockwise order around w: e, active tree edge = stack[m_depth[w]], m_parentEdge[w] */
					m_returnSide[e] = RIGHT;
					m_angularDistance[e] = ang_e_d;
				} else {
					m_returnSide[e] = LEFT;
					m_angularDistance[e] = w->degree() - ang_e_d;
				}

				/* go to next edge of v */
				stack.back() = adj->succ();
			} else {
				/* go to next edge of v */
				stack.back() = adj->succ();
			}
		} else {
			stack.pop_back();

			if (!stack.empty()) {
				/* w is done now, m_parentEdge[w] inherits properties of outgoing edges of w */
				const auto w = stack.back()->twinNode();
				const auto p = stack.back()->theEdge(); // = m_parentEdge[w]->theEdge()

				m_lowpoint[p] = w;
				m_returnSide[p] = 0;
				m_angularDistance[p] = 0;
				for (const auto a : w->adjEntries) {
					/* only consider outgoing edges of w */
					if (!a->isSource()) {
						continue;
					}
					const auto e = a->theEdge();

					if (m_depth[m_lowpoint[e]] < m_depth[m_lowpoint[p]]) {
						m_lowpoint[p] = m_lowpoint[e];
						m_returnSide[p] = m_returnSide[e];
						m_angularDistance[p] = m_angularDistance[e];
					} else if (m_depth[m_lowpoint[e]] == m_depth[m_lowpoint[p]]) {
						m_returnSide[p] |= m_returnSide[e];
						m_angularDistance[p] = std::max(m_angularDistance[e], m_angularDistance[e]);
					}
				}

				/* go to next edge of v */
				stack.back() = stack.back()->succ();
			}
		}
	}
}

void FourBlockTreeBuilder::orderTriangles() {
	/* sort edges e using Counting-/RadixSort on */
	/* 1. e->source() */
	/* 2. -m_depth[m_lowpoint[e]] */
	/* 3. m_angularDistance[e] */
	const auto byAngularDistance =
			countingSort(m_g.edges, [&](edge e) { return m_angularDistance[e]; });
	const auto byDepthLowpoint = countingSort(byAngularDistance,
			[&](edge e) { return m_g.numberOfNodes() - m_depth[m_lowpoint[e]]; });
	NodeArray<std::vector<edge>> edge_order(m_g);
	for (const auto e : byDepthLowpoint) {
		edge_order[e->source()].push_back(e);
	}

	std::vector<unsigned char> untraversedEdges(m_sepTriangles.size(),
			3); // number of edges for each triangle that have not yet been traversed
	EdgeArray<std::vector<size_t>> triangleIndices(
			m_g); // indices into untraversedEdges and m_sepTriangles
	for (size_t i = 0; i < m_sepTriangles.size(); ++i) {
		const auto& t = m_sepTriangles[i];
		for (const auto e : t) {
			triangleIndices[e].push_back(i);
		}
	}

	/* find order of triangles */
	/* indices into m_sepTriangles are also indices into these two arrays */
	std::vector<size_t> internalAngle(m_sepTriangles.size(), 0);
	std::vector<size_t> foundWithEdge(m_sepTriangles.size(), 0);

	/* second DFS */
	std::vector<std::pair<node, std::vector<edge>::const_iterator>> stack = {
			{m_root, edge_order[m_root].cbegin()}};
	size_t edge_id = 0;
	while (!stack.empty()) {
		const auto v = stack.back().first;
		auto& it = stack.back().second;
		if (it == edge_order[v].cend()) {
			stack.pop_back();
		} else {
			++edge_id;
			const auto e = *it++;
			const auto w = e->target();
			if (m_isTreeEdge[e]) {
				stack.emplace_back(w, edge_order[w].cbegin());
			}

			std::vector<size_t> completeTriangles; // elements are indices into m_sepTriangles
			for (const auto i : triangleIndices[e]) {
				if (--untraversedEdges[i] == 0) {
					completeTriangles.push_back(i);
					foundWithEdge[i] = edge_id;

					/* make e first edge in triangle */
					auto& t = m_sepTriangles[i];
					if (e == t[1]) {
						std::swap(t[0], t[1]);
					} else if (e == t[2]) {
						std::swap(t[0], t[2]);
					}

					/* "direct" t clockwise */
					if (t[1]->isIncident(w) ^ (m_returnSide[e] == LEFT)) {
						std::swap(t[1], t[2]);
					}
				}
			}
			if (completeTriangles.size() > 1) {
				/* e must be a back edge, and all triangles are on the same side of e */
				const bool trianglesRight = m_returnSide[e] == LEFT;

				const auto ix_e = m_indices[e->adjTarget()];

				for (const auto i : completeTriangles) {
					const auto& t = m_sepTriangles[i];
					const auto a = [&]() {
						/* find adjEntry a on t with a->theNode() == w and a->twinNode() != v */
						for (const auto side : t) {
							for (const auto adj : {side->adjSource(), side->adjTarget()}) {
								if (adj->theNode() == w && adj->twinNode() != v) {
									return adj;
								}
							}
						}
						/* unreachable */
						std::terminate();
					}();
					const auto ang_a_e = (w->degree() + ix_e - m_indices[a]) % w->degree();
					if (trianglesRight) {
						internalAngle[i] = ang_a_e;
					} else {
						internalAngle[i] = w->degree() - ang_a_e;
					}
				}
			}
		}
	}

	std::vector<size_t> tmp;
	tmp.reserve(m_sepTriangles.size());
	for (size_t i = 0; i < m_sepTriangles.size(); ++i) {
		tmp.push_back(i);
	}
	const auto byInternalAngle = countingSort(tmp, [&](size_t i) { return internalAngle[i]; });
	const auto byFoundWithEdge =
			countingSort(byInternalAngle, [&](size_t i) { return foundWithEdge[i]; });
	std::vector<triangle_t> res;
	for (const size_t i : byFoundWithEdge) {
		res.push_back(m_sepTriangles[i]);
	}

	m_sepTriangles = std::move(res);
}

std::unique_ptr<FourBlockTree> FourBlockTreeBuilder::buildTree() {
	std::vector<std::unique_ptr<FourBlockTree>> blocks;

	/* used to set FourBlockTree.parent */
	AdjEntryArray<FourBlockTree*> parentInv(m_g, nullptr);

	/* isInner[v] => v is part of an inner block that is about to be copied from m_g */
	/* initialized once to avoid linear cost for each block */
	NodeArray<bool> isInner(m_g, false);

	/* the ID the edge will have in b */
	/* initialized once to avoid linear cost for each block */
	EdgeArray<size_t> edgeIdsForB(m_g, 0);

	for (const auto& t : m_sepTriangles) {
		/* u,v,w in clockwise order around t */
		const auto [uv, vw, wu] = t;
		const auto v = uv->commonNode(vw);
		const auto w = vw->commonNode(wu);
		const auto u = wu->commonNode(uv);

		/* create new block b */
		blocks.push_back(std::make_unique<FourBlockTree>());
		auto& b = *blocks.back();

		/* split m_g at t */
		const auto u_ = m_g.newNode();
		const auto v_ = m_g.newNode();
		const auto w_ = m_g.newNode();
		const auto vw_ = m_g.newEdge(v_, w_);
		const auto wu_ = m_g.newEdge(w_, u_);
		const auto uv_ = m_g.newEdge(u_, v_);
		m_originalNodes[u_] = m_originalNodes[u];
		m_originalNodes[v_] = m_originalNodes[v];
		m_originalNodes[w_] = m_originalNodes[w];

		/* move edges incident to v */
		auto adj = vw->getAdj(v)->cyclicSucc();
		while (adj != uv->getAdj(v)) {
			const auto tmp = adj->cyclicSucc();
			if (adj->isSource()) {
				m_g.moveSource(adj->theEdge(), uv_->adjTarget(), Direction::before);
			} else {
				m_g.moveTarget(adj->theEdge(), uv_->adjTarget(), Direction::before);
			}
			adj = tmp;
		}

		/* move edges incident to w */
		adj = wu->getAdj(w)->cyclicSucc();
		while (adj != vw->getAdj(w)) {
			const auto tmp = adj->cyclicSucc();
			if (adj->isSource()) {
				m_g.moveSource(adj->theEdge(), vw_->adjTarget(), Direction::before);
			} else {
				m_g.moveTarget(adj->theEdge(), vw_->adjTarget(), Direction::before);
			}
			adj = tmp;
		}

		/* move edges incident to u */
		adj = uv->getAdj(u)->cyclicSucc();
		while (adj != wu->getAdj(u)) {
			const auto tmp = adj->cyclicSucc();
			if (adj->isSource()) {
				m_g.moveSource(adj->theEdge(), wu_->adjTarget(), Direction::before);
			} else {
				m_g.moveTarget(adj->theEdge(), wu_->adjTarget(), Direction::before);
			}
			adj = tmp;
		}

		/* set parent pointers to b */
		for (const auto a : {uv->getAdj(v), vw->getAdj(w), wu->getAdj(u)}) {
			auto& c = parentInv[a];
			if (c) {
				c->parent = &b;
			}
			c = &b;
		}

		/* list all nodes to be copied to b */
		std::vector<node> innerNodes = {u_, v_, w_}; // nodes in inner block
		isInner[u_] = isInner[v_] = isInner[w_] = true;
		size_t b_m = 0; // number of edges in inner block
		size_t nextEdgeIdForB = 0;
		for (size_t i = 0; i < innerNodes.size(); ++i) {
			for (const auto a : innerNodes[i]->adjEntries) {
				/* count edges */
				if (a->isSource()) {
					++b_m;
					edgeIdsForB[a] = nextEdgeIdForB++;
				}

				/* add neighbor x if not already added */
				const auto x = a->twinNode();
				if (!isInner[x]) {
					isInner[x] = true;
					innerNodes.push_back(x);
				}
			}
		}

		/* copy innerNodes with their edges to b */
		b.originalNodes.init(*b.g);
		b.parent = nullptr;
		b.parentFace = nullptr;
		const auto bDummy = b.g->newNode();
		std::vector<edge> bEdges;
		for (size_t i = 0; i < b_m; ++i) {
			bEdges.push_back(b.g->newEdge(bDummy, bDummy));
		}
		for (const auto x : innerNodes) {
			/* copy node */
			const auto b_v = b.g->newNode();
			b.originalNodes[b_v] = m_originalNodes[x];

			for (const auto a : x->adjEntries) {
				/* copy adjEntry */
				const auto b_e = bEdges[edgeIdsForB[a]];
				if (a->isSource()) {
					b.g->moveSource(b_e, b_v);
				} else {
					b.g->moveTarget(b_e, b_v);
				}

				/* set c.parent and c.parentFace for children c of b */
				const auto c = parentInv[a];
				if (c) {
					const auto b_a = a->isSource() ? b_e->adjSource() : b_e->adjTarget();
					c->parent = &b;
					c->parentFace = b_a;
				}
			}

			if (x == v_) {
				b.externalFace = b_v->lastAdj();
			}
		}
		b.g->delNode(bDummy);
	}

	/* root of 4-block tree */
	auto resP = std::make_unique<FourBlockTree>();
	{
		auto& res = *resP;
		/* list all nodes to be copied to res */
		std::vector<node> innerNodes = {m_root}; // nodes in outermost block
		isInner[m_root] = true;
		size_t b_m = 0; // number of edges in outermost block
		size_t nextEdgeIdForB = 0;
		for (size_t i = 0; i < innerNodes.size(); ++i) {
			for (const auto a : innerNodes[i]->adjEntries) {
				/* count edges */
				if (a->isSource()) {
					++b_m;
					edgeIdsForB[a] = nextEdgeIdForB++;
				}

				/* add neighbor w if not already added */
				const auto w = a->twinNode();
				if (!isInner[w]) {
					isInner[w] = true;
					innerNodes.push_back(w);
				}
			}
		}

		/* copy innerNodes with their edges to res */
		res.originalNodes.init(*res.g);
		res.parent = nullptr;
		res.parentFace = nullptr;
		const auto bDummy = res.g->newNode();
		std::vector<edge> bEdges;
		for (size_t i = 0; i < b_m; ++i) {
			bEdges.push_back(res.g->newEdge(bDummy, bDummy));
		}
		for (const auto v : innerNodes) {
			/* copy node */
			const auto b_v = res.g->newNode();
			res.originalNodes[b_v] = m_originalNodes[v];

			for (const auto a : v->adjEntries) {
				/* copy adjEntry */
				const auto b_e = bEdges[edgeIdsForB[a]];
				if (a->isSource()) {
					res.g->moveSource(b_e, b_v);
				} else {
					res.g->moveTarget(b_e, b_v);
				}

				/* set res.externalFace */
				if (a == m_externalFace) {
					res.externalFace = a->isSource() ? b_e->adjTarget() : b_e->adjSource();
				}

				/* set c.parent and c.parentFace for children c of res */
				const auto c = parentInv[a];
				if (c) {
					const auto b_a = a->isSource() ? b_e->adjSource() : b_e->adjTarget();
					c->parent = &res;
					c->parentFace = b_a;
				}
			}
		}
		res.g->delNode(bDummy);
	}

	/* move inner blocks to their parents */
	for (auto& b : blocks) {
		auto& p = *b->parent;
		p.children.push_back(std::move(b));
	}

	return resP;
}

FourBlockTreeBuilder::FourBlockTreeBuilder(Graph& g, NodeArray<node>& originalNodes,
		adjEntry externalFace)
	: m_g(g)
	, m_originalNodes(originalNodes)
	, m_externalFace(externalFace)
	, m_root(externalFace->theNode())
	, m_indices(g, 0)
	, m_depth(g, 0)
	, m_parentEdge(g, nullptr)
	, m_isTreeEdge(g, false)
	, m_lowpoint(g, nullptr)
	, m_returnSide(g, 0)
	, m_angularDistance(g, 0) { }

std::unique_ptr<FourBlockTree> FourBlockTreeBuilder::call() {
	populateIndices();
	populateSepTriangles();
	firstDfs();
	orderTriangles();
	return buildTree();
}

std::unique_ptr<FourBlockTree> FourBlockTree::construct(const Graph& g, adjEntry externalFace) {
	OGDF_ASSERT(externalFace != nullptr);
	OGDF_ASSERT(externalFace->graphOf() == &g);
	OGDF_ASSERT(g.numberOfNodes() * 3 == g.numberOfEdges() + 6 && "g must be triangulated");
	OGDF_ASSERT(isSimpleUndirected(g));
	OGDF_ASSERT(isPlanar(g));

	externalFace = externalFace->twin();

	/* set up working copy of g and populate originalNodes */
	Graph copy;
	NodeArray<node> originalNodes(copy, nullptr);
	EdgeArray<edge> edgeCopies(g, nullptr);
	const node dummySource = copy.newNode();
	const node dummyTarget = copy.newNode();
	for (const edge e : g.edges) {
		edgeCopies[e] = copy.newEdge(dummySource, dummyTarget);
	}
	for (const node v : g.nodes) {
		const node v_ = copy.newNode();
		originalNodes[v_] = v;
		for (const adjEntry a : v->adjEntries) {
			if (a->isSource()) {
				copy.moveSource(edgeCopies[a->theEdge()], v_);
			} else {
				copy.moveTarget(edgeCopies[a->theEdge()], v_);
			}
			if (a == externalFace) {
				externalFace = a->isSource() ? edgeCopies[a->theEdge()]->adjSource()
											 : edgeCopies[a->theEdge()]->adjTarget();
			}
		}
	}
	copy.delNode(dummySource);
	copy.delNode(dummyTarget);

	/* establish externalFace == externalFace->theNode()->firstAdj() */
	const node root = externalFace->theNode();
	while (externalFace != root->firstAdj()) {
		const adjEntry adj = root->firstAdj();
		if (adj->isSource()) {
			copy.moveSource(adj->theEdge(), root->lastAdj(), Direction::after);
		} else {
			copy.moveTarget(adj->theEdge(), root->lastAdj(), Direction::after);
		}
	}

	FourBlockTreeBuilder builder(copy, originalNodes, externalFace);
	return builder.call();
}
