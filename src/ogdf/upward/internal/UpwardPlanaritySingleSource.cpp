/** \file
 * \brief Implements class UpwardPlanaritySingleSource, which represents
 * the upward planarity testing and embedding algorithm for single-source digraphs.
 *
 * Reference: "Optimal upward planarity testing of single-source
 *  digraphs" P. Bertolazzi, G. Di Battista, C. Mannino, and
 *  R. Tamassia, SIAM J.Comput., 27(1) Feb. 1998, pp. 132-169
 *
 * \author Carsten Gutwenger
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
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/decomposition/PertinentGraph.h>
#include <ogdf/decomposition/SPQRTree.h>
#include <ogdf/decomposition/Skeleton.h>
#include <ogdf/decomposition/StaticPlanarSPQRTree.h>
#include <ogdf/decomposition/StaticSkeleton.h>
#include <ogdf/upward/ExpansionGraph.h>
#include <ogdf/upward/FaceSinkGraph.h>
#include <ogdf/upward/UpwardPlanarity.h>
#include <ogdf/upward/internal/UpwardPlanaritySingleSource.h>

#include <ostream>

namespace ogdf {


class UpwardPlanaritySingleSource::SkeletonInfo {
public:
	SkeletonInfo() { }

	explicit SkeletonInfo(const Skeleton& S)
		: m_degInfo(S.getGraph()), m_containsSource(S.getGraph(), false) { }

	void init(const Skeleton& S) {
		m_degInfo.init(S.getGraph());
		m_containsSource.init(S.getGraph(), false);
	}

	EdgeArray<DegreeInfo> m_degInfo;
	EdgeArray<bool> m_containsSource;
	ConstCombinatorialEmbedding m_E;
	FaceSinkGraph m_F;
	SList<face> m_externalFaces;
};

class UpwardPlanaritySingleSource::ConstraintRooting {
public:
	explicit ConstraintRooting(const SPQRTree& T);

	// constrains a Q-node vQ associated with real edge e to be directed
	// leaving vQ
	void constrainRealEdge(edge e);

	// constrains a tree edge e in original SPQR-tree T to be directed
	// leaving src; returns false iff e was already constrained to be directed
	// in opposite direction
	bool constrainTreeEdge(edge e, node src);

	// if a roting satisfying all constraints exits, return a real edge at
	// which the SPQR-tree can be rooted, otherwise 0 is returned
	edge findRooting();

	void outputConstraints(std::ostream& os);

	ConstraintRooting& operator=(const ConstraintRooting&) = delete;

private:
	bool checkEdge(edge e, node parent, EdgeArray<bool>& checked);

	Graph m_tree; // tree with Q-nodes
	const SPQRTree& m_T; // original SPQR-tree

	// edge in m_tree corresponding to real edge
	EdgeArray<edge> m_realToConstraint;
	// node in m_tree corresponding to node in original SPQR-tree T
	NodeArray<node> m_internalNode;
	// edge in m_tree corresponding to edge in original SPQR-tree T
	EdgeArray<edge> m_treeToConstraint;
	// true <=> edge in m_tree has been constrained
	EdgeArray<bool> m_isConstrained;
};

// constructor
// builds copy of SPQR-tree with Q-nodes
UpwardPlanaritySingleSource::ConstraintRooting::ConstraintRooting(const SPQRTree& T)
	: m_T(T), m_isConstrained(m_tree, false) {
	const Graph& GT = T.tree();

	m_internalNode.init(GT);

	for (node v : GT.nodes) {
		m_internalNode[v] = m_tree.newNode();
	}

	m_treeToConstraint.init(GT);

	for (edge e : GT.edges) {
		m_treeToConstraint[e] =
				m_tree.newEdge(m_internalNode[e->source()], m_internalNode[e->target()]);
	}

	// create Q-nodes and adjacent edges
	const Graph& G = T.originalGraph();
	m_realToConstraint.init(G);

	for (edge e : G.edges) {
		node qNode = m_tree.newNode();
		node internal = m_internalNode[T.skeletonOfReal(e).treeNode()];
		// An edge adjacent with a Q-node qNode is always directed leaving
		// qNode because this is the only constraint used for those edges.
		// If we have to constrain such an edge, we simply mark it constrained
		// because it is already directed correctly.
		m_realToConstraint[e] = m_tree.newEdge(qNode, internal);
	}
}

// output for debugging
void UpwardPlanaritySingleSource::ConstraintRooting::outputConstraints(std::ostream& os) {
	const Graph& G = m_T.originalGraph();
	const Graph& GT = m_T.tree();

	os << "constrained edges in tree:\n";
	os << "real edges:";

	for (edge e : G.edges) {
		if (m_isConstrained[m_realToConstraint[e]]) {
			os << " " << e;
		}
	}

	os << "\ntree edges:";
	for (edge e : GT.edges) {
		if (m_isConstrained[m_treeToConstraint[e]]) {
			if (m_internalNode[e->source()] == m_treeToConstraint[e]->source()) {
				os << " " << e->source() << "->" << e->target();
			} else {
				os << " " << e->target() << "->" << e->source();
			}
		}
	}
	os << std::endl;
}

// constrains a real edge to be directed leaving its Q-node
// (such a Q-node cannot be selected as root node)
void UpwardPlanaritySingleSource::ConstraintRooting::constrainRealEdge(edge e) {
	m_isConstrained[m_realToConstraint[e]] = true;
}

// constrains a tree edge to be directed leaving src
// if it was already constrained to be directed in opposite direction
// false is returned and no rooting of the tree is possible
bool UpwardPlanaritySingleSource::ConstraintRooting::constrainTreeEdge(edge e, node src) {
	OGDF_ASSERT(src == e->source() || src == e->target());

	edge eTree = m_treeToConstraint[e];
	node srcTree = m_internalNode[src];

	// tree edge in wrong direction ?
	if (srcTree != eTree->source()) {
		// if it was already constriant we have a contradiction
		if (m_isConstrained[eTree]) {
			return false;
		}

		m_tree.reverseEdge(eTree);
	}

	m_isConstrained[eTree] = true;
	return true;
}

// find a rooting of the tree satisfying all constraints
// if such a rooting exists, a root edge is returned (an edge in the
// original graph), otherwise 0 is returned
edge UpwardPlanaritySingleSource::ConstraintRooting::findRooting() {
	EdgeArray<bool> checked(m_tree, false);

	// we check each constrained edge
	// such a constraint constrains the complete subtree rooted at
	// e->source(), hence procedure checkEdge() recursively checks
	// all (not-checked) edges in this subtree
	for (edge e : m_tree.edges) {
		if (m_isConstrained[e] && !checkEdge(e, e->target(), checked)) {
			return nullptr;
		}
	}

	// if all constraints are checked, we can find a rooting iff there is
	// a Q-node whose (only) adjacent edge was not constrained
	const Graph& G = m_T.originalGraph();
	for (edge e : G.edges) {
		edge eQ = m_realToConstraint[e];

		if (!checked[eQ]) {
			return e;
		}
	}

	return nullptr;
}

bool UpwardPlanaritySingleSource::ConstraintRooting::checkEdge(edge e, node parent,
		EdgeArray<bool>& checked) {
	if (checked[e]) {
		return e->target() == parent;
	}

	if (e->target() != parent) {
		if (m_isConstrained[e]) {
			return false;
		} else {
			m_tree.reverseEdge(e);
		}
	}

	checked[e] = true;
	node child = e->source();

	for (adjEntry adj : child->adjEntries) {
		edge eAdj = adj->theEdge();
		if (eAdj != e && !checkEdge(eAdj, child, checked)) {
			return false;
		}
	}

	return true;
}

// tests if a single-source digraph G is upward-planar
// computes sorted adjacency lists of an upward-planar embedding if
// embed is true
bool UpwardPlanaritySingleSource::testAndFindEmbedding(const Graph& G, bool embed,
		NodeArray<SListPure<adjEntry>>& adjacentEdges) {
	// trivial cases
	if (G.empty()) {
		return true;
	}

	if (!isAcyclic(G)) {
		return false;
	}

	// build expansion graph of G
	ExpansionGraph exp(G);

	// determine single source; if not present (or several sources) test fails
	node sG;
	if (!hasSingleSource(G, sG)) {
		return false;
	}

	// If embed is true we compute also the list of adjacency entries for each
	// node of G in an upward planar embedding.
	// Function testBiconnectedComponent() iterates over all biconnected
	// components of exp starting with the components adjacent to the single
	// source in exp.
	return testBiconnectedComponent(exp, sG, -1, embed, adjacentEdges);
}

// embeds single-source digraph G upward-planar
// also computes a planar st-augmentation of G and returns the list of
// augmented nodes and edges if augment is true
void UpwardPlanaritySingleSource::embedAndAugment(Graph& G,
		NodeArray<SListPure<adjEntry>>& adjacentEdges, bool augment, node& superSink,
		SList<edge>& augmentedEdges) {
	for (node vG : G.nodes) {
		G.sort(vG, adjacentEdges[vG]);
	}

	// the following tests check if the assigned embedding is upward planar
#ifdef OGDF_DEBUG
	G.consistencyCheck();
#endif
	OGDF_ASSERT(G.representsCombEmbedding());

#ifdef OGDF_DEBUG
	if (!augment) {
		CombinatorialEmbedding E(G);
		SList<face> externalFaces;
		OGDF_ASSERT(UpwardPlanarity::isUpwardPlanar_singleSource_embedded(E, externalFaces));
	}
#endif

	if (augment) {
#ifdef OGDF_DEBUG
		bool isUpwardPlanar =
#endif
				UpwardPlanarity::upwardPlanarAugment_singleSource_embedded(G, superSink,
						augmentedEdges);
		OGDF_ASSERT(isUpwardPlanar);
	}
}

// performs the actual test (and computation of sorted adjacency lists) for
// each biconnected component
// Within this procedure all biconnected components in which sG is (the
// original vertex) of the single-source are handled, for all other vertices
// in the biconnected component the procedure is called recursively.
// This kind of traversal is crucial for the embedding algorithm
bool UpwardPlanaritySingleSource::testBiconnectedComponent(ExpansionGraph& exp, node sG,
		int parentBlock, bool embed, NodeArray<SListPure<adjEntry>>& adjacentEdges) {
	SListConstIterator<int> itBlock;
	for (itBlock = exp.adjacentComponents(sG).begin(); itBlock.valid(); ++itBlock) {
		int i = *itBlock;
		if (i == parentBlock) {
			continue;
		}

		exp.init(i);

		// cannot construct SPQR-tree of graph with a single edge so treat
		// it as special case here
		if (exp.numberOfNodes() == 2) {
			edge eST = exp.original(exp.firstEdge());
			if (embed) {
				node src = eST->source();
				node tgt = eST->target();

				for (edge e : exp.edges) {
					edge eG = exp.original(e);
					adjacentEdges[src].pushBack(eG->adjSource());
					adjacentEdges[tgt].pushFront(eG->adjTarget());
				}
			}

			bool testOk = testBiconnectedComponent(exp, eST->target(), i, embed, adjacentEdges);
			if (!testOk) {
				return false;
			}

			continue;
		}

		// test whether expansion graph is planar
		if (!isPlanar(exp)) {
			return false;
		}

		// construct SPQR-tree T of exp with embedded skeleton graphs
		StaticPlanarSPQRTree T(exp);
		const Graph& tree = T.tree();

		// skeleton info maintains precomputed information, i.e., degrees of
		// nodes in expansion graph of virtual edges, if expansion graph
		// constains the single source
		NodeArray<SkeletonInfo> skInfo(tree);
		for (node vT : tree.nodes) {
			skInfo[vT].init(T.skeleton(vT));
		}

		// the single source in exp
		node s = exp.copy(sG);
		OGDF_ASSERT(exp.original(s) == sG);

		// precompute information maintained in skInfo
		// for each virtual edge e of a skeleton, determine its in- and
		// outdegree in the pertinent digraph of e; also determine if the
		// pertinent digraph of e contains the source
		computeDegreesInPertinent(T, s, skInfo, T.rootNode());

		OGDF_HEAVY_ASSERT(checkDegrees(T, s, skInfo));

		// For each R-node vT:
		//   compute its sT-skeleton, test whether the sT-skeleton is upward-
		//     planar
		//   mark the virtual edges of skeleton(vT) whose endpoints are on the
		//     external face in some upward-drawing of the sT-skeleton of vT
		//   for each unmarked edge e of skeleton(vT), constrain the tree
		//     edge associated with e to be directed towards vT
		//   if the source is not in skeleton(vT), let wT be the node neighbour
		//     of vT whose pertinet digraph contains the source, and constrain
		//     the tree edge (vT.wT) to be directed towards wT
		//
		// Determine whether T can be rooted at a Q-node in such a way that
		// orienting edges from children to parents satisfies the constraints
		// above. Procedure directSkeletons() returns true iff such a rooting
		// exists
		edge eRoot = directSkeletons(T, skInfo);


		if (eRoot == nullptr) {
			return false;
		}

#ifdef OGDF_DEBUG
		exp.consistencyCheck();
#endif

		if (embed) {
			T.rootTreeAt(eRoot);

			embedSkeleton(exp, T, skInfo, T.rootNode(), true);

			T.embed(exp);
#ifdef OGDF_DEBUG
			exp.consistencyCheck();
#endif

			OGDF_ASSERT(exp.representsCombEmbedding());
			CombinatorialEmbedding E(exp);

			FaceSinkGraph F(E, s);

			// find possible external faces (the faces in T containing s)
			//externalFaces.clear();
			SList<face> externalFaces;
			F.possibleExternalFaces(externalFaces);

			OGDF_ASSERT(!externalFaces.empty());

			face extFace = externalFaces.front();

			NodeArray<face> assignedFace(exp, nullptr);
			assignSinks(F, extFace, assignedFace);

			adjEntry adj1 = nullptr;
			for (adjEntry adj : s->adjEntries) {
				if (E.leftFace(adj) == extFace) {
					adj1 = adj;
					break;
				}
			}

			OGDF_ASSERT(adj1 != nullptr);

			// handle (single) source
			adjacentEdges[sG].pushBack(exp.original(adj1->theEdge())->adjSource());

			for (adjEntry adj = adj1->cyclicSucc(); adj != adj1; adj = adj->cyclicSucc()) {
				adjacentEdges[sG].pushBack(exp.original(adj->theEdge())->adjSource());
			}

			// handle internal vertices
			for (edge e : exp.edges) {
				if (exp.original(e)) {
					continue;
				}

				node vG = exp.original(e->source());
				adj1 = e->adjSource();
				for (adjEntry adj = adj1->cyclicSucc(); adj != adj1; adj = adj->cyclicSucc()) {
					adjacentEdges[vG].pushBack(exp.original(adj->theEdge())->adjTarget());
				}

				adj1 = e->adjTarget();
				for (adjEntry adj = adj1->cyclicSucc(); adj != adj1; adj = adj->cyclicSucc()) {
					adjacentEdges[vG].pushBack(exp.original(adj->theEdge())->adjSource());
				}
			}

			// handle sinks
			for (node v : exp.nodes) {
				if (v->outdeg() > 0) {
					continue;
				}
				node vG = exp.original(v);

				adj1 = nullptr;
				for (adjEntry adj : v->adjEntries) {
					if (E.leftFace(adj) == assignedFace[v]) {
						adj1 = adj;
						break;
					}
				}

				OGDF_ASSERT(adj1 != nullptr);

				adjacentEdges[vG].pushBack(exp.original(adj1->theEdge())->adjTarget());

				for (adjEntry adj = adj1->cyclicSucc(); adj != adj1; adj = adj->cyclicSucc()) {
					adjacentEdges[vG].pushBack(exp.original(adj->theEdge())->adjTarget());
				}
			}
		}

		// for each cut-vertex, process all components different from i
		SListPure<node> origNodes;
		for (node v : exp.nodes) {
			node vG = exp.original(v);

			if (vG && vG != sG) {
				origNodes.pushBack(vG);
			}
		}

		SListConstIterator<node> itV;
		for (itV = origNodes.begin(); itV.valid(); ++itV) {
			bool testOk = testBiconnectedComponent(exp, *itV, i, embed, adjacentEdges);
			if (!testOk) {
				return false;
			}
		}
	}

	return true;
}

// checks if precomputed in-/outdegrees in pertinet graphs are correctly
// (for debugging only)
bool UpwardPlanaritySingleSource::checkDegrees(SPQRTree& T, node s, NodeArray<SkeletonInfo>& skInfo) {
	const Graph& tree = T.tree();

	for (node vT : tree.nodes) {
		T.rootTreeAt(vT);

		const Skeleton& S = T.skeleton(vT);
		const Graph& M = S.getGraph();

		for (edge e : M.edges) {
			node wT = S.twinTreeNode(e);
			if (wT == nullptr) {
				continue;
			}

			PertinentGraph P;
			T.pertinentGraph(wT, P);

			Graph& Gp = P.getGraph();

			if (P.referenceEdge()) {
				Gp.delEdge(P.referenceEdge());
			}

			node x = nullptr, y = nullptr;
			for (node v : Gp.nodes) {
				if (P.original(v) == S.original(e->source())) {
					x = v;
				}
				if (P.original(v) == S.original(e->target())) {
					y = v;
				}
			}

			OGDF_ASSERT(x != nullptr);
			OGDF_ASSERT(y != nullptr);

			const DegreeInfo& degInfo = skInfo[vT].m_degInfo[e];
			if (x->indeg() != degInfo.m_indegSrc) {
				return false;
			}
			if (x->outdeg() != degInfo.m_outdegSrc) {
				return false;
			}
			if (y->indeg() != degInfo.m_indegTgt) {
				return false;
			}
			if (y->outdeg() != degInfo.m_outdegTgt) {
				return false;
			}

			bool contSource = false;
			for (node v : Gp.nodes) {
				if (v != x && v != y && P.original(v) == s) {
					contSource = true;
				}
			}

			if (skInfo[vT].m_containsSource[e] != contSource) {
				return false;
			}
		}
	}

	return true;
}

// precompute information: in-/outdegrees in pertinent graph, contains
// pertinent graph the source?
void UpwardPlanaritySingleSource::computeDegreesInPertinent(const SPQRTree& T, node s,
		NodeArray<SkeletonInfo>& skInfo, node vT) {
	const Skeleton& S = T.skeleton(vT);
	const Graph& M = S.getGraph();
	EdgeArray<DegreeInfo>& degInfo = skInfo[vT].m_degInfo;
	EdgeArray<bool>& containsSource = skInfo[vT].m_containsSource;

	// recursively compute in- and outdegrees at virtual edges except for
	// the reference edge of S; additionally compute if source is contained
	// in pertinent (and no endpoint of reference edge)
	for (adjEntry adj : vT->adjEntries) {
		edge eT = adj->theEdge();
		node wT = eT->target();
		if (wT != vT) {
			computeDegreesInPertinent(T, s, skInfo, wT);
		}
	}

	edge eRef = S.referenceEdge();
	node src = eRef->source();
	node tgt = eRef->target();

	bool contSource = false;
	for (node v : M.nodes) {
		if (v != src && v != tgt && S.original(v) == s) {
			contSource = true;
		}
	}

	// the in- and outdegree at real edges is obvious ...
	// m_containsSource is set to false by default
	for (edge e : M.edges) {
		if (!S.isVirtual(e)) {
			degInfo[e].m_indegSrc = 0;
			degInfo[e].m_outdegSrc = 1;
			degInfo[e].m_indegTgt = 1;
			degInfo[e].m_outdegTgt = 0;
		} else if (e != eRef) {
			contSource |= containsSource[e];
		}
	}


	if (vT == T.rootNode()) {
		return; // no virtual reference edge
	}


	// compute in- and outdegree of poles of reference edge in pertinent graph
	// of reference edge

	int indegSrc = 0;
	int outdegSrc = 0;
	{
		for (adjEntry adj : src->adjEntries) {
			edge e = adj->theEdge();
			if (e == eRef) {
				continue;
			}
			if (e->source() == src) {
				indegSrc += degInfo[e].m_indegSrc;
				outdegSrc += degInfo[e].m_outdegSrc;
			} else {
				indegSrc += degInfo[e].m_indegTgt;
				outdegSrc += degInfo[e].m_outdegTgt;
			}
		}
	}

	int indegTgt = 0;
	int outdegTgt = 0;
	{
		for (adjEntry adj : tgt->adjEntries) {
			edge e = adj->theEdge();
			if (e == eRef) {
				continue;
			}
			if (e->source() == tgt) {
				indegTgt += degInfo[e].m_indegSrc;
				outdegTgt += degInfo[e].m_outdegSrc;
			} else {
				indegTgt += degInfo[e].m_indegTgt;
				outdegTgt += degInfo[e].m_outdegTgt;
			}
		}
	}

	// set degrees at reference edge
	node srcOrig = S.original(src);
	degInfo[eRef].m_indegSrc = srcOrig->indeg() - indegSrc;
	degInfo[eRef].m_outdegSrc = srcOrig->outdeg() - outdegSrc;

	node tgtOrig = S.original(tgt);
	degInfo[eRef].m_indegTgt = tgtOrig->indeg() - indegTgt;
	degInfo[eRef].m_outdegTgt = tgtOrig->outdeg() - outdegTgt;

	containsSource[eRef] = !contSource && s != S.original(src) && s != S.original(tgt);

	// set degres at twin edge of reference edge
	node wT = S.twinTreeNode(eRef);
	DegreeInfo& degInfoTwin = skInfo[wT].m_degInfo[S.twinEdge(eRef)];
	degInfoTwin.m_indegSrc = indegSrc;
	degInfoTwin.m_outdegSrc = outdegSrc;
	degInfoTwin.m_indegTgt = indegTgt;
	degInfoTwin.m_outdegTgt = outdegTgt;

	skInfo[wT].m_containsSource[S.twinEdge(eRef)] = contSource;
}

// by default, corresponding virtual edges should be oriented in the
// same directions, i.e., the original nodes of source/target are equal
// this procedure serves to assert that this really holds!
bool UpwardPlanaritySingleSource::virtualEdgesDirectedEqually(const SPQRTree& T) {
	for (node v : T.tree().nodes) {
		const Skeleton& S = T.skeleton(v);
		const Graph& M = S.getGraph();

		for (edge e : M.edges) {
			edge eTwin = S.twinEdge(e);

			if (eTwin == nullptr) {
				continue;
			}

			const Skeleton& STwin = T.skeleton(S.twinTreeNode(e));

			if (S.original(e->source()) != STwin.original(eTwin->source())) {
				return false;
			}

			if (S.original(e->target()) != STwin.original(eTwin->target())) {
				return false;
			}
		}
	}

	return true;
}

// compute sT-skeletons
// test for upward-planarity, build constraints for rooting, and find a
// rooting of the tree satisfying all constraints
// returns true iff such a rooting exists
edge UpwardPlanaritySingleSource::directSkeletons(SPQRTree& T, NodeArray<SkeletonInfo>& skInfo) {
	const Graph& tree = T.tree();
	ConstraintRooting rooting(T);

	// we assume that corresponding virtual edges are directed equally before,
	// i.e., the original nodes of the sources are equal and the original nodes
	// of the targets are equal
	OGDF_ASSERT(virtualEdgesDirectedEqually(T));

	for (node vT : tree.nodes) {
		const StaticSkeleton& S = *dynamic_cast<StaticSkeleton*>(&T.skeleton(vT));
		const Graph& M = S.getGraph();

		for (edge e : M.edges) {
			edge eTwin = S.twinEdge(e);
			if (eTwin == nullptr) {
				continue;
			}

			const DegreeInfo& degInfo = skInfo[vT].m_degInfo[e];
			bool contSource = skInfo[vT].m_containsSource[e];

			node u = e->source();
			node v = e->target();
			// K = pertinent subgraph of e
			// K^0 = K - {u,v}

			const DegreeInfo& degInfoTwin = skInfo[S.twinTreeNode(e)].m_degInfo[eTwin];

			bool uTwinIsSource = (degInfoTwin.m_indegSrc == 0);
			bool vTwinIsSource = (degInfoTwin.m_indegTgt == 0);

			// Rule 1
			//  u and v are sources of K
			if (degInfo.m_indegSrc == 0 && degInfo.m_indegTgt == 0) {
				// replace e by a peak
				T.replaceSkEdgeByPeak(vT, e);

				// Rule 2
				// u is a source of K and v is a sink of K
			} else if (degInfo.m_indegSrc == 0 && degInfo.m_outdegTgt == 0) {
				// s not in K^0 ?
				if (!contSource) {
					// a directed edge (u,v)
					T.directSkEdge(vT, e, u);
				} else {
					// a peak
					T.replaceSkEdgeByPeak(vT, e);
				}

				// Rule 2'
				// v is a source of K and u is a sink of K
			} else if (degInfo.m_indegTgt == 0 && degInfo.m_outdegSrc == 0) {
				// s not in K^0 ?
				if (!contSource) {
					// a directed edge (v,u)
					T.directSkEdge(vT, e, v);
				} else {
					// a peak
					T.replaceSkEdgeByPeak(vT, e);
				}


				// Rule 3
				// u is a source of K and v is an internal vertex of K
			} else if (degInfo.m_indegSrc == 0 && degInfo.m_indegTgt > 0 && degInfo.m_outdegTgt > 0) {
				// v is a source of G-K and s not in K^0
				if (vTwinIsSource && !contSource) {
					// a directed edge (u,v)
					T.directSkEdge(vT, e, u);
				} else {
					// a peak
					T.replaceSkEdgeByPeak(vT, e);
				}

				// Rule 3'
				// v is a source of K and u is an internal vertex of K
			} else if (degInfo.m_indegTgt == 0 && degInfo.m_indegSrc > 0 && degInfo.m_outdegSrc > 0) {
				// u is a source of G-K and s not in K^0
				if (uTwinIsSource && !contSource) {
					// a directed edge (v,u)
					T.directSkEdge(vT, e, v);
				} else {
					// a peak
					T.replaceSkEdgeByPeak(vT, e);
				}

				// Rule 4
				// u and v ar not source of K
			} else {
				// u is a source of G-K
				if (uTwinIsSource) {
					// a directed edge (u,v)
					T.directSkEdge(vT, e, u);
				} else {
					// a directed edge (v,u)
					T.directSkEdge(vT, e, v);
				}
			}
		}


		// at this point, the sT-skeleton of vT is computed

		// if the source is not in skeleton(vT), let eWithSource be the virtual
		// edge in skeleton(vT) whose expansion graph contains the source
		edge eWithSource = nullptr;
		for (edge e : M.edges) {
			if (skInfo[vT].m_containsSource[e]) {
				eWithSource = e;
				break;
			}
		}

		// constrain the tree edge associated with eWithSource to be directed
		// leaving vT
		if (eWithSource && !rooting.constrainTreeEdge(S.treeEdge(eWithSource), vT)) {
			return nullptr;
		}

		// test whether the sT-skeletonof vT is upward planar and determine
		// the possible external faces
		if (!isAcyclic(M)) {
			return nullptr;
		}

		// for S- and P-nodes, we know already that they are upward planar
		if (T.typeOf(vT) != SPQRTree::NodeType::RNode) {
			continue;
		}

		//FaceSinkGraph &F = skInfo[vT].m_F;
		//ConstCombinatorialEmbedding &E = skInfo[vT].m_E;
		SList<face>& externalFaces = skInfo[vT].m_externalFaces;

		if (!initFaceSinkGraph(M, skInfo[vT])) {
			return nullptr;
		}


		// mark the edges of the skeleton of vT whose endpoints are on the
		// external face of some upward drawing of the sT-skeleton of vT
		EdgeArray<bool> marked(M, false);

		for (face f : externalFaces) {
			for (adjEntry adj : f->entries) {
				marked[adj] = true;
			}
		}

		// for each unmarked edge e of the skeleton of vT, constrain the tree
		// edge associated with e to be directed towards vT
		for (edge e : M.edges) {
			if (marked[e]) {
				continue;
			}

			edge eT = S.treeEdge(e);
			edge eR = S.realEdge(e);
			if (eR) {
				rooting.constrainRealEdge(eR);
			} else if (eT && !rooting.constrainTreeEdge(eT, S.twinTreeNode(e))) {
				return nullptr;
			}
		}
	}


	// determine whether T can be rooted at a Q-node in such a way that
	// orienting edges from children to parents satisfies the constraints
	// above
	//rooting.outputConstraints(std::cout);
	edge eRoot = rooting.findRooting();

	//std::cout << "\nroot edge: " << eRoot << std::endl;

	return eRoot;
}

// initializes embedding and face-sink graph in skeleton info for
// embedded skeleton graph
// determines the possible external faces and returns true if skeleton
// graph is upward planar
bool UpwardPlanaritySingleSource::initFaceSinkGraph(const Graph& M, SkeletonInfo& skInfo) {
	ConstCombinatorialEmbedding& E = skInfo.m_E;
	FaceSinkGraph& F = skInfo.m_F;
	SList<face>& externalFaces = skInfo.m_externalFaces;

	E.init(M);

	node s;
	hasSingleSource(M, s);
	OGDF_ASSERT(s != nullptr);

	F.init(E, s);

	// find possible external faces (the faces in T containing s)
	F.possibleExternalFaces(externalFaces);

	return !externalFaces.empty();
}

// embeds skeleton(vT) and mirrors it if necessary such that the external
// face is to the left of the reference edge iff extFaceIsLeft = true
void UpwardPlanaritySingleSource::embedSkeleton(Graph& G, StaticPlanarSPQRTree& T,
		NodeArray<SkeletonInfo>& skInfo, node vT, bool extFaceIsLeft) {
	StaticSkeleton& S = *dynamic_cast<StaticSkeleton*>(&T.skeleton(vT));
	Graph& M = S.getGraph();

	edge eRef = S.referenceEdge();

	// We still have to embed skeletons of P-nodes
	if (T.typeOf(vT) == SPQRTree::NodeType::PNode) {
		node lowerPole = eRef->source();

		SListPure<adjEntry> lowerAdjs, upperAdjs;

		adjEntry adjRef = eRef->adjSource();
		adjEntry adjRefTwin = adjRef->twin();

		for (adjEntry adj : lowerPole->adjEntries) {
			// ignore reference edge
			if (adj == adjRef) {
				continue;
			}

			adjEntry adjTwin = adj->twin();
			if (S.original(adjTwin->theNode()) == nullptr) { // peak
				lowerAdjs.pushFront(adj);
				upperAdjs.pushBack(adjTwin->cyclicSucc()->twin());

			} else { // non-peak
				lowerAdjs.pushBack(adj);
				upperAdjs.pushFront(adjTwin);
			}
		}

		// adjacency entries in lowerAdjs are now sorted: peaks, non-peaks
		// (without reference edge)

		// is reference edge a peak
		bool isRefPeak = (S.original(eRef->target()) == nullptr);
		adjEntry adjRefUpper = isRefPeak ? adjRefTwin->cyclicSucc()->twin() : adjRefTwin;

		if (isRefPeak) {
			// lowerAdjs: ref. peak, other peaks, non-peaks
			lowerAdjs.pushFront(adjRef);
			upperAdjs.pushBack(adjRefUpper);
		} else {
			// lowerAdjs: peaks, non-peaks, ref. non-peak
			lowerAdjs.pushBack(adjRef);
			upperAdjs.pushFront(adjRefUpper);
		}


		M.sort(lowerPole, lowerAdjs);
		M.sort(adjRefUpper->theNode(), upperAdjs);
	}


	OGDF_ASSERT(M.representsCombEmbedding());

	// we still have to compute the face-sink graph for S and P nodes
	if (T.typeOf(vT) != SPQRTree::NodeType::RNode) {
#ifdef OGDF_DEBUG
		bool testOk =
#endif
				initFaceSinkGraph(M, skInfo[vT]);
		OGDF_ASSERT(testOk);
	}

	// variables stored in skeleton info
	ConstCombinatorialEmbedding& E = skInfo[vT].m_E;
	FaceSinkGraph& F = skInfo[vT].m_F;
	SList<face>& externalFaces = skInfo[vT].m_externalFaces;


	// determine a possible external face extFace adjacent to eRef
	face fLeft = E.leftFace(eRef->adjSource());
	face fRight = E.rightFace(eRef->adjSource());
	face extFace = nullptr;

	SListConstIterator<face> it;
	for (it = externalFaces.begin(); it.valid(); ++it) {
		face f = *it;
		if (f == fLeft || f == fRight) {
			extFace = f;
		}
	}

	OGDF_ASSERT(extFace != nullptr);

	// possibly we have to mirror the embedding of S
	bool mirrorEmbedding = (extFaceIsLeft != (extFace == fLeft));

	// assign sinks to faces (a sink t is assigned to the face above t in
	// an upward planar drawing)
	NodeArray<face> assignedFace(M, nullptr);
	assignSinks(F, extFace, assignedFace);


	// Traverse SPQR-tree by depth-first search.
	// This kind of traversal is required for correctly embedding the
	// skeletons (see paper).
	for (edge e : M.edges) {
		edge eT = S.treeEdge(e);
		if (eT == nullptr) {
			continue;
		}

		node wT = eT->target();
		if (wT == vT) {
			continue;
		}

		bool eExtFaceLeft = true;
		node p = e->target();
		if (S.original(p) == nullptr) {
			const Skeleton& S2 = T.skeleton(wT);
			node vOrig = S2.original(S2.referenceEdge()->source());
			adjEntry adj = e->adjSource();

			if (S.original(e->source()) != vOrig) {
				adj = adj->twin()->cyclicSucc()->twin();
			}

			if (assignedFace[p] == E.rightFace(adj)) {
				eExtFaceLeft = true;
			} else {
				OGDF_ASSERT(assignedFace[p] == E.leftFace(adj));
				eExtFaceLeft = false;
			}
			if (mirrorEmbedding) {
				eExtFaceLeft = !eExtFaceLeft;
			}
		}

		// recursive call
		embedSkeleton(G, T, skInfo, wT, eExtFaceLeft);
	}

	if (mirrorEmbedding) {
		T.reverse(vT);
	}

#ifdef OGDF_HEAVY_DEBUG
	M.consistencyCheck();
#endif


	// We have to unsplit all peaks in skeleton S again, since the embedding
	// procedure of StaticPlanarSPQRTree does not support such modified skeletons
	// (only reversing skeleton edges is allowed).
	node v, vSucc;
	for (v = M.firstNode(); v; v = vSucc) {
		vSucc = v->succ();

		if (S.original(v) == nullptr) {
			OGDF_ASSERT(v->indeg() == 2);
			OGDF_ASSERT(v->outdeg() == 0);

			edge e1 = v->firstAdj()->theEdge();
			edge e2 = v->lastAdj()->theEdge();

			if (S.realEdge(e1) != nullptr || S.twinEdge(e1) != nullptr) {
				M.reverseEdge(e2);
			} else {
				OGDF_ASSERT(S.realEdge(e2) != nullptr || S.twinEdge(e2) != nullptr);
				M.reverseEdge(e1);
			}

			M.unsplit(v);
		}
	}
}

// assigns each sink to a face such that there exists an upward planar drawing
// with exteriour fac extFace in which the assigned face of each sink t is
// above t
void UpwardPlanaritySingleSource::assignSinks(FaceSinkGraph& F, face extFace,
		NodeArray<face>& assignedFace) {
	// find the representative h of face extFace in face-sink graph F
	node h = nullptr;
	for (node v : F.nodes) {
		if (F.originalFace(v) == extFace) {
			h = v;
			break;
		}
	}
	OGDF_ASSERT(h != nullptr);

	// find all roots of trees in F different from (the unique tree T)
	// rooted at h
	SListPure<node> roots;
	for (node v : F.nodes) {
		node vOrig = F.originalNode(v);
		if (vOrig != nullptr && vOrig->indeg() > 0 && vOrig->outdeg() > 0) {
			roots.pushBack(v);
		}
	}

	// recursively assign faces to sinks
	dfsAssignSinks(F, h, nullptr, assignedFace);

	SListConstIterator<node> itRoots;
	for (itRoots = roots.begin(); itRoots.valid(); ++itRoots) {
		dfsAssignSinks(F, *itRoots, nullptr, assignedFace);
	}
}

node UpwardPlanaritySingleSource::dfsAssignSinks(FaceSinkGraph& F, // face-sink graph
		node v, // current node
		node parent, // its parent
		NodeArray<face>& assignedFace) {
	bool isFace = (F.originalFace(v) != nullptr);
	node vf = nullptr;

	// we perform a dfs-traversal (underlying graph is a tree)
	for (adjEntry adj : v->adjEntries) {
		node w = adj->twinNode();

		if (w == parent) {
			continue;
		}

		if (isFace) {
			assignedFace[F.originalNode(w)] = F.originalFace(v);
		}

		dfsAssignSinks(F, w, v, assignedFace);
	}

	return vf;
}

}
