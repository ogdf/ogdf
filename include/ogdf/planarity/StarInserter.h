/** \file
 * \brief Declaration of class StarInserter.
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

#pragma once

#include <memory>
#include <unordered_map>
#include <ogdf/basic/DualGraph.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/comparer.h>

namespace ogdf {

using PredecessorMap = std::unordered_map<node, std::unique_ptr<NodeArray<edge>>>;

/**
 * Orders edges such that they do not cross each other when embeddded as
 * insertion paths. In particular, multiple edges e1...ek splitting the same
 * edge are ordered such that no additional crossings arise.
 *
 * The comparison uses the insertion paths of the edges (overall forming a tree)
 * as a basis: These insertion paths pass through several faces; the face that
 * serves as the lowest common ancestor of both edges is observed.
 * The edges are ordered clockwise around the lca face, starting at the edge to
 * the parent of the lca in the insertion path tree.
 */
class EdgeOrderComparer {
private:
	node m_origNode; //!< Common (original) node of compared edges.

	node m_rootDualNode; //!< Dual node, root of the insertion path tree.

	//! Insertion path tree, given by several insertion paths.
	/**
	 * The array is indexed by the copy node corresponding to the opposite of
	 * (the copy of) #m_origNode relative to the inserted edge.
	 * An insertion path is given by the predecessor edge of each node on a
	 * path from #m_rootDualNode to one leaf the insertion path tree.
	 */
	const PredecessorMap &m_predecessors;

	const GraphCopy &m_graphCopy; //!< Planarization.

	const DynamicDualGraph *m_dualGraph; //!< Dual graph of #m_graphCopy.

public:
	/**
	 * Creates a new EdgeOrderComparer.
	 *
	 * @param origNode Common (original) node of compared edges.
	 * @param rootDualNode Dual node, root of the insertion path tree.
	 * @param predecessors Insertion path tree given by several insertion paths.
	 * @param graphCopy Planarization.
	 * @param dualGraph Dual graph of #m_graphCopy.
	 */
	EdgeOrderComparer(
			node origNode,
			node rootDualNode,
			const PredecessorMap &predecessors,
			const GraphCopy &graphCopy,
			const DynamicDualGraph *dualGraph)
	: m_origNode(origNode)
	, m_rootDualNode(rootDualNode)
	, m_predecessors(predecessors)
	, m_graphCopy(graphCopy)
	, m_dualGraph(dualGraph)
	{}

	/**
	 * Returns the lowest common ancestor of \p w1 and \p w2 in the insertion
	 * path tree given by #m_predecessors.
	 *
	 * @param w1 node incident to the first leaf of the insertion path tree
	 * @param w2 node incident to the second leaf of the insertion path tree
	 * @param parentOfLCA is assigned the edge from the parent of the lca to the
	 * lca in the insertion path tree
	 * @return lowest common ancestor of \p w1 and \p w2 in the insertion path
	 * tree rooted at \p root
	 */
	node findLCAInInsertionPathTree(node w1, node w2, edge &parentOfLCA) const {
		const NodeArray<edge> &preds1 {*m_predecessors.at(w1)};
		const NodeArray<edge> &preds2 {*m_predecessors.at(w2)};

		// Go backwards from root until predecessor path diverges.
		node lastDualNode1 {m_rootDualNode};
		node lastDualNode2 {m_rootDualNode};
		edge edgeToChild {nullptr};
		node lca {nullptr};
		while (lastDualNode1 == lastDualNode2) {
			// Remember current lca.
			parentOfLCA = edgeToChild;
			lca = lastDualNode1;

			// Go further down.
			edgeToChild = preds1[lastDualNode1];
			if (preds1[lastDualNode1] == nullptr ||
			    preds2[lastDualNode2] == nullptr) {
				// If either node has no child in the insertion path tree,
				// the lca cannot be lower.
				break;
			} else {
				lastDualNode1 = preds1[lastDualNode1]->opposite(lastDualNode1);
				lastDualNode2 = preds2[lastDualNode2]->opposite(lastDualNode2);
			}
		}

		return lca;
	}

	/**
	 * Compares two edges as described for EdgeOrderComparer.
	 *
	 * @param e1 first edge to compare
	 * @param e2 second edge to compare
	 */
	int compare(const edge &e1, const edge &e2) const {
		node w1 {m_graphCopy.copy(e1->opposite(m_origNode))};
		node w2 {m_graphCopy.copy(e2->opposite(m_origNode))};
		if (w1 == w2) {
			return 0;
		}

		// Find lowest common ancestor in predecessor tree.
		edge parentOfLCA {nullptr};
		node lcaDualNode {findLCAInInsertionPathTree(w1, w2, parentOfLCA)};

		// w1 and w2 must have a common ancestor.
		OGDF_ASSERT(lcaDualNode != nullptr);

		// Get the adjEntry of the primal edge of parentOfLCA which is incident
		// to the face corresponding to lcaDualNode.
		// This adjEntry marks the start of the cyclic order around lcaDualNode.
		face lcaFace {m_dualGraph->primalFace(lcaDualNode)};
		adjEntry parentAdj {nullptr};
		if (parentOfLCA == nullptr) {
			OGDF_ASSERT(lcaDualNode == m_rootDualNode);
			// The lca is the root of the insertion path tree, hence it has no
			// predecessor. Just use any adjEntry to start the cyclic order.
			parentAdj = lcaFace->firstAdj();
		} else {
			// The lca is not the root of the insertion path tree.
			// Determine the correct adjEntry of the primal of its parent-edge.
			edge primalParentOfLCA {m_dualGraph->primalEdge(parentOfLCA)};
			adjEntry sourceAdj {primalParentOfLCA->adjSource()};
			adjEntry targetAdj {primalParentOfLCA->adjTarget()};

			const ConstCombinatorialEmbedding &embedding =
				m_dualGraph->getPrimalEmbedding();
			if (embedding.rightFace(sourceAdj) == lcaFace) {
				parentAdj = sourceAdj;
			} else {
				OGDF_ASSERT(embedding.rightFace(targetAdj) == lcaFace);
				parentAdj = targetAdj;
			}
		}

		// Traverse adjEntries of face corrsponding to lcaDualNode clockwise,
		// starting at parent edge of lcaDualNode in the insertion path tree.
		auto cyclicNextAdj = [&parentAdj](adjEntry adj) {
			adjEntry nextAdj = adj->faceCycleSucc();
			return nextAdj == parentAdj ? nullptr : nextAdj;
		};

		edge pred1 {(*m_predecessors.at(w1))[lcaDualNode]};
		edge pred2 {(*m_predecessors.at(w2))[lcaDualNode]};
		for (adjEntry adj {parentAdj}; adj != nullptr; adj = cyclicNextAdj(adj)) {
			edge primalEdge {adj->theEdge()};
			node primalNode {adj->theNode()};
			// First, if lcaFace has no pred in the insertion path of w1/w2,
			// w1/w2 is adjacent to lcaFace.
			// Check whether the node of the current adj is w1/w2.
			if (pred1 == nullptr && primalNode == w1) {
				return -1;
			}
			if (pred2 == nullptr && primalNode == w2) {
				return 1;
			}

			// If w1/w2 is not adjacent to lcaFace, check whether the insertion
			// path from lcaFace to the face adjacent to w1/w2 crosses the
			// current edge.
			edge dualEdge = m_dualGraph->dualEdge(primalEdge);
			if (dualEdge == pred1) {
				OGDF_ASSERT(dualEdge != pred2);
				return -1;
			}
			if (dualEdge == pred2) {
				OGDF_ASSERT(dualEdge != pred1);
				return 1;
			}
		}

		// Something went terribly wrong: The LCA of the insertion paths of w1
		// and w2 is not actually part of the insertion paths of w1 and w2?
		OGDF_ASSERT(false);
		return 0;
	}

	OGDF_AUGMENT_COMPARER(edge)
};


/**
 * @ingroup ga-insert
 *
 * Inserts a star (a vertex and its incident edges) optimally into an
 * embedding.
 */
class OGDF_EXPORT StarInserter
{
private:
	GraphCopy *m_graphCopy; //!< Planarization.

	CombinatorialEmbedding *m_combEmbedding; //!< Embedding of #m_graphCopy.

	DynamicDualGraph *m_dual; //!< Dual graph of #m_combEmbedding.

	//! Maps new faces (created after the insertion of edge paths) to old
	//! (non-split) faces. Computed insertion paths refer to the old faces.
	FaceArray<face> *m_newToOldFace;

	//! Maps each edge e originally contained in #m_graphCopy to the edge in the
	//! same chain which should be split next if the dual of e is part of an
	//! insertion path.
	EdgeArray<edge> *m_edgeInChainToSplit;

	//! Maps the parts of a chain in #m_graphCopy to the respective copy edge
	//! contained in #m_graphCopy before any insertion paths were embedded.
	EdgeArray<edge> *m_originalEdge;

	//! Returns the primal face of \p dualNode which existed before any new
	//! edges were inserted.
	inline face oldPrimalFace(node dualNode) {
		return (*m_newToOldFace)[m_dual->primalFace(dualNode)];
	}

	/**
	 * Computes the optimal dual node to insert \p origNode and the insertion
	 * paths of its incident edges.
	 *
	 * @warning The insertion paths, i.e. \p predecessors, for different nodes
	 * might cross each other.
	 *
	 * @param origNode original node that is inserted
	 * @param pCostOrig points to the cost of each edge in the original graph.
	 * @param predecessors is assigned the insertion path of each edge incident
	 * to \p origNode. The map is indexed by copy nodes corresponding to the
	 * opposite of \p origNode for each original edge. The predecessors are
	 * given as edges in the dual graph, starting at the returned dual node.
	 * @return the node of #m_dual in which \p origNode should be inserted.
	 */
	node getOptimalDualNode(node origNode,
			const EdgeArray<int> *pCostOrig,
			PredecessorMap &predecessors);

	/**
	 * Modify the insertion paths \p predecessors such that they do not cross
	 * each other anymore, i.e. such that they form an arborescence with
	 * \p optimalDualNode as its root.
	 *
	 * @param origNode original node that is inserted
	 * @param optimalDualNode dual node that \p origNode is inserted in
	 * @param predecessors is assigned the predecessor edge of each node in the
	 * insertion path arborescence. The predecessors will not cross each other.
	 */
	void makePredsConsistent(node origNode,
			node optimalDualNode,
			PredecessorMap &predecessors);

	/**
	 * Returns the adjEntry of \p primalEdgeToSplit whose left face corresponds
	 * to \p leftDualNode.
	 *
	 * @param primalEdgeToSplit primal edge whose adjEntry should be returned
	 * @param leftDualNode dual node corresponding to the primal face left of the
	 * returned adjEntry
	 * @return the adjEntry of \p primalEdgeToSplit whose left face corresponds
	 * to \p leftDualNode.
	 */
	adjEntry getCrossedAdjEntry(edge primalEdgeToSplit, node leftDualNode);

	/**
	 * Returns the adjEntry of \p primalNode whose right face corresponds to \p
	 * rightDualNode (if \p otherPrimaNode is isolated) or is incident to
	 * \p otherPrimalNode (otherwise).
	 *
	 * @pre \p primalNode and \p otherPrimalNode have a common face (or one of
	 * \p primalNode and \p otherPrimalNode is isolated). \p primalNode is
	 * incident to \p rightDualNode.
	 *
	 * @param primalNode primal node whose adjEntry should be returned
	 * @param rightDualNode dual node corresponding to the face which is right
	 * of the returned adjEntry
	 * @param otherPrimalNode primal node which shares a common face with \p
	 * primal node (or is isolated)
	 * @return adjEntry of \p primalNode whose right face is incident to
	 * \p primalEdge, or \c nullptr if \p primalNode is isolated.
	 * @returns adjEntry of \p primalNode whose right face corresponds to
	 * \p rightDualNode (if \p otherPrimaNode is isolated) or is incident to
	 * \p otherPrimalNode (otherwise), or \c nullptr if \p primalNode is
	 * isolated.
	 */
	adjEntry getAdjEntry(node primalNode, node rightDualNode, node otherPrimalNode);

	/**
	 * Returns the adjEntry of \p primalNode whose right face is incident to
	 * \p primalEdge.
	 *
	 * @pre \p primalNode and \p primalEdge have the common face corersponding
	 * to \p rightDualNode (or \p primalNode is isolated).
	 *
	 * @param primalNode primal node whose adjEntry should be returned
	 * @param rightDualNode dual node corresponding to the face which is right
	 * of the returned adjEntry
	 * @param primalEdge primal edge incident to the face that the returned
	 * adjEntry belongs to
	 * @param first determines whether the first or last possible adjEntry in
	 * clockwise order around \p primalNode is chosen if there are multiple ones
	 * which have the correct right face.
	 * @return adjEntry of \p primalNode whose right face is incident to
	 * \p primalEdge, or \c nullptr if \p primalNode is isolated.
	 */
	adjEntry getAdjEntry(node primalNode, node rightDualNode, edge primalEdge, bool first);

	/**
	 * Collects adjEntries which can be passed to
	 * GraphCopy::insertEdgePathEmbedded to embed the edge \p insertedNode -
	 * \p w.
	 *
	 * If \p insertedNode is isolated, a new dummy edge with \p insertedNode as
	 * its endpoint is created and returned. The edge has to be deleted after
	 * the insertion path is embedded.
	 *
	 * @param w the primal copy node opposite of \p insertedNode for the edge
	 * which is currently embedded
	 * @param insertedNode inserted primal node
	 * @param optimalDualNode dual node that \p insertedNode is inserted in
	 * @param predecessors insertion path tree given by several insertion paths.
	 * @param crossedEdges is assigned a List of adjEntries to be passed to
	 * GraphCopy::insertEdgePathEmbedded.
	 * @return newly created dummy edge incident to \p insertedNode (which has
	 * to be deleted later), or \c nullptr if no dummy edge was created.
	 */
	edge collectAdjEntries(
		node w,
		node insertedNode,
		node optimalDualNode,
		const PredecessorMap &predecessors,
		List<adjEntry> &crossedEdges);

	/**
	 * Transfers the adjEntries from the List \p crossedEdges to the SList
	 * \p finalCrossedEdges such that they can be passed to
	 * GraphCopy::insertEdgePathEmbedded.
	 *
	 * Depending on \p startAtSource, the list may be reversed in the process
	 * and the contained adjEntries for crossed edges may be reversed as well.
	 *
	 * @param crossedEdges initial list of adjEntries
	 * @param finalCrossedEdges is assigned the list of adjEntries that can be
	 * passed to GraphCopy::insertEdgePathEmbedded
	 * @param startAtSource whether the order of adjEntries in \p crossedEdges
	 * already conforms to the direction of the edge to embed (i.e. whether
	 * \p crossedEdges is ordered such that the first adjEntry starts at the
	 * source of the edge to embed)
	 */
	void transferCrossedEdges(
		const List<adjEntry> &crossedEdges,
		SList<adjEntry> &finalCrossedEdges,
		bool startAtSource);

	/**
	 * Initialize all member variables.
	 *
	 * @param graphCopy initial planarization
	 * @param dualGraph dual graph of \p graphCopy
	 */
	void initMemberData(GraphCopy &graphCopy, DynamicDualGraph &dualGraph);

	/**
	 * Update member variables, in particular #m_newToOldFace and
	 * #m_edgeInChainToSplit, after an edge was inserted.
	 *
	 * @param origEdge the inserted edge
	 * @param startAtSource whether the source of \p origEdge is the inserted
	 * origNode or not
	 */
	void updateMemberData(edge origEdge, bool startAtSource);

public:
	//! Creates a StarInserter instance with default settings.
	StarInserter()
		: m_combEmbedding{nullptr}
		, m_dual{nullptr}
	{ }

	//! Creates a StarInserter instance with the same settings as \p inserter.
	/**
	 * @param inserter StarInserter to be copied
	 */
	StarInserter(const StarInserter &inserter) { }

	//! Destructor.
	~StarInserter() { }

	//! Assignment operator, copies option settings only.
	StarInserter &operator=(const StarInserter &inserter);

	//! Inserts the node \p origNode and its incident edges into \p graphCopy.
	/**
	 * @param graphCopy is the input planarized representation and will also
	 * be assigned the result.
	 * @param dualGraph is the dual graph of \p graphCopy.
	 * @param origNode is the original node (in the original graph of
	 * \p graphCopy) to be inserted.
	 * @param pCostOrig points to an edge array containing the costs of original
	 * edges; edges in \p graphCopy without an original edge have zero costs.
	 * May be nullptr, in which case all edges have cost 1.
	 */
	virtual void call(GraphCopy &graphCopy,
		DynamicDualGraph &dualGraph,
		node origNode,
		const EdgeArray<int> *pCostOrig);
};

}
