/** \file
 * \brief Declaration of graph copy classes.
 *
 * \author Carsten Gutwenger, Max Ilsen, Simon D. Fink
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

#include <ogdf/basic/DualGraph.h>
#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/NodeArray.h>
#include <ogdf/basic/SList.h>

namespace ogdf {

template<bool>
class FaceSet;

class OGDF_EXPORT GraphCopyBase : public Graph {
protected:
	const Graph* m_pGraph = nullptr; //!< The original graph.
	NodeArray<node> m_vOrig; //!< The corresponding node in the original graph.
	EdgeArray<edge> m_eOrig; //!< The corresponding edge in the original graph.
	NodeArray<node> m_vCopy; //!< The corresponding node in the graph copy.
	bool m_linkCopiesOnInsert =
			true; //!< Whether \c insert(getOriginalGraph()) will set \c copy and \c original

public:
	//! Constructs a GraphCopyBase associated with no graph.
	GraphCopyBase() = default;

	GraphCopyBase(const GraphCopyBase& other) = delete;

	GraphCopyBase(GraphCopyBase&& other) noexcept = delete;

	GraphCopyBase& operator=(const GraphCopyBase& other) = delete;

	GraphCopyBase& operator=(GraphCopyBase&& other) noexcept = delete;

	//! Re-initializes the copy using \p G, creating copies for all nodes and edges in \p G.
	void init(const Graph& G) {
		Graph::clear();
		setOriginalGraph(&G);
		insert(G);
	}

	//! Re-initializes the copy using \p G (which might be null), creating copies for all nodes and edges in \p G.
	void init(const Graph* G) {
		Graph::clear();
		setOriginalGraph(G);
		if (G != nullptr) {
			insert(*G);
		}
	}

	//! Re-initializes the copy using \p G, but does not create any nodes or edges.
	OGDF_DEPRECATED("setOriginalGraph() should be used instead.")

	void createEmpty(const Graph& G) { setOriginalGraph(&G); }

	//! Re-initializes the copy using \p G (which might be null), but does not create any nodes or edges.
	virtual void setOriginalGraph(const Graph* G) = 0;

	//! Re-initializes the copy using \p G, but does not create any nodes or edges.
	void setOriginalGraph(const Graph& G) { setOriginalGraph(&G); };

	const Graph* getOriginalGraph() const { return m_pGraph; }

	//! Removes all nodes and edges from this copy but does not break the link with the original graph.
	void clear() override = 0;

	//! Returns a reference to the original graph.
	const Graph& original() const {
		OGDF_ASSERT(m_pGraph != nullptr);
		return *m_pGraph;
	}

	/**
	 * \brief Returns the node in the original graph corresponding to \p v.
	 * @param v is a node in the graph copy.
	 * \return the corresponding node in the original graph, or 0 if no
	 *         such node exists.
	 */
	node original(node v) const { return m_vOrig[v]; }

	/**
	 * \brief Returns the edge in the original graph corresponding to \p e.
	 * @param e is an edge in the graph copy.
	 * \return the corresponding edge in the original graph, or 0 if no
	 *         such edge exists.
	 */
	edge original(edge e) const { return m_eOrig[e]; }

	/**
	 * Returns the adjacency entry in the original graph corresponding to \p adj.
	 *
	 * Note that this method does not pay attention to reversed edges.
	 * Given a source (target) adjacency entry, the source (target) adjacency entry of the
	 * original edge is returned.
	 *
	 * @param adj is an adjacency entry in the copy graph.
	 * \return the corresponding adjacency entry in the original graph.
	 */
	adjEntry original(adjEntry adj) const {
		edge f = m_eOrig[adj->theEdge()];
		return adj->isSource() ? f->adjSource() : f->adjTarget();
	}

	/**
	 * \brief Returns the node in the graph copy corresponding to \p v.
	 * @param v is a node in the original graph.
	 * \return the corresponding node in the graph copy,
	 * or \c nullptr if it doesn't exist.
	 */
	node copy(node v) const { return m_vCopy[v]; }

	/**
	 * \brief Returns the edge in the graph copy corresponding to \p e.
	 * Has to be defined by the implemented GraphCopyBase subclass.
	 *
	 * @param e is an edge in the original graph.
	 * \return the corresponding edge in the graph copy,
	 * or \c nullptr if it doesn't exist.
	 */
	virtual edge copy(edge e) const = 0;

	/**
	 * \brief Returns the adjacency entry in the graph copy corresponding to \p adj.
	 * Has to be defined by the implemented GraphCopyBase subclass.
	 *
	 * @param adj is an adjacency entry in the original graph.
	 * \return the corresponding adjacency entry in the graph copy,
	 * or \c nullptr if it doesn't exist.
	 */
	virtual adjEntry copy(adjEntry adj) const = 0;

	/**
	 * \brief Returns true iff \p v has no corresponding node in the original graph.
	 * @param v is a node in the graph copy.
	 */
	bool isDummy(node v) const { return m_vOrig[v] == nullptr; }

	/**
	 * \brief Returns true iff \p e has no corresponding edge in the original graph.
	 * @param e is an edge in the graph copy.
	 */
	bool isDummy(edge e) const { return m_eOrig[e] == nullptr; }

	/**
	 * \brief Returns true iff \p adj->theEdge() has no corresponding edge in the original graph.
	 * @param adj is an adjEntry in the graph copy.
	 */
	bool isDummy(adjEntry adj) const { return isDummy(adj->theEdge()); }

	/**
	 * \brief Creates a new node in the graph copy with original node \p vOrig.
	 * \warning You have to make sure that the original node makes sense, in
	 *   particular that \p vOrig is not the original node of another node in the copy.
	 */
	node newNode(node vOrig) {
		OGDF_ASSERT(vOrig != nullptr);
		OGDF_ASSERT(vOrig->graphOf() == m_pGraph);
		node v = Graph::newNode();
		m_vCopy[m_vOrig[v] = vOrig] = v;
		return v;
	}

	using Graph::newNode;

	/**
	 * \brief Removes node \p v.
	 *
	 * \param v is a node in the graph copy.
	 */
	void delNode(node v) override {
		node vOrig = m_vOrig[v];
		Graph::delNode(v);
		if (vOrig != nullptr) {
			m_vCopy[vOrig] = nullptr;
		}
	}

	//! Sets the embedding of the graph copy to the embedding of the original graph.
	/**
	 * Edges that have no copy in this graph will be left out, while dummy edges that are not present
	 * in the original graph will be moved to the end of their nodes' adjacency lists.
	 */
	virtual void setOriginalEmbedding() = 0;

	bool getLinkCopiesOnInsert() const { return m_linkCopiesOnInsert; }

	//! Whether \c insert(getOriginalGraph()) will automatically set \c copy and \c original
	/**
	 * Whether the inserted elements should automatically have assigned \c copy and \c original
	 * values when calling insert() with nodes and edges from getOriginalGraph().
	 * Note that this also applies to elements inserted when calling init().
	 *
	 * @param linkCopiesOnInsert When true, \c copy and \c original will be automatically set for
	 *   elements from the original graph. When false, all inserted elements (no matter from which
	 *   Graph) will be dummies.
	 */
	void setLinkCopiesOnInsert(bool linkCopiesOnInsert) {
		m_linkCopiesOnInsert = linkCopiesOnInsert;
	}

protected:
	void* preInsert(bool copyEmbedding, bool copyIDs, bool notifyObservers, bool edgeFilter,
			NodeArray<node>& nodeMap, EdgeArray<edge>& edgeMap, int* newNodes,
			int* newEdges) override;

	void nodeInserted(void* userData, node original, node copy) override;
};

/**
 * \brief Copies of graphs with mapping between nodes and edges
 *
 * @ingroup graphs
 *
 * The class GraphCopySimple represents a copy of a graph and
 * maintains a mapping between the nodes and edges of the original
 * graph to the copy and vice versa.
 *
 * New nodes and edges can be added to the copy; the counterpart
 * of those nodes and edges is 0 indicating that there is no counterpart.
 * This class <b>does not</b> support splitting of edges in such a way
 * that both edges resulting from the split are mapped to the same
 * original edge; this feature is provided by GraphCopy.
 */
class OGDF_EXPORT GraphCopySimple : public GraphCopyBase {
public:
	EdgeArray<edge> m_eCopy; //!< The corresponding edge in the graph copy.

	explicit GraphCopySimple() : GraphCopySimple(nullptr) { }

	explicit GraphCopySimple(const Graph& G) : GraphCopySimple(&G) { }

	explicit GraphCopySimple(const Graph* G) : GraphCopyBase() {
		if (G) {
			init(*G);
		}
	}

	GraphCopySimple(const GraphCopySimple& other) : GraphCopyBase() { *this = other; }

	GraphCopySimple& operator=(const GraphCopySimple& other);

	using GraphCopyBase::setOriginalGraph;

	//! Re-initializes the copy using \p G (which might be null), but does not create any nodes or edges.
	void setOriginalGraph(const Graph* G) override;

	//! Removes all nodes and edges from this copy but does not break the link with the original graph.
	void clear() override;

	using GraphCopyBase::copy;

	/**
	 * \brief Returns the edge in the graph copy corresponding to \p e.
	 * @param e is an edge in the original graph.
	 * \return the corresponding edge in the graph copy,
	 * or \c nullptr if it doesn't exist.
	 */
	edge copy(edge e) const override { return m_eCopy[e]; }

	/**
	 * Returns the adjacency entry in the graph copy corresponding to \p adj.
	 *
	 * Note that this method does not pay attention to reversed edges.
	 * Given a source (target) adjacency entry, the source (target) adjacency entry of the
	 * copy edge is returned.
	 *
	 * @param adj is an adjacency entry in the original graph.
	 * \return the corresponding adjacency entry in the graph copy,
	 * or \c nullptr if it doesn't exist.
	 */
	adjEntry copy(adjEntry adj) const override {
		edge f = m_eCopy[adj->theEdge()];
		if (f == nullptr) {
			return nullptr;
		}
		return adj->isSource() ? f->adjSource() : f->adjTarget();
	}

	/**
	 * \brief Creates a new edge in the graph copy with original edge \p eOrig.
	 * \warning You have to make sure that the original edge makes sense, in
	 *   particular that \p eOrig is not the original edge of another edge in the copy.
	 */
	edge newEdge(edge eOrig) {
		OGDF_ASSERT(eOrig != nullptr);
		OGDF_ASSERT(eOrig->graphOf() == m_pGraph);

		edge e = Graph::newEdge(m_vCopy[eOrig->source()], m_vCopy[eOrig->target()]);
		m_eCopy[m_eOrig[e] = eOrig] = e;

		return e;
	}

	using Graph::newEdge;

	/**
	 * \brief Removes edge \p e.
	 *
	 * \param e is an edge in the graph copy.
	 */
	void delEdge(edge e) override;

	void setOriginalEmbedding() override;

protected:
	void edgeInserted(void* userData, edge original, edge copy) override;
};

/**
 * \brief Copies of graphs supporting edge splitting
 *
 * @ingroup graphs
 *
 * The class GraphCopy represents a copy of a graph and
 * maintains a mapping between the nodes and edges of the original
 * graph to the copy and vice versa.
 *
 * New nodes and edges can be added to the copy; the counterpart
 * of those nodes and edges is 0 indicating that there is no counterpart.
 * GraphCopy also support splitting of edges in such a way
 * that both edges resulting from the split are mapped to the same
 * original edge, and each edge of the original graph is mapped to a list
 * of edges. Furthermore, it is allowed to reverse edges in the graph copy.
 *
 * <h3>Do's and Dont's</h3>
 * Here is a short summary, what can be done with GraphCopy, and what should not
 * be done. The following operations are safely supported:
 *   - Splitting of edges such that an original edge is represented by a path
 *     in the copy (split(), unsplit()).
 *   - Reversing edges in the copy (Graph::reverseEdge(), Graph::reverseAllEdges()).
 *   - Reinsertion of original edges such that they are represented by a path
 *     in the copy (insertEdgePath(), insertEdgePathEmbedded(), removeEdgePath(),
 *     removeEdgePathEmbedded()).
 *   - Inserting and removing dummy edges in the copy which are not associated
 *     with edges in the original graph.
 *
 * The following operations are <b>not supported</b> and are thus dangerous:
 *   - Any modifications on the original graph, since the copy will not be
 *     notified.
 *   - Moving the source or target node of an edge in the copy to a different node.
 *   - Removing edges in the graph copy that belong to a path representing an
 *     original edge.
 *   - ... (better think first!)
 */
class OGDF_EXPORT GraphCopy : public GraphCopyBase {
protected:
	EdgeArray<ListIterator<edge>> m_eIterator; //!< The position of copy edge in the list.
	EdgeArray<List<edge>> m_eCopy; //!< The corresponding list of edges in the graph copy.

public:
	explicit GraphCopy() : GraphCopy(nullptr) { }

	explicit GraphCopy(const Graph& G) : GraphCopy(&G) { }

	explicit GraphCopy(const Graph* G) : GraphCopyBase() {
		if (G) {
			init(*G);
		}
	}

	GraphCopy(const GraphCopy& other) : GraphCopyBase() { *this = other; }

	GraphCopy& operator=(const GraphCopy& other);

	using GraphCopyBase::setOriginalGraph;

	//! Associates the graph copy with \p G, but does not create any nodes or edges.
	/**
	 * By using this method, the graph copy gets associated with \p G.
	 * This does not modify the existing nodes and edges, but they become dummies.
	 *
	 * The following code snippet shows a typical application of this functionality:
	 * \code
	 *   GraphCopy GC;
	 *   GC.setOriginalGraph(&G);
	 *
	 *   // compute connected components of G
	 *   Graph::CCsInfo ccs = Graph::CCsInfo(graph);
	 *
	 *   // create and use graph copy for each connected component separately
	 *   NodeArray<edge> nodeCopy(G);
	 *   EdgeArray<edge> edgeCopy(G);
	 *   for (int i = 0; i < ccs.numberOfCCs(); i++) {
	 *     nodeCopy.init();
	 *     edgeCopy.init();
	 *     GC.clear();
	 *     GC.insert(ccs, i, nodeCopy, edgeCopy);
	 *     ...
	 *   }
	 * \endcode
	 * @param G is the graph of which this graph copy shall be a copy.
	 */
	void setOriginalGraph(const Graph* G) override;

	//! Removes all nodes and edges from this copy but does not break the link with the original graph.
	void clear() override;

	/**
	 * @name Mapping between original graph and copy
	 */
	//! @{

	/**
	 * \brief Returns the list of edges coresponding to edge \p e.
	 * \param e is an edge in the original graph.
	 * \return the corresponding list of edges in the graph copy.
	 */
	const List<edge>& chain(edge e) const { return m_eCopy[e]; }

	using GraphCopyBase::copy;

	/**
	 * \brief Returns the first edge in the list of edges corresponding to edge \p e.
	 * @param e is an edge in the original graph.
	 * \return the first edge in the corresponding list of edges in
	 * the graph copy or nullptr if it does not exist.
	 */
	edge copy(edge e) const override { return m_eCopy[e].empty() ? nullptr : m_eCopy[e].front(); }

	/**
	 * \brief Returns the adjacency entry in the copy graph corresponding to \p adj.
	 *
	 * Note that this method does not pay attention to reversed edges.
	 * Given a source (target) adjacency entry, the first (last) source (target) adjacency entry of the
	 * copy chain is returned.
	 *
	 * @param adj is an adjacency entry in the copy graph.
	 * \return the corresponding adjacency entry in the original graph or nullptr if it does not exist.
	 */
	adjEntry copy(adjEntry adj) const override {
		edge e = adj->theEdge();
		if (e == nullptr) {
			return nullptr;
		} else if (adj->isSource()) {
			return m_eCopy[e].front()->adjSource();
		} else {
			return m_eCopy[e].back()->adjTarget();
		}
	}

	/**
	 * \brief Returns true iff edge \p e has been reversed.
	 * @param e is an edge in the original graph.
	 */
	bool isReversed(edge e) const { return e->source() != original(copy(e)->source()); }

	/**
	 * \brief Returns true iff \p e is reversed w.r.t. the original edge of \a e.
	 * This method should be used, if the copy edge is split and \p e is part of the chain of the original edge.
	 * This method assumes that the list of copy edges forms a chain
	 * \param e is an edge in the graphcopy
	 */
	bool isReversedCopyEdge(edge e) const;
	//@}

	/**
	 * @name Creation and deletion of nodes and edges
	 */
	//! @{

	/**
	 * \brief Removes edge e and clears the list of edges corresponding to \p e's original edge.
	 *
	 * \pre The list of edges corresponding to \p e's original edge contains only \a e.
	 * \param e is an edge in the graph copy.
	 */
	void delEdge(edge e) override;

	/**
	 * \brief Splits edge \p e. See Graph::split() for details.
	 * Both resulting edges have the same original edge.
	 * @param e is an edge in the graph copy.
	 */
	edge split(edge e) override;

	/**
	 * \brief Undoes a previous split operation.
	 * The two edges \p eIn and \p eOut are merged to a single edge \a eIn.
	 * \pre The vertex \a u that was created by the previous split operation has
	 *      exactly one incoming edge \p eIn and one outgoing edge \p eOut.
	 * @param eIn is an edge (*,\a u) in the graph copy.
	 * @param eOut is an edge (\a u,*) in the graph copy.
	 */
	void unsplit(edge eIn, edge eOut) override;

	//! Creates a new edge (\a v,\a w) with original edge \a eOrig.
	edge newEdge(edge eOrig);

	using Graph::newEdge;

	//! sets eOrig to be the corresponding original edge of eCopy and vice versa
	/**
	 * @param eOrig is the original edge
	 * @param eCopy is the edge copy
	 */
	void setEdge(edge eOrig, edge eCopy);

	void setOriginalEmbedding() override;

	//! Removes all crossing nodes which are actually only two "touching" edges.
	void removePseudoCrossings();

	//! Returns whether there are two edges in the GraphCopy that cross each
	//! other multiple times.
	bool hasSameEdgesCrossings() const;

	//! Returns whether the GraphCopy contains at least one crossing of two
	//! adjacent edges.
	bool hasAdjacentEdgesCrossings() const;

	//! Returns whether the GraphCopy contains crossings that will result in a
	//! non-simple drawing.
	/**
	 * This method will return true iff the GraphCopy contains
	 * a) a crossing of two adjacent edges (see hasAdjacentEdgesCrossings()), or
	 * b) two edges crossing each other multiple times (see hasSameEdgesCrossings()).
	 *
	 * @warning Crossings of an edge with itself are currently not detected.
	 */
	inline bool hasNonSimpleCrossings() const {
		return hasAdjacentEdgesCrossings() || hasSameEdgesCrossings();
	};

	/**
	 * Removes all non-simple cossings involving edges from \p edgesToCheck
	 * (see hasNonSimpleCrossings() for a definition of non-simple crossings).
	 *
	 * @warning Crossings of an edge with itself are currently not removed.
	 *
	 * @param edgesToCheck edges from which non-simple crossings should be
	 * removed.
	 * @param dualGraph points to the dual graph of *this. Can be nullptr if
	 * only the GraphCopy should be changed.
	 */
	void removeNonSimpleCrossings(SListPure<edge>& edgesToCheck,
			DynamicDualGraph* dualGraph = nullptr);

	/**
	 * Removes all non-simple cossings (see hasNonSimpleCrossings() for a
	 * definition of non-simple crossings).
	 *
	 * @warning Crossings of an edge with itself are currently not removed.
	 *
	 * @param dualGraph points to the dual graph of *this. Can be nullptr if
	 * only the GraphCopy should be changed.
	 */
	inline void removeNonSimpleCrossings(DynamicDualGraph* dualGraph = nullptr) {
		SListPure<edge> edgesToCheck;
		m_pGraph->allEdges(edgesToCheck);
		removeNonSimpleCrossings(edgesToCheck, dualGraph);
	};

	/**
	 * Removes all non-simple cossings involving edges incident to \p origNode
	 * (see hasNonSimpleCrossings() for a definition of non-simple crossings).
	 *
	 * @warning Crossings of an edge with itself are currently not removed.
	 *
	 * @param origNode original node from whose incident edges non-simple
	 * crossings are removed.
	 * @param dualGraph points to the dual graph of *this. Can be nullptr if
	 * only the GraphCopy should be changed.
	 */
	inline void removeNonSimpleCrossings(node origNode, DynamicDualGraph* dualGraph = nullptr) {
		SListPure<edge> edgesToCheck;
		for (adjEntry adj : origNode->adjEntries) {
			edgesToCheck.pushBack(adj->theEdge());
		}
		removeNonSimpleCrossings(edgesToCheck, dualGraph);
	}

	//! Re-inserts edge \p eOrig by "crossing" the edges in \p crossedEdges.
	/**
	 * Let \a v and \a w be the copies of the source and target nodes of \p eOrig.
	 * Each edge in \p crossedEdges is split creating a sequence
	 * \a u_1, ..., \a u_k of new nodes, and additional edges are inserted creating
	 * a path \a v, \a u_1, ..., \a u_k, \a w.
	 * @param eOrig is an edge in the original graph and becomes the original edge of
	 *        all edges in the path \a v, \a u_1, ..., \a u_k, \a w.
	 * @param crossedEdges are edges in the graph copy.
	 */
	void insertEdgePath(edge eOrig, const SList<adjEntry>& crossedEdges);

	//! Special version (for FixedEmbeddingUpwardEdgeInserter only).
	void insertEdgePath(node srcOrig, node tgtOrig, const SList<adjEntry>& crossedEdges);

	//! Removes the complete edge path for edge \p eOrig.
	/**
	 * \@param eOrig is an edge in the original graph.
	 */
	void removeEdgePath(edge eOrig);

	//! Inserts crossings between two copy edges.
	/**
	 * This method is used in TopologyModule.
	 *
	 * Let \p crossingEdge = (\a a, \a b) and \p crossedEdge = (\a v, \a w).
	 * Then \p crossedEdge is split creating two edges \p crossedEdge = (\a v, \a u)
	 * and (\a u, \a w), \p crossingEdge is removed and replaced by two new edges
	 * \a e1  = (\a a, \a u) and \a e2 = (\a u, \a b).
	 * Finally it sets \p crossingEdge to \a e2 and returns (\a u, \a w).
	 *
	 * @param crossingEdge is the edge that is replaced by two new edges.
	 *                     Note that this parameter will be modified to equal the newly created edge (\a u, \a b).
	 * @param crossedEdge is the edge that gets split.
	 * @param rightToLeft is used as follows: If set to true, \p crossingEdge will cross
	 *        \p crossedEdge from right to left, otherwise from left to right.
	 * @return the rear edge resulting from the split operation: (\a u, \a w)
	 */
	edge insertCrossing(edge& crossingEdge, edge crossedEdge, bool rightToLeft);
	//@}

	/**
	 * @name Combinatorial Embeddings
	 */
	//! @{

	//! Creates a new edge with original edge \p eOrig in an embedding \p E.
	/**
	 * Let \a w be the node whose adjacency list contains \a adjTgt. The original
	 * edge \p eOrig must connect the original nodes of \p v and \a w. If \p eOrig =
	 * (original(\p v),original(\a w)), then the created edge is (\p v,\a w), otherwise
	 * it is (\a w,\p v). The new edge \a e must split a face in \p E, such that \a e
	 * comes after \p adj in the adjacency list of \p v and at the end of the adjacency
	 * list of \p v.
	 *
	 * @param v is a node in the graph copy.
	 * @param adj is an adjacency entry in the graph copy.
	 * @param eOrig is an edge in the original graph.
	 * @param E is an embedding of the graph copy.
	 * @return the created edge.
	 */
	edge newEdge(node v, adjEntry adj, edge eOrig, CombinatorialEmbedding& E);

	//! Re-inserts edge \p eOrig by "crossing" the edges in \p crossedEdges in embedding \p E.
	/**
	 * Let \a v and \a w be the copies of the source and target nodes of \p eOrig,
	 * and let \a e_0, \a e_1, ..., \a e_k, \a e_{k+1} be the sequence of edges corresponding
	 * to the adjacency entries in \p crossedEdges. The first edge needs to be incident
	 * to \a v and the last to \a w; the edges \a e_1, ..., \a e_k are each split
	 * creating a sequence \a u_1, ..., \a u_k of new nodes, and additional edges
	 * are inserted creating a path \a v, \a u_1, ..., \a u_k, \a w.
	 *
	 * The following figure illustrates, which adjacency entries need to be in the list
	 * \p crossedEdges. It inserts an edge connecting \a v and \a w by passing through
	 * the faces \a f_0, \a f_1, \a f_2; in this case, the list \p crossedEdges must contain
	 * the adjacency entries \a adj_0, ..., \a adj_3 (in this order).
	 * \image html insertEdgePathEmbedded.png
	 *
	 * @param eOrig is an edge in the original graph and becomes the original edge of
	 *        all edges in the path \a v, \a u_1, ..., \a u_k, \a w.
	 * @param E is an embedding of the graph copy.
	 * @param crossedEdges are a list of adjacency entries in the graph copy.
	 */
	void insertEdgePathEmbedded(edge eOrig, CombinatorialEmbedding& E,
			const SList<adjEntry>& crossedEdges);

	void insertEdgePathEmbedded(edge eOrig, CombinatorialEmbedding& E, DynamicDualGraph& dual,
			const SList<adjEntry>& crossedEdges);

	//! Removes the complete edge path for edge \p eOrig while preserving the embedding.
	/**
	 * If an endpoint of \p eOrig has degree 1, the node is also deleted (since
	 * removing the edge path would otherwise isolated the node, resulting in a
	 * broken embedding).
	 *
	 * @param E is an embedding of the graph copy.
	 * @param eOrig is an edge in the original graph.
	 * @param newFaces is assigned the set of new faces resulting from joining faces
	 *        when removing edges.
	 */
	void removeEdgePathEmbedded(CombinatorialEmbedding& E, edge eOrig, FaceSet<false>& newFaces);

	void removeEdgePathEmbedded(CombinatorialEmbedding& E, DynamicDualGraph& dual, edge eOrig);

	//! @}
	/**
	 * @name Miscellaneous
	 */
	//! @{

#ifdef OGDF_DEBUG
	//! Asserts that this copy is consistent.
	void consistencyCheck() const;
#endif

	//! @}

protected:
	void edgeInserted(void* userData, edge original, edge copy) override;


	void removeUnnecessaryCrossing(adjEntry adjA1, adjEntry adjA2, adjEntry adjB1, adjEntry adjB2);

	/**
	 * Removes the pseudo crossing that adj belongs to.
	 * In comparison to removeUnnecessaryCrossing(adjEntry, adjEntry, adjEntry,
	 * adjEntry), this method allows to change a dual graph as well.
	 *
	 * \pre adj->theNode() is a crossing with two incoming and two outgoing
	 * edges. adj and its successor must be part of the same original edge; the
	 * same holds for the next two successors respectively.
	 */
	void removeUnnecessaryCrossing(adjEntry adj, DynamicDualGraph* dualGraph);

	/**
	 * Removes the crossing of the two adjacent edges adj1->theEdge() and
	 * adj2->theEdge().
	 *
	 * \pre adj1 and adj2 are successive adjEntries of the same node, pointing
	 * towards the common node of both of their original edges.
	 */
	void removeAdjacentEdgesCrossing(adjEntry adj1, adjEntry adj2, DynamicDualGraph* dualGraph);

	/**
	 * Removes the two crossings given by the adjEntries, assuming that both
	 * crossings stem from the same two edges.
	 *
	 * \pre adjFirstCrossing1 and adjFirstCrossing2 as well as
	 * adjSecondCrossing1 and adjSecondCrossing2 are successive adjEntries of
	 * the same node respectively, such that the former point towards the latter
	 * and vice versa.
	 */
	void removeSameEdgesCrossing(adjEntry adjFirstCrossing1, adjEntry adjFirstCrossing2,
			adjEntry adjSecondCrossing1, adjEntry adjSecondCrossing2, DynamicDualGraph* dualGraph);

	/**
	 * Swaps the original edges from adjCopy1 up to the common node of
	 * {adjCopy1, adjCopy2} with the original edges from adjCopy2 up to the same
	 * common node. Can be used to fix crossings of adjacent edges.
	 *
	 * \pre Both adjCopy1 and adjCopy2 must point towards a common original node
	 * at the end of their chains.
	 */
	void swapOriginalEdgesAtCrossing(adjEntry adjCopy1, adjEntry adjCopy2,
			DynamicDualGraph* dual = nullptr);

	/**
	 * Swaps the original edges from adjFirstCrossing1 up to
	 * adjSecondCrossing1->theNode() with the original edges from
	 * adjFirstCrossing2 up to adjSecondCrossing2->theNode().
	 * Can be used to fix multiple crossings of the same two edges.
	 */
	void swapOriginalEdgesBetweenCrossings(adjEntry adjFirstCrossing1, adjEntry adjFirstCrossing2,
			adjEntry adjSecondCrossing1, adjEntry adjSecondCrossing2,
			DynamicDualGraph* dual = nullptr);

	/**
	 * Sets the original edges from adjCopy1 up to vCopy to eOrig2, and the
	 * original edges from adjCopy2 up to vCopy to eOrig1.
	 */
	void setOriginalEdgeAlongCrossings(adjEntry adjCopy1, adjEntry adjCopy2, node vCopy,
			edge eOrig1, edge eOrig2);
};

}
