/** \file
 * \brief Declaration of CombinatorialEmbedding and face.
 *
 * Enriches graph by the notion of faces
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

#pragma once

#include <ogdf/basic/AdjEntryArray.h>

namespace ogdf {

typedef FaceElement *face;

// Definition of iterator and container types for adjacency entries in a face
// These declarations are just internal representations
namespace internal {

//! Forward iterator for adjacency entries in a face
class FaceAdjIterator {

	adjEntry m_adj;
	adjEntry m_adjFirst;

public:
	FaceAdjIterator() : m_adj(nullptr), m_adjFirst(nullptr) { }
	FaceAdjIterator(adjEntry adj) : m_adj(adj), m_adjFirst(adj) { }
	FaceAdjIterator(adjEntry adjFirst, adjEntry adj) : m_adj(adj), m_adjFirst(adjFirst) { }

	bool operator==(const FaceAdjIterator &other) const {
		return m_adj == other.m_adj;
	}

	bool operator!=(const FaceAdjIterator &other) const {
		return m_adj != other.m_adj;
	}

	FaceAdjIterator &operator=(FaceAdjIterator &other) {
		m_adjFirst = other.m_adjFirst;
		m_adj = other.m_adj;
		return *this;
	}

	adjEntry operator*() const { return m_adj; }

	FaceAdjIterator &operator++() {
		m_adj = m_adj->faceCycleSucc();
		if (m_adj == m_adjFirst)
			m_adj = nullptr;
		return *this;
	}
};


//! Container for the adjacency entries in a face.
/**
 * The entries are not stored explicitly (in a list), but implicitly by the cyclic ordering of the adjacency lists
 * in the underlying graph and by storing the first adjacency entry in the face.
 */
class FaceAdjContainer {

	friend class ogdf::FaceElement;
	friend class ogdf::ConstCombinatorialEmbedding;
	friend class ogdf::CombinatorialEmbedding;

	adjEntry m_adjFirst;

	FaceAdjContainer() : m_adjFirst(nullptr) { }
	FaceAdjContainer(adjEntry adjFirst) : m_adjFirst(adjFirst) { }

public:
	typedef FaceAdjIterator iterator;

	iterator begin() const { return iterator(m_adjFirst); }
	iterator end() const { return iterator(); }
};

} // end namespace internal


/**
 * \brief Faces in a combinatorial embedding.
 */
class OGDF_EXPORT FaceElement : private internal::GraphElement
{
	friend class ConstCombinatorialEmbedding;
	friend class CombinatorialEmbedding;
	friend class internal::GraphList<FaceElement>;

	//adjEntry m_adjFirst; //!< The first adjacency element in the face.
	int m_id;   //!< The index of the face.
	int m_size; //!< The size of the face.

#ifdef OGDF_DEBUG
	const ConstCombinatorialEmbedding *m_pEmbedding;
#endif

	// constructor
#ifdef OGDF_DEBUG
	FaceElement(const ConstCombinatorialEmbedding *pEmbedding,
		adjEntry adjFirst,
		int id) :
		m_id(id), m_size(0), m_pEmbedding(pEmbedding), entries(adjFirst) { }
#else
	//! Creates a face with given first adjacency element \a adjFirst and face index \a id.
	FaceElement(adjEntry adjFirst, int id) :
		m_id(id), m_size(0), entries(adjFirst) { }
#endif

public:
	//! Container maintaining the adjacency entries in the face.
	internal::FaceAdjContainer entries;

	//! Returns the index of the face.
	int index() const { return m_id; }

	//! Returns the first adjacency element in the face.
	adjEntry firstAdj() const { return entries.m_adjFirst; }

	//! Returns the size of the face, i.e., the number of edges in the face.
	int size() const { return m_size; }

	//! Returns the successor in the list of all faces.
	face succ() const { return static_cast<face>(m_next); }

	//! Returns the predecessor in the list of all faces.
	face pred() const { return static_cast<face>(m_prev); }

	//! Returns the successor of \a adj in the list of all adjacency elements in the face.
	adjEntry nextFaceEdge(adjEntry adj) const {
		adj = adj->faceCycleSucc();
		return (adj != entries.m_adjFirst) ? adj : nullptr;
	}

#ifdef OGDF_DEBUG
	const ConstCombinatorialEmbedding *embeddingOf() const { return m_pEmbedding; }
#endif

	//! Standard Comparer
	static int compare(const FaceElement& x, const FaceElement& y) { return x.m_id - y.m_id; }
	OGDF_AUGMENT_COMPARER(FaceElement)

		OGDF_NEW_DELETE
}; // class FaceElement


class FaceArrayBase;
template<class T>class FaceArray;


/**
 * \brief Combinatorial embeddings of planar graphs.
 *
 * @ingroup graphs
 *
 * Maintains a combinatorial embedding of an embedded graph, i.e., the set of
 * faces. A combinatorial embedding is defined by the (cyclic) order of the
 * adjacency entries around a vertex; more precisely, the adjacency list
 * gives the cyclic order of the adjacency entries in clockwise order.
 * Each adjacency entry \a adj is contained in exactly one face, the face
 * to the right of \a adj. The list of adjacency entries defining a face is given
 * in clockwise order for internal faces, and in counter-clockwise order for the
 * external face.
 *
 * <H3>Thread Safety</H3>
 * The class Graph allows shared access of threads to const methods only.
 * If one thread executes a non-const method, shared access is no longer thread-safe.
 *
 * \see CombinatorialEmbedding provides additional functionality for modifying
 *      the embedding.
 */
class OGDF_EXPORT ConstCombinatorialEmbedding
{
protected:
	const Graph *m_cpGraph; //!< The associated graph.

	int m_faceIdCount; //!< The index assigned to the next created face.
	int m_faceArrayTableSize; //!< The current table size of face arrays.

	AdjEntryArray<face> m_rightFace; //!< The face to which an adjacency entry belongs.
	face m_externalFace; //! The external face.

	mutable ListPure<FaceArrayBase*> m_regFaceArrays; //!< The registered face arrays.

#ifndef OGDF_MEMORY_POOL_NTS
	mutable std::mutex m_mutexRegArrays; //!< The critical section for protecting shared acces to register/unregister methods.
#endif

public:
	//! Provides a bidirectional iterator to a face in a combinatorial embedding.
	typedef internal::GraphIterator<face> face_iterator;

	//! The container containing all face objects.
	internal::GraphObjectContainer<FaceElement> faces;

	/** @{
	 * \brief Creates a combinatorial embedding associated with no graph.
	 */
	ConstCombinatorialEmbedding();

	/**
	 * \brief Creates a combinatorial embedding of graph \a G.
	 *
	 * \pre Graph \a G must be embedded, i.e., the adjacency lists of its nodes
	 *      must define an embedding.
	 */
	explicit ConstCombinatorialEmbedding(const Graph &G);


	//! Copy constructor.
	ConstCombinatorialEmbedding(const ConstCombinatorialEmbedding &C);

	//! Assignment operator.
	ConstCombinatorialEmbedding &operator=(const ConstCombinatorialEmbedding &C);

	//! Destructor
	virtual ~ConstCombinatorialEmbedding();

	//! Returns whether the embedding is associated with a graph.
	bool valid() const { return m_cpGraph != nullptr; }

	/** @} @{
	 * \brief Returns the associated graph of the combinatorial embedding.
	 *
	 * \pre The associated graph exists. See #valid().
	 */
	const Graph &getGraph() const {
		OGDF_ASSERT(valid());
		return *m_cpGraph;
	}

	//! Returns associated graph
	operator const Graph &() const { return getGraph(); }

	/** @} @{
	 * \brief Returns the first face in the list of all faces.
	 */
	face firstFace() const { return faces.head(); }

	//! Returns the last face in the list of all faces.
	face lastFace() const { return faces.tail(); }

	//! Returns the number of faces.
	int numberOfFaces() const { return faces.size(); }

	/** @} @{
	 * \brief Returns the face to the right of \a adj, i.e., the face containing \a adj.
	 * @param adj is an adjecency element in the associated graph.
	 */
	face rightFace(adjEntry adj) const { return m_rightFace[adj]; }

	/**
	 * \brief Returns the face to the left of \a adj, i.e., the face containing the twin of \a adj.
	 * @param adj is an adjacency element in the associated graph.
	 */
	face leftFace(adjEntry adj) const { return m_rightFace[adj->twin()]; }

	/** @} @{
	 * \brief Returns the largest used face index.
	 */
	int maxFaceIndex() const { return m_faceIdCount-1; }

	//! Returns the table size of face arrays associated with this embedding.
	int faceArrayTableSize() const { return m_faceArrayTableSize; }

	/** @} @{
	 * Returns a random face.
	 * \c nullptr is returned if no feasible face exists.
	 *
	 * @see chooseIteratorFrom
	 */
	face chooseFace(std::function<bool(face)> includeFace = [](face) { return true; }, bool isFastTest = true) const;

	//! Returns a face of maximal size.
	face maximalFace() const;

	/** @} @{
	 * \brief Returns the external face.
	 */
	face externalFace() const {
		return m_externalFace;
	}

	/**
	 * \brief Sets the external face to \a f.
	 * @param f is a face in this embedding.
	 */
	void setExternalFace(face f) {
		OGDF_ASSERT(f->embeddingOf() == this);
		m_externalFace = f;
	}

	bool isBridge(edge e) const {
		return m_rightFace[e->adjSource()] == m_rightFace[e->adjTarget()];
	}

	/** @} @{
	 * \brief Initializes the embedding for graph \a G.
	 *
	 * \pre Graph \a G must be embedded, i.e., the adjacency lists of its nodes
	 *      must define an embedding.
	 */
	void init(const Graph &G);

	void init();

	//! Computes the list of faces.
	void computeFaces();


	/** @} @{
	 * \brief Checks the consistency of the data structure.
	 */
	bool consistencyCheck();


	/** @} @{
	 * \brief Registers the face array \a pFaceArray.
	 *
	 * This method is only used by face arrays.
	 */
	ListIterator<FaceArrayBase*> registerArray(FaceArrayBase *pFaceArray) const;

	/**
	 * \brief Unregisters the face array identified by \a it.
	 *
	 * This method is only used by face arrays.
	 */
	void unregisterArray(ListIterator<FaceArrayBase*> it) const;

	//! Move the registration \a it of a node array to \a pFaceArray (used with move semantics for face arrays).
	void moveRegisterArray(ListIterator<FaceArrayBase*> it, FaceArrayBase *pFaceArray) const;

	/** @} */

protected:
	//! Create a new face.
	face createFaceElement(adjEntry adjFirst);

	//! Reinitialize associated face arrays.
	void reinitArrays();

}; // class ConstCombinatorialEmbedding



/**
 * \brief Combinatorial embeddings of planar graphs with modification functionality.
 *
 * @ingroup graphs
 *
 * Maintains a combinatorial embedding of an embedded graph, i.e., the set of
 * faces, and provides method for modifying the embedding, e.g., by inserting edges.
 *
 * <H3>Thread Safety</H3>
 * The class Graph allows shared access of threads to const methods only.
 * If one thread executes a non-const method, shared access is no longer thread-safe.
 */
class OGDF_EXPORT CombinatorialEmbedding : public ConstCombinatorialEmbedding
{
	friend class GraphCopy;

	Graph *m_pGraph; //!< The associated graph.

	// the following methods are explicitly deleted
	// It is not clear which meaning copying of a comb. embedding should
	// have since we only store a pointer to the topology (Graph)
	CombinatorialEmbedding(const CombinatorialEmbedding &) = delete;
	CombinatorialEmbedding &operator=(const CombinatorialEmbedding &) = delete;

public:
	/** @{
	 * \brief Creates a combinatorial embedding associated with no graph.
	 */
	CombinatorialEmbedding() : ConstCombinatorialEmbedding() {
		m_pGraph = nullptr;
	}

	/**
	 * \brief Creates a combinatorial embedding of graph \a G.
	 *
	 * \pre Graph \a G must be embedded, i.e., the adjacency lists of its nodes
	 *      must define an embedding.
	 */
	explicit CombinatorialEmbedding(Graph &G) : ConstCombinatorialEmbedding(G) {
		m_pGraph = &G;
	}

	//@}
	/**
	 * @name Access to the associated graph
	 */
	//@{

	//! Returns the associated graph.
	const Graph &getGraph() const
	{
		OGDF_ASSERT(valid());
		return *m_cpGraph;
	}

	Graph &getGraph()
	{
		OGDF_ASSERT(valid());
		return *m_pGraph;
	}

	operator const Graph &() const { return getGraph(); }

	operator Graph &() { return getGraph(); }


	//@}
	/**
	 * @name Initialization
	 */
	//@{

	/**
	 * \brief Initializes the embedding for graph \a G.
	 *
	 * \pre Graph \a G must be embedded, i.e., the adjacency lists of its nodes
	 *      must define an embedding.
	 */
	void init(Graph &G) {
		ConstCombinatorialEmbedding::init(G);
		m_pGraph = &G;
	}

	/**
	 * \brief Removes all nodes, edges, and faces from the graph and the embedding.
	 */
	void clear();


	//@}
	/**
	 * @name Update of embedding
	 */
	//@{

	/**
	 * \brief Splits edge \a e=(\a v,\a w) into \a e=(\a v,\a u) and \a e'=(\a u,\a w) creating a new node \a u.
	 * @param e is the edge to be split; \a e is modified by the split.
	 * \return the edge \a e'.
	 */
	edge split(edge e);

	/**
	 * \brief Undoes a split operation.
	 * @param eIn is the edge (\a v,\a u).
	 * @param eOut is the edge (\a u,\a w).
	 */
	void unsplit(edge eIn, edge eOut);

	/**
	 * \brief Splits a node while preserving the order of adjacency entries.
	 *
	 * This method splits a node \a v into two nodes \a vl and \a vr. Node
	 * \a vl receives all adjacent edges of \a v from \a adjStartLeft until
	 * the edge preceding \a adjStartRight, and \a vr the remaining nodes
	 * (thus \a adjStartRight is the first edge that goes to \a vr). The
	 * order of adjacency entries is preserved. Additionally, a new edge
	 * (\a vl,\a vr) is created, such that this edge is inserted before
	 * \a adjStartLeft and \a adjStartRight in the the adjacency lists of
	 * \a vl and \a vr.
	 *
	 * Node \a v is modified to become node \a vl, and node \a vr is returned.
	 *
	 * @param adjStartLeft is the first entry that goes to the left node.
	 * @param adjStartRight is the first entry that goes to the right node.
	 * \return the newly created node.
	 */
	node splitNode(adjEntry adjStartLeft, adjEntry adjStartRight);

	/**
	 * \brief Contracts edge \a e.
	 * @param e is an edge is the associated graph.
	 * @return the node resulting from the contraction.
	 */
	node contract(edge e);

	/**
	 * Splits a face by inserting a new edge.
	 *
	 * Creates a new edge from the node of \c adjSrc to the one of \c adjTgt.
	 * Note that this can also be achieved by inserting an edge
	 * in the underlying graph directly and calling #computeFaces again.
	 * In contrast, this operation achieves constant running time.
	 *
	 * \pre \c adjSrc and \c adjTgt are distinct AdjEntries, belonging to the same face.
	 * \return The new edge.
	 */
	edge splitFace(adjEntry adjSrc, adjEntry adjTgt);

	/**
	 * Inserts a new edge similarly to #splitFace without having to call #computeFaces again.
	 *
	 * Creates a new edge from the degree 0 node \c v to the node of \c adjTgt.
	 * The face that \c adjTgt belongs to is split.
	 *
	 * \return The new edge.
	 */
	edge addEdgeToIsolatedNode(node v, adjEntry adjTgt);

	/**
	 * Inserts a new edge similarly to #splitFace without having to call #computeFaces again.
	 *
	 * Creates a new edge from the node of \c adjSrc to the degree 0 node \c v.
	 * The face that \c adjSrc belongs to is split.
	 *
	 * \return The new edge.
	 */
	edge addEdgeToIsolatedNode(adjEntry adjSrc, node v);

	/**
	 * \brief Removes edge e and joins the two faces adjacent to \a e.
	 * @param e is an edge in the associated graph.
	 * \return the resulting (joined) face.
	 */
	face joinFaces(edge e);

	//! Reverses edges \a e and updates embedding.
	void reverseEdge(edge e);

	//! Moves a bridge in the graph.
	void moveBridge(adjEntry adjBridge, adjEntry adjBefore);

	//! Removes degree-1 node \a v.
	void removeDeg1(node v);

	//! Update face information after inserting a merger in a copy graph.
	void updateMerger(edge e, face fRight, face fLeft);


	/** @} */

protected:
	/**
	 * \brief Joins the two faces adjacent to \a e but does not remove edge \e e.
	 * @param e is an edge in the associated graph.
	 * \return the resulting (joined) face.
	 */
	face joinFacesPure(edge e);

private:
	/**
	 * Inserts a new edge similarly to #splitFace without having to call #computeFaces again.
	 *
	 * \param adj The adjacency entry belonging to the face that we want to insert the new edge into
	 * \param v The degree 0 node.
	 * \param adjSrc whether v will be the target node.
	 * \return The new edge.
	 */
	edge addEdgeToIsolatedNode(adjEntry adj, node v, bool adjSrc);
}; // class CombinatorialEmbedding

} // end namespace ogdf
