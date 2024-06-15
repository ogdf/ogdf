#pragma once

#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/Graph_d.h>
#include <ogdf/basic/NodeArray.h>
#include <ogdf/cluster/sync_plan/basic/RegisteredMultiArray.h>

using namespace ogdf;


template<typename T>
using NA = NodeArray<T, true>;
template<typename T>
using EA = EdgeArray<T, true>;

template<typename Key2, typename Value>
using NodeMultiArray = RegisteredMultiArray<node, Key2, Value, NA>;
template<typename Key2, typename Value>
using EdgeMultiArray = RegisteredMultiArray<edge, Key2, Value, EA>;

class OverlappingGraphCopies;

class OverlappingGraphCopy : public Graph { // TODO common interface with GraphCopyBase
	friend class OverlappingGraphCopies;

	OverlappingGraphCopies* m_pOGC; //!< The master instance.
	NodeArray<node> m_vOrig; //!< The corresponding node in the original graph.
	EdgeArray<edge> m_eOrig; //!< The corresponding edge in the original graph.

	void unmap();

public:
	explicit OverlappingGraphCopy(OverlappingGraphCopies& mPOgc)
		: m_pOGC(&mPOgc), m_vOrig(*this, nullptr), m_eOrig(*this, nullptr) { }

	~OverlappingGraphCopy() override { unmap(); }

	OGDF_NO_MOVE(OverlappingGraphCopy)

	OGDF_NO_COPY(OverlappingGraphCopy)

	//! Returns a reference to the master instance.
	const OverlappingGraphCopies& master() const { return *m_pOGC; }

	//! Returns a reference to the original graph.
	const Graph& original() const;

	/**
     * \brief Returns the node in the original graph corresponding to \p v.
     * @param v is a node in the graph copy.
     * \return the corresponding node in the original graph, or 0 if no
     *         such node exists.
     */
	node original(node v) const {
		node ov = m_vOrig[v];
		OGDF_ASSERT(ov == nullptr || ov->graphOf() == &original());
		return ov;
	}

	/**
     * \brief Returns the edge in the original graph corresponding to \p e.
     * @param e is an edge in the graph copy.
     * \return the corresponding edge in the original graph, or 0 if no
     *         such edge exists.
     */
	edge original(edge e) const {
		edge oe = m_eOrig[e];
		OGDF_ASSERT(oe == nullptr || oe->graphOf() == &original());
		return oe;
	}

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
		edge f = original(adj->theEdge());
		return adj->isSource() ? f->adjSource() : f->adjTarget();
	}

	/**
     * \brief Returns the node in the graph copy corresponding to \p v.
     * @param v is a node in the original graph.
     * \return the corresponding node in the graph copy,
     * or \c nullptr if it doesn't exists.
     */
	node copy(node v) const;

	/**
     * \brief Returns the edge in the graph copy corresponding to \p e.
     * @param e is an edge in the original graph.
     * \return the corresponding edge in the graph copy,
     * or \c nullptr if it doesn't exists.
     */
	edge copy(edge e) const;

	/**
     * Returns the adjacency entry in the graph copy corresponding to \p adj.
     *
     * Note that this method does not pay attention to reversed edges.
     * Given a source (target) adjacency entry, the source (target) adjacency entry of the
     * copy edge is returned.
     *
     * @param adj is an adjacency entry in the original graph.
     * \return the corresponding adjacency entry in the graph copy,
     * or \c nullptr if it doesn't exists.
     */
	adjEntry copy(adjEntry adj) const {
		edge f = copy(adj->theEdge());
		if (f == nullptr) {
			return nullptr;
		}
		return adj->isSource() ? f->adjSource() : f->adjTarget();
	}

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
     * \brief Creates a new node in the graph copy with original node \p vOrig.
     * \warning You have to make sure that the original node makes sense, in
     *   particular that \p vOrig is not the original node of another node in the copy.
     */
	node newNode(node vOrig);

	using Graph::newNode;

	/**
     * \brief Creates a new edge in the graph copy with original edge \p eOrig.
     * \warning You have to make sure that the original edge makes sense, in
     *   particular that \p eOrig is not the original edge of another edge in the copy.
     */
	edge newEdge(edge eOrig);

	using Graph::newEdge;

	/**
     * \brief Removes edge \p e.
     *
     * \param e is an edge in the graph copy.
     */
	void delEdge(edge e) override;

	/**
     * \brief Removes node \p v.
     *
     * \param v is a node in the graph copy.
     */
	void delNode(node v) override;

	void clear() override {
		unmap();
		Graph::clear();
	}

	void breakLinkForMasterDeconstruction() {
		m_pOGC = nullptr;
		m_vOrig.init();
		m_eOrig.init();
	}
};

class OverlappingGraphCopies {
	friend class OverlappingGraphCopy;

	const Graph* G;
	mutable NodeMultiArray<const OverlappingGraphCopy*, node> node_copies;
	mutable EdgeMultiArray<const OverlappingGraphCopy*, edge> edge_copies;

public:
	explicit OverlappingGraphCopies(const Graph& G) : G(&G), node_copies(G), edge_copies(G) { }

	OGDF_NO_COPY(OverlappingGraphCopies)

	OGDF_NO_MOVE(OverlappingGraphCopies)
};
