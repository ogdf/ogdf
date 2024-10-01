/** \file
 * \brief Multiple GraphCopies that contain different, overlapping parts of the same original graph. TODO should be moved to a central location.
 *
 * \author Simon D. Fink <ogdf@niko.fink.bayern>
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

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/Observer.h>
#include <ogdf/basic/basic.h>
#include <ogdf/cluster/sync_plan/basic/GraphMultiArray.h>

namespace ogdf {

class OverlappingGraphCopies;

// room for improvement: common interface with GraphCopyBase, add auto-linking ::insert method
//! Version of GraphCopySimple that may efficiently share some overlap with other instances of the same original Graph via its OverlappingGraphCopies manager.
class OGDF_EXPORT OverlappingGraphCopy : public Graph {
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

//! The manager class for multiple OverlappingGraphCopy instances of the same graph.
/**
 * This is similar to using multiple GraphCopySimple instances for a single graph, but is more efficient
 * storage-wise and easily allows enumerating in which copies a node or edge occurs.
 */
class OGDF_EXPORT OverlappingGraphCopies {
	friend class OverlappingGraphCopy;
	using NA = NodeMultiArray<const OverlappingGraphCopy*, node>;
	using EA = EdgeMultiArray<const OverlappingGraphCopy*, edge>;

	const Graph* m_G;
	NA m_node_copies;
	EA m_edge_copies;

public:
	explicit OverlappingGraphCopies(const Graph& G)
		: m_G(&G), m_node_copies(G), m_edge_copies(G) { }

	OGDF_NO_COPY(OverlappingGraphCopies)

	OGDF_NO_MOVE(OverlappingGraphCopies)

	const NA::EntryType& copies(node n) const { return m_node_copies.get_all(n); }

	const EA::EntryType& copies(edge e) const { return m_edge_copies.get_all(e); }

	const Graph* constGraph() const { return m_G; }
};

}
