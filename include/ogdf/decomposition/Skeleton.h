/** \file
 * \brief Declaration of class Skeleton.
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

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/memory.h>

namespace ogdf {

class SPQRTree;

//! %Skeleton graphs of nodes in an SPQR-tree.
/**
 * @ingroup dc-helper
 *
 * The class Skeleton maintains the skeleton \a S of a node \a vT in an SPQR-tree
 * \a T. We call \a T the owner tree of \a S and \a vT the corresponding tree node. Let
 * \a G be the original graph of \a T.
 *
 * \a S consists of an undirected multi-graph \a M. If the owner tree of \a S is
 * a PlanarSPQRTree, then \a M represents a combinatorial embedding.
 * The vertices in \a M correspond to vertices in \a G. The edges in \a M are of
 * two types: Real edges correspond to edges in \a G and virtual edges correspond
 * to tree-edges in \a T having \a vT as an endpoint.
 *
 * Let \a e in \a M be a virtual edge and \a eT be the corresponding tree-edge.
 * Then there exists exactly one edge \a e' in another skeleton \a S' of \a T that
 * corresponds to \a eT as well. We call \a e' the twin edge of \a e.
 */
class OGDF_EXPORT Skeleton {
public:
	// constructor

	//! Creates a skeleton \a S with owner tree \a T and corresponding node \p vT.
	/**
	 * \pre \p vT is a node in \a T
	 *
	 * \remarks Skeletons are created by the constructor of SPQRTree
	 *          and can be accessed with the skeleton() function of SPQRTree.
	 */
	explicit Skeleton(node vT) : m_treeNode(vT) { }

	// destructor
	virtual ~Skeleton() { }

	//! Returns the owner tree \a T.
	virtual const SPQRTree& owner() const = 0;

	//! Returns the corresponding node in the owner tree \a T to which \a S belongs.
	node treeNode() const { return m_treeNode; }

	//! Returns the reference edge of \a S in \a M.
	/**
	 * The reference edge is a real edge if \a S is the skeleton of the root node
	 * of \a T, and a virtual edge otherwise.
	 */
	edge referenceEdge() const { return m_referenceEdge; }

	//! Returns a reference to the skeleton graph \a M.
	const Graph& getGraph() const { return m_M; }

	//! Returns a reference to the skeleton graph \a M.
	Graph& getGraph() { return m_M; }

	//! Returns the vertex in the original graph \a G that corresponds to \p v.
	/**
	 * \pre \p v is a node in \a M.
	 */
	virtual node original(node v) const = 0;

	//! Returns true iff \p e is a virtual edge.
	/**
	 * \pre \p e is an edge in \a M
	 */
	virtual bool isVirtual(edge e) const = 0;

	//! Returns the real edge that corresponds to skeleton edge \p e
	/**
	 * If \p e is virtual edge, then 0 is returned.
	 * \pre \p e is an edge in \a M
	 */
	virtual edge realEdge(edge e) const = 0;

	//! Returns the twin edge of skeleton edge \p e.
	/**
	 * If \p e is not a virtual edge, then 0 is returned.
	 * \pre \p e is an edge in \a M
	 */
	virtual edge twinEdge(edge e) const = 0;

	//! Returns the tree node in T containing the twin edge of skeleton edge \p e.
	/**
	 * If \p e is not a virtual edge, then 0 is returned.
	 * \pre \p e is an edge in \a M
	 */
	virtual node twinTreeNode(edge e) const = 0;

	OGDF_NEW_DELETE

protected:
	node m_treeNode; //!< corresp. tree node in owner tree
	edge m_referenceEdge; //!< reference edge
	Graph m_M; //!< actual skeleton graph
};

}
