/** \file
 * \brief Declaration of class StaticSkeleton.
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
#include <ogdf/decomposition/Skeleton.h>

namespace ogdf {

class StaticSPQRTree;
class SPQRTree;

//! %Skeleton graphs of nodes in a static SPQR-tree.
/**
 * @ingroup dc-helper
 *
 * The class StaticSkeleton maintains the skeleton \a S of a node \a vT in a static SPQR-tree
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
class OGDF_EXPORT StaticSkeleton : public Skeleton {
	friend class StaticSPQRTree;

public:
	// constructor

	//! Creates a skeleton \a S with owner tree \a T and corresponding node \p vT.
	/**
	 * \pre \p vT is a node in \a T
	 *
	 * \remarks Skeletons are created by the constructor of SPQRTree
	 *          and can be accessed with the skeleton() function of SPQRTree.
	 */
	StaticSkeleton(const StaticSPQRTree* T, node vT);

	// destructor
	~StaticSkeleton() { }

	//! Returns the owner tree \a T.
	const SPQRTree& owner() const override;

	//! Returns the vertex in the original graph \a G that corresponds to \p v.
	/**
	 * \pre \p v is a node in \a M.
	 */
	node original(node v) const override { return m_orig[v]; }

	//! Returns true iff \p e is a virtual edge.
	/**
	 * \pre \p e is an edge in \a M
	 */
	bool isVirtual(edge e) const override { return m_real[e] == nullptr; }

	//! Returns the real edge that corresponds to skeleton edge \p e
	/**
	 * If \p e is virtual edge, then 0 is returned.
	 * \pre \p e is an edge in \a M
	 */
	edge realEdge(edge e) const override { return m_real[e]; }

	//! Returns the twin edge of skeleton edge \p e.
	/**
	 * If \p e is not a virtual edge, then 0 is returned.
	 * \pre \p e is an edge in \a M
	 */
	edge twinEdge(edge e) const override;

	//! Returns the tree node in T containing the twin edge of skeleton edge \p e.
	/**
	 * If \p e is not a virtual edge, then 0 is returned.
	 * \pre \p e is an edge in \a M
	 */
	node twinTreeNode(edge e) const override;

	//! Returns the tree edge which is associated with skeleton edge \p e.
	/**
	 * If \p e is not a virtual edge, then 0 is retuned.
	 * \pre \p e is an edge in \a M
	 */
	edge treeEdge(edge e) const { return m_treeEdge[e]; }

	OGDF_MALLOC_NEW_DELETE

protected:
	const StaticSPQRTree* m_owner; //!< owner tree
	NodeArray<node> m_orig; //!< corresp. original node
	EdgeArray<edge> m_real; //!< corresp. original edge (0 if virtual)
	EdgeArray<edge> m_treeEdge; //!< corresp. tree edge (0 if real)
};

}
