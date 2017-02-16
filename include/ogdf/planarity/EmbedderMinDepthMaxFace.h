/** \file
 * \brief Computes an embedding of a graph with minimum depth and
 * maximum external face.
 *
 * \author Thorsten Kerkhof
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

#include <ogdf/module/EmbedderModule.h>
#include <ogdf/decomposition/BCTree.h>
#include <ogdf/planarity/embedder/MDMFLengthAttribute.h>

namespace ogdf {

//! Planar graph embedding with minimum block-nesting depth and maximum external face.
/**
 * @ingroup ga-planembed
 *
 * See the paper "Graph Embedding with Minimum Depth and Maximum External Face"
 * by C. Gutwenger and P. Mutzel (2004) for details.
 */
class OGDF_EXPORT EmbedderMinDepthMaxFace : public EmbedderModule
{
public:
	//constructor:
	EmbedderMinDepthMaxFace() { }

	/**
	 * \brief Call embedder algorithm.
	 * \param G is the original graph. Its adjacency list has to be  changed by the embedder.
	 * \param adjExternal is an adjacency entry on the external face and has to be set by the embedder.
	 */
	virtual void doCall(Graph& G, adjEntry& adjExternal) override;

private:
	using MDMFLengthAttribute = embedder::MDMFLengthAttribute;

	/**
	 * \brief Bottom-up-traversal of bcTree computing the values \a m_{cT, bT}
	 * for all edges \a (cT, bT) in the BC-tree. The length of each vertex
	 * \a v != \a c in \p bT is set to 1 if \a v in M_{bT} and to 0 otherwise.
	 *
	 * \param bT is a block vertex in the BC-tree.
	 * \param cH is a vertex in the original graph \a G.
	 * \return Minimum depth of an embedding of \p bT with \p cH on the external
	 *    face.
	 */
	int md_bottomUpTraversal(const node& bT, const node& cH);

	/**
	 * \brief Top-down-traversal of BC-tree. The minimum depth of the BC-tree-node
	 * bT is calculated and before calling the function recursively for all
	 * children of bT in the BC-tree, the nodeLength of the cut-vertex which bT
	 * and the child have in common is computed. The length of each node is set to
	 * 1 if it is in M_B and 0 otherwise, except for |M_B| = 1, than it is set to
	 * 1 if it is in M2 with m2 = max{m_B(v) : v in V_B, v != c} and
	 * M2 = {c in V_B \ {v} | m_B(c) = m2}.
	 *
	 * \param bT is a block vertex in the BC-tree.
	 */
	void md_topDownTraversal(const node& bT);

	/**
	 * \brief Bottom up traversal of BC-tree.
	 *
	 * \param bT is the BC-tree node treated in this function call.
	 * \param cH is the block node which is related to the cut vertex which is
	 *   parent of bT in BC-tree.
	 */
	int mf_constraintMaxFace(const node& bT, const node& cH);

	/**
	 * \brief Top down traversal of BC-tree.
	 *
	 * \param bT is the tree node treated in this function call.
	 * \param bT_opt is assigned a block node in BC-tree which contains a face which
	 *   cann be expanded to a maximum face.
	 * \param ell_opt is the size of a maximum face.
	 */
	void mf_maximumFaceRec(const node& bT, node& bT_opt, int& ell_opt);

	/**
	 * \brief Computes the adjacency list for all nodes in a block and calls
	 * recursively the function for all blocks incident to nodes in bT.
	 *
	 * \param bT is the tree node treated in this function call.
	 */
	void embedBlock(const node& bT);

	/**
	 * \brief Computes the adjacency list for all nodes in a block and calls
	 * recursively the function for all blocks incident to nodes in bT.
	 *
	 * \param bT is the tree node treated in this function call.
	 * \param cT is the parent cut vertex node of bT in the BC-tree. cT is 0 if bT
	 *   is the root block.
	 * \param after is the adjacency entry of the cut vertex, after which bT has to
	 *   be inserted.
	 */
	void embedBlock(const node& bT, const node& cT, ListIterator<adjEntry>& after);

private:
	/** the BC-tree of G */
	BCTree* pBCTree;

	/** an adjacency entry on the external face */
	adjEntry* pAdjExternal;

	/** saving for each node in the block graph its length */
	NodeArray<int> md_nodeLength;

	/** an array containing the minimum depth of each block */
	NodeArray<int> md_minDepth;

	/** an array saving the length for each edge in the BC-tree */
	EdgeArray<int> md_m_cB;

	/** M_B = {cH in B | m_B(cH) = m_B} with m_B = max{m_B(c) : c in B}
	 *  and m_B(c) = max( {0} cup {m_{c, B'} | c in B', B' != B} ). */
	NodeArray< List<node> > md_M_B;

	/** M2 is empty, if |M_B| != 1, otherwise M_B = {cH}
	 *  M2 = {cH' in V_B \ {v} | m_B(cH') = m2} with
	 *  m2 = max{m_B(vH) : vH in V_B, vH != cH}. */
	NodeArray< List<node> > md_M2;

	/** is saving for each node of the block graph its length */
	NodeArray<int> mf_nodeLength;

	/** is saving for each node of the block graph its cstrLength */
	NodeArray<int> mf_cstrLength;

	/** an array containing the maximum face size of each block */
	NodeArray<int> mf_maxFaceSize;

	/** is saving for each node of the block graph its length */
	NodeArray<MDMFLengthAttribute> mdmf_nodeLength;

	/** is saving for each edge of the block graph its length */
	EdgeArray<MDMFLengthAttribute> mdmf_edgeLength;

	/** saves for every node of G the new adjacency list */
	NodeArray< List<adjEntry> > newOrder;

	/** treeNodeTreated saves for all block nodes in the
	 *  BC-tree if it has already been treated or not. */
	NodeArray<bool> treeNodeTreated;
};

} // end namespace ogdf
