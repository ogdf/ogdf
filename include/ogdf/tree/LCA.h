/*
 * $Revision: 3550 $
 *
 * last checkin:
 *   $Author: beyer $
 *   $Date: 2013-06-07 14:16:24 +0200 (Fr, 07. Jun 2013) $
 ***************************************************************/

/** \file
 * \brief The Sparse Table Algorithm for the Least Common Ancestor
 *  problem as proposed by Bender and Farach-Colton.
 *
 * \author Matthias Woste
 * Porting to OGDF Graph class based trees by Stephan Beyer
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.txt in the root directory of the OGDF installation for details.
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
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/

#ifndef LCA_OGDF_H_
#define LCA_OGDF_H_

#include <ogdf/basic/Graph.h>

namespace ogdf {

/*!
 * \brief This class implements the <O(n log n), O(1)>-time "sparse table" LCA algorithm
 * by Bender and Farach-Colton.
 *
 * This implementation is based on:
 *
 * (M. Bender, M. Farach-Colton, The LCA problem revisited, LATIN '00, volume 1776 of LNCS,
 * pages 88-94, Springer, 2000)
 */
class LCA {
public:
	/*!
	 * \brief Builds the LCA data structure for an arborescence
	 * @param G a tree
	 * @param root root node of the tree
	 * \pre Each node in G is reachable from the root via a unique directed path, that is, G is an arborescence.
	 */
	LCA(const Graph &G, node root);

	/*!
	 * \brief Returns the LCA of two nodes. If the two nodes are the same, the node itself is defined to be the LCA.
	 * @param u First node
	 * @param v Second node
	 * @return The LCA of u and v
	 */
	node call(node u, node v) const;

	//! \brief Returns the level of a node. The level of the root is 0.
	int level(node v) const
	{
		return m_level[m_representative[v]];
	}

private:
	const node m_root; //!< the root of the tree
	const int m_n; //!< number of nodes in graph
	const int m_len; //!< length of the RMQ array (always 2 m_n - 1)
	const int m_rangeJ; //!< always floor(log(m_len)) (size of a row in the table)
	Array<node> m_euler; //!< Euler[i] is i-th node visited in Euler Tour
	NodeArray<int> m_representative; //!< Euler[Representative[v]] = v
	Array<int> m_level; //!< L[i] is distance of node E[i] from root
	Array<int> m_table; //!< preprocessed M[i,j] array

	/*!
	 * \brief Performs an Euler tour (actually a DFS with virtual back-edges) through the underlying tree
	 *  and fill Euler tour and Level arrays.
	 */
	void dfs(const Graph &G, node root);

	/*!
	 * \brief Fills the O(n log n)-space matrix with data on which basis the LCA values can be computed.
	 */
	void buildTable();

	//! \brief Access the sparse table at [i, j] for i = 0..m_len-1, j = 1..m_rangeJ
	const int &sparseTable(int i, int j) const
	{
		return m_table[i * m_rangeJ + j - 1];
	}
	int &sparseTable(int i, int j)
	{
		return m_table[i * m_rangeJ + j - 1];
	}

	/*!
	 * \brief Returns the internal index pointing to the LCA between two nodes
	 * @param u first node
	 * @param v first node
	 * @return Internal index pointing to LCA
	 */
	int rmq(int u, int v) const;
};

} // end namespace ogdf
#endif // LCA_OGDF_H_
