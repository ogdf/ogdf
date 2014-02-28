/*
 * $Revision: 3931 $
 *
 * last checkin:
 *   $Author: beyer $
 *   $Date: 2014-02-20 14:56:42 +0100 (Do, 20. Feb 2014) $
 ***************************************************************/

/** \file
 * \brief Declares class AdjacencyOracle.
 *
 * \author Rene Weiskircher
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


#ifdef _MSC_VER
#pragma once
#endif

#ifndef OGDF_ADJACENCY_ORACLE_H
#define OGDF_ADJACENCY_ORACLE_H

#include <ogdf/basic/NodeArray.h>
#include <ogdf/basic/Array2D.h>

namespace ogdf {

//! Tells you in constant time if two nodes are adjacent
/**
 * AdjacencyOracle is initialized with a Graph and returns for
 * any pair of nodes in constant time if they are adajcent.
 */
class AdjacencyOracle {
public:
	//! The constructor for the class, needs time O(n + m)
	AdjacencyOracle(const Graph &G);
	//! The destructor
	~AdjacencyOracle() { }
	//! This returns true if the two nodes are adjacent in G, false otherwise
	bool adjacent(const node, const node) const;
private:
	NodeArray<int> m_nodeNum; //!< The internal number given to each node
	Array2D<bool> m_adjacencyMatrix; //!< A 2D-array where the entry is true if the nodes with the corresponding number are adjacent
};

}
#endif
