/** \file
 * \brief Applies the node coloring approximation specified by Johnson.
 *
 * \author Jan-Niklas Buckow
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
#include <ogdf/graphalg/NodeColoringModule.h>

namespace ogdf {

/**
 * Approximation algorithms for the node coloring problem in graphs.
 * This class implements the approximation given by Johnson which tries
 * to find heuristically large independent sets one other the other.
 */
class OGDF_EXPORT NodeColoringJohnson : public NodeColoringModule {
public:
	/**
	 * Declares procedure to find the minimum degree nodes.
	 */
	enum class MinDegreeProcedure {
		smallestIndex, ///< Use the node with the smallest index
		minDegreeOriginal, ///< Use the node which has the smallest degree in the original graph
		maxDegreeOriginal, ///< Use the node which has the largest degree in the original graph
		minDegreeSubgraph, ///< Use the node which has the smallest degree in the current subgraph
		maxDegreeSubgraph ///< Use the node which has the largest degree in the current subgraph
	};

	/**
	 * The constructor.
	 * Initializes the procedure of finding minimum degree nodes with the smallest index.
	 */
	NodeColoringJohnson() : m_minDegreeProcedure(MinDegreeProcedure::maxDegreeSubgraph) { }

	/**
	 * Sets the procedure of finding minimum degree nodes.
	 * @param minDegreeProcedure The desired minimum degree finding procedure
	 */
	inline void setMinDegreeProcedure(MinDegreeProcedure minDegreeProcedure) {
		m_minDegreeProcedure = minDegreeProcedure;
	}

	virtual NodeColor call(const Graph& graph, NodeArray<NodeColor>& colors,
			NodeColor start = 0) override;

private:
	MinDegreeProcedure m_minDegreeProcedure;
};
}
