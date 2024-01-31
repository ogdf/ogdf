/** \file
 * \brief Utility class representing a cycle in the Blossom algorithm.
 *
 * \author Joshua Sangmeister
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

#include <tuple>
#include <unordered_set>
#include <vector>

namespace ogdf {
namespace matching_blossom {

class OGDF_EXPORT Cycle {
private:
	std::unordered_set<node> m_nodes;

	std::vector<edge> m_edgeOrder;

	//! Get the indices of the edges before \p nodesToFind in edge order. For the startNode,
	//! `size() - 1` is returned.
	std::vector<long> indexOf(std::vector<node> nodesToFind);

public:
	Cycle(edge startEdge);

	/* Getters */

	const std::vector<edge>& edgeOrder();

	const std::unordered_set<node>& nodes();

	/* End of getters */

	//! Use this method to add edges in cycle order.
	void addEdge(edge e);

	//! The first node of the cycle in edge order.
	node startNode();

	//! Get the indices of the edges before \p u and \p v in edge order. For the startNode,
	//! `size() - 1` is returned.
	std::tuple<long, long> indexOf(node u, node v);

	//! Get the index of the edge before \p u in edge order. If \p u is the startNode,
	//! `size() - 1` is returned.
	long indexOf(node u);

	//! Whether the cycle contains the node \p v or not.
	bool contains(node v);
};

}
}
