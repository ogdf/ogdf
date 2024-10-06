/** \file
 * \brief Implementation of a utility class representing a cycle in the Blossom
 * algorithm.
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

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/basic.h>
#include <ogdf/graphalg/matching_blossom/Cycle.h>

#include <cstddef>
#include <tuple>
#include <unordered_set>
#include <vector>

namespace ogdf {
namespace matching_blossom {

Cycle::Cycle(edge startEdge) { addEdge(startEdge); }

const std::vector<edge>& Cycle::edgeOrder() { return m_edgeOrder; }

const std::unordered_set<node>& Cycle::nodes() { return m_nodes; }

void Cycle::addEdge(edge e) {
	node commonNode;
	if (m_edgeOrder.size()) {
		commonNode = e->commonNode(m_edgeOrder[m_edgeOrder.size() - 1]);
		OGDF_ASSERT(commonNode);
	} else {
		m_nodes.insert(e->source());
		commonNode = e->source();
	}
	m_nodes.insert(e->opposite(commonNode));
	m_edgeOrder.push_back(e);
}

node Cycle::startNode() { return m_edgeOrder[0]->commonNode(m_edgeOrder[m_edgeOrder.size() - 1]); }

std::vector<long> Cycle::indexOf(std::vector<node> nodesToFind) {
	std::vector<long> indices(nodesToFind.size(), -1);
	for (size_t i = 0; i < nodesToFind.size(); ++i) {
		if (nodesToFind[i] == startNode()) {
			// if node is equal to cycle start node, skip the search and set indices to the last edge
			indices[i] = m_edgeOrder.size() - 1;
		}
	}
	for (size_t i = 0; i < m_edgeOrder.size(); ++i) {
		for (size_t j = 0; j < nodesToFind.size(); ++j) {
			if (indices[j] < 0 && m_edgeOrder[i]->isIncident(nodesToFind[j])) {
				indices[j] = i;
			}
		}
	}
	return indices;
}

std::tuple<long, long> Cycle::indexOf(node u, node v) {
	auto indices = indexOf({u, v});
	return std::make_tuple(indices[0], indices[1]);
}

long Cycle::indexOf(node u) {
	std::vector<node> nodesToFind = {u};
	return indexOf(nodesToFind)[0];
}

bool Cycle::contains(node v) { return m_nodes.find(v) != m_nodes.end(); }

}
}
