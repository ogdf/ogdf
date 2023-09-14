/** \file
 *
 * \author Marcel Radermacher
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

#ifdef OGDF_INCLUDE_CGAL

#	include <ogdf/geometric/cr_min/graph/OGDFGraphWrapper.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace graph {

class Path {
private:
	using Graph = OGDFGraphWrapper;
	using Node = typename Graph::Node;
	using Edge = typename Graph::Edge;
	std::vector<Node> m_nodes;
	std::vector<Edge> m_edges;

public:
	inline const std::vector<Node>& nodes() const { return m_nodes; }

	inline const std::vector<Edge>& edges() const { return m_edges; }

	inline bool is_reversed(unsigned int edge_i) const {
		return m_edges[edge_i]->source() != m_nodes[edge_i];
	}

	void push_back(const Node& src, const Edge& e) {
		OGDF_ASSERT(m_nodes.empty() || e->isIncident(m_nodes.back()));
		OGDF_ASSERT(m_nodes.empty() || m_nodes.back() == src);
		m_edges.push_back(e);
		if (!m_nodes.empty()) {
			m_nodes.push_back(e->opposite(m_nodes.back()));
		} else {
			m_nodes.push_back(src);
			m_nodes.push_back(e->opposite(src));
		}
	}

	void clear() {
		m_nodes.clear();
		m_edges.clear();
	}

	inline bool empty() const { return m_nodes.empty(); }

	void reverse() {
		std::reverse(m_nodes.begin(), m_nodes.end());
		std::reverse(m_edges.begin(), m_edges.end());
	}

	void print() {
		if (!empty()) {
			for (unsigned int i = 0; i < m_edges.size(); ++i) {
				std::cout << m_nodes[i] << ", " << m_edges[i] << "->";
			}
			std::cout << m_nodes.back() << std::endl;
		} else {
			std::cout << "path is empty" << std::endl;
		}
	}
};


}
}
}
}

#endif
