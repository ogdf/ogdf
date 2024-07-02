/** \file
 * \brief A simple solver for 2-SAT instances.
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
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/simple_graph_alg.h>

#include <utility>
#include <vector>

namespace ogdf {

#ifdef OGDF_DEBUG

//! In debug mode, twosat_var is a class instead of a simple int to prevent unintened use of the default 0-value instead of TwoSAT_Var_Undefined
class OGDF_EXPORT twosat_var {
	int m_val;

public:
	OGDF_DEPRECATED("uninitialized value")

	twosat_var() : m_val(0) { }

	explicit twosat_var(int val) : m_val(val) { }

	explicit operator int() { return m_val; }

	bool operator==(const twosat_var& rhs) const { return m_val == rhs.m_val; }

	bool operator!=(const twosat_var& rhs) const { return m_val != rhs.m_val; }
};

const twosat_var TwoSAT_Var_Undefined(-1);

#else
using twosat_var = int;
const twosat_var TwoSAT_Var_Undefined = -1;
#endif

//! A simple solver for TwoSAT instances, representing the instance as implication graph and solving it via its strongly-connected components.
class OGDF_EXPORT TwoSAT : protected Graph {
	std::vector<bool> m_assignment;
	std::vector<node> m_node_map;

public:
	TwoSAT() {};

	twosat_var newVariable() {
		node pos = newNode(), neg = newNode();
		OGDF_ASSERT(pos->index() % 2 == 0);
		OGDF_ASSERT(pos->index() + 1 == neg->index());
		m_assignment.push_back(false);
		OGDF_ASSERT(m_node_map.size() == pos->index());
		m_node_map.push_back(pos);
		OGDF_ASSERT(m_node_map.size() == neg->index());
		m_node_map.push_back(neg);
		return (twosat_var)(pos->index() / 2);
	}

	int variableCount() const {
		OGDF_ASSERT(nodes.size() % 2 == 0);
		return nodes.size() / 2;
	}

	int newClause(twosat_var var1, bool sign1, twosat_var var2, bool sign2) {
		auto [p1, n1] = getNodes(var1, !sign1);
		auto [p2, n2] = getNodes(var2, !sign2);
		edge e1 = newEdge(n1, p2);
		edge e2 = newEdge(n2, p1);
		OGDF_ASSERT(e1->index() % 2 == 0);
		OGDF_ASSERT(e1->index() + 1 == e2->index());
		return e1->index();
	}

	int clauseCount() const {
		OGDF_ASSERT(edges.size() % 2 == 0);
		return edges.size() / 2;
	}

	bool solve() {
		NodeArray<int> components(*this, -1);
		strongComponents(*this, components);

		OGDF_ASSERT(nodes.size() % 2 == 0);
		OGDF_ASSERT(edges.size() % 2 == 0);
		for (auto it = nodes.begin(); it != nodes.end(); it++) {
			node pos = *(it);
			it++;
			node neg = *(it);
			OGDF_ASSERT(pos->index() % 2 == 0);
			OGDF_ASSERT(pos->index() + 1 == neg->index());

			if (components[pos] == components[neg]) {
				return false;
			}
			m_assignment[pos->index() / 2] = (components[pos] > components[neg]);
		}
		return true;
	}

	bool getAssignment(twosat_var var) const { return m_assignment[(int)var]; }

	void clear() override {
		Graph::clear();
		m_assignment.clear();
	}

	using const_iterator = typename std::vector<bool>::const_iterator;
	using const_reverse_iterator = typename std::vector<bool>::const_reverse_iterator;

	const_iterator begin() const { return m_assignment.begin(); }

	const_iterator end() const { return m_assignment.end(); }

	const_reverse_iterator rbegin() const { return m_assignment.rbegin(); }

	const_reverse_iterator rend() const { return m_assignment.rend(); }

protected:
	std::pair<node, node> getNodes(twosat_var var, bool flip = false) const {
		int idx = getNodeIndex(var, true);
		node pos = m_node_map[idx];
		OGDF_ASSERT(pos->index() == idx);
		OGDF_ASSERT(pos->index() / 2 == (int)var);
		OGDF_ASSERT(pos->index() % 2 == 0);
		node neg = m_node_map[idx + 1];
		OGDF_ASSERT(neg->index() == idx + 1);
		if (flip) {
			return std::pair(neg, pos);
		} else {
			return std::pair(pos, neg);
		}
	}

	node getNode(twosat_var var, bool sign) const {
		int idx = getNodeIndex(var, sign);
		node n = m_node_map[idx];
		OGDF_ASSERT(n->index() == idx);
		return n;
	}

	int getNodeIndex(twosat_var var, bool sign) const {
		int idx = ((int)var) * 2 + (sign ? 0 : 1);
		OGDF_ASSERT(0 <= idx);
		OGDF_ASSERT(idx <= Graph::maxNodeIndex());
		return idx;
	}
};

}
