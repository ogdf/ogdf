/** \file
 * \brief Linear time layout algorithm for free trees (RadialTreeLayout).
 *
 * Based on chapter 3.1.1 Radial Drawings of Graph Drawing by
 * Di Battista, Eades, Tamassia, Tollis
 *
 * \author Carsten Gutwenger, Mirko H. Wagner
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

#include <ogdf/basic/Queue.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/tree/RadialTreeLayout.h>

namespace ogdf {

RadialTreeLayout::RadialTreeLayout()
	: m_levelDistance(30), m_selectRoot(RootSelectionType::Center) { }

RadialTreeLayout::RadialTreeLayout(const RadialTreeLayout& tl)
	: m_levelDistance(tl.m_levelDistance), m_selectRoot(tl.m_selectRoot) { }

RadialTreeLayout& RadialTreeLayout::operator=(const RadialTreeLayout& tl) {
	m_levelDistance = tl.m_levelDistance;
	m_selectRoot = tl.m_selectRoot;

	return *this;
}

void RadialTreeLayout::call(GraphAttributes& AG) {
	const Graph& tree = AG.constGraph();
	if (tree.numberOfNodes() < 2) {
		return;
	}

	OGDF_ASSERT(isArborescence(tree));

	OGDF_ASSERT(m_levelDistance > 0);

	// determine root of tree (m_root)
	FindRoot(tree);

	// compute m_level[v], m_parent[v], m_numLevels
	ComputeLevels(tree);

	// computes diameter of each node
	ComputeDiameters(AG);

	// partitions the nodes of each level by their parents
	ComputeGroupings(tree);

	// computes m_angle[v] and m_wedge[v]
	ComputeAngles(tree);

	// computes final coordinates of nodes
	ComputeCoordinates(AG);
}

void RadialTreeLayout::FindRoot(const Graph& G) {
	switch (m_selectRoot) {
	case RootSelectionType::Source:
		for (node v : G.nodes) {
			if (v->indeg() == 0) {
				m_root = v;
			}
		}
		break;

	case RootSelectionType::Sink:
		for (node v : G.nodes) {
			if (v->outdeg() == 0) {
				m_root = v;
			}
		}
		break;

	case RootSelectionType::Center: {
		NodeArray<int> degree(G);
		Queue<node> leaves;

		for (node v : G.nodes) {
			if ((degree[v] = v->degree()) == 1) {
				leaves.append(v);
			}
		}

		node v = nullptr;
		while (!leaves.empty()) {
			v = leaves.pop();

			for (adjEntry adj : v->adjEntries) {
				node u = adj->twinNode();
				if (--degree[u] == 1) {
					leaves.append(u);
				}
			}
		}
		m_root = v;
	} break;
	}
}

void RadialTreeLayout::ComputeLevels(const Graph& G) {
	m_parent.init(G, nullptr);
	m_level.init(G, -1);
	m_relWidth.init(G, 0);

	Queue<node> Q;

	Q.append(m_root);
	m_parent[m_root] = nullptr;
	m_level[m_root] = 0;

	int maxLevel = 0;

	while (!Q.empty()) {
		node v = Q.pop();
		int levelV = m_level[v];

		for (adjEntry adj : v->adjEntries) {
			node u = adj->twinNode();
			if (u == m_parent[v]) {
				continue;
			}

			Q.append(u);
			m_parent[u] = v;
			m_level[u] = maxLevel = levelV + 1;
		}
	}

	m_numLevels = maxLevel + 1;
}

void RadialTreeLayout::ComputeDiameters(GraphAttributes& AG) {
	const Graph& G = AG.constGraph();

	m_diameter.init(G, -1);
	m_nodes.init(m_numLevels);
	m_maxDiameter.init(m_numLevels);
	m_maxDiameter.fill(0);

	for (node v : G.nodes) {
		int i = m_level[v];
		m_nodes[i].pushBack(v);

		double w = AG.width(v);
		double h = AG.height(v);

		m_diameter[v] = sqrt(w * w + h * h);

		if (m_diameter[v] > m_maxDiameter[i]) {
			m_maxDiameter[i] = m_diameter[v];
		}
	}
}

void RadialTreeLayout::ComputeAngles(const Graph& G) {
	m_angle.init(G);
	m_wedge.init(G);
	m_absWidth.init(G, 0);

	m_radius.init(m_numLevels);

	Queue<node> Q;
	Q.append(m_root);
	m_angle[m_root] = 0;
	m_wedge[m_root] = 2 * Math::pi;
	m_radius[0] = 0;

	int iProcessed = 0;

	while (!Q.empty()) {
		node v = Q.pop();
		node p = m_parent[v];

		// nothing to do if v is a leaf
		if (v->degree() == 1 && p) {
			continue;
		}

		int i = m_level[v];
		if (i + 1 > iProcessed) {
			m_radius[i + 1] =
					m_radius[i] + 0.5 * (m_maxDiameter[i + 1] + m_maxDiameter[i]) + m_levelDistance;

			for (node w : m_nodes[i]) {
				if (w->degree() == 1 && m_parent[w]) {
					continue;
				}

				for (node u : m_children[w]) {
					m_absWidth[w] += m_relWidth[u] * (i + 1);
				}
				// radius needed to allow for the given wedge to have this width
				m_radius[i + 1] = max(m_absWidth[w] / m_wedge[w], m_radius[i + 1]);
			}
			iProcessed = i + 1;
		}

		double offset = m_angle[v] - 0.5 * m_wedge[v];

		for (node u : m_children[v]) {
			if (u->degree() != 1) {
				Q.append(u);
			}

			double desiredWedge = m_relWidth[u] * (i + 1) * m_wedge[v] / m_absWidth[v];
			double allowedWedge = 2 * acos(m_radius[i] / m_radius[i + 1]);
			m_wedge[u] = min(desiredWedge, allowedWedge);
			m_angle[u] = offset + 0.5 * m_wedge[u];

			offset += m_wedge[u];
		}
	}
}

// group children of nodes on level i (i.e. nodes of level i+1)
void RadialTreeLayout::ComputeGroupings(const Graph& G) {
	m_children.init(G);

	NodeArray<int> degree(G);
	Queue<node> leaves;

	for (int i = 0; i < m_numLevels; i++) {
		for (node v : m_nodes[i]) {
			node p = m_parent[v];

			if ((degree[v] = v->degree()) == 1 && p) {
				leaves.append(v);
				continue;
			}

			adjEntry adj = v->firstAdj();
			adjEntry adjStop = v->firstAdj();
			if (p) {
				while (adj->twinNode() != p) {
					adj = adj->cyclicSucc();
				}
				adjStop = adj;
				adj = adj->cyclicSucc();
			}

			do {
				node u = adj->twinNode();
				m_children[v].pushBack(u);
				adj = adj->cyclicSucc();
			} while (adj != adjStop);
		}
	}

	while (!leaves.empty()) {
		node v = leaves.pop();
		m_relWidth[v] = 0;

		for (adjEntry adj : v->adjEntries) {
			node u = adj->twinNode();
			if (u != m_parent[v]) {
				m_relWidth[v] += m_relWidth[u];
			} else if (--degree[u] == (v == m_root ? 0 : 1)) {
				leaves.append(u);
			}
		}
		if (v != m_root) {
			// the subtree at v is either as wide as v's own width
			// or as the sum of its subtrees widths
			m_relWidth[v] = max(m_relWidth[v], (m_diameter[v] + m_levelDistance) / m_level[v]);
		}
	}
}

void RadialTreeLayout::ComputeCoordinates(GraphAttributes& AG) {
	const Graph& G = AG.constGraph();

	for (node v : G.nodes) {
		double r = m_radius[m_level[v]];
		double alpha = m_angle[v];

		AG.x(v) = r * cos(alpha);
		AG.y(v) = r * sin(alpha);
	}

	AG.clearAllBends();
}
}
