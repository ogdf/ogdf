/** \file
 * \brief TODO Document
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
#include <ogdf/basic/Graph.h>

using namespace ogdf;

class DFS {
	Graph G;

	NodeArray<int> discovery;
	NodeArray<int> finish;
	NodeArray<node> predecessor;

	int time = 0;

public:
	explicit DFS(const Graph& g)
		: G(g), discovery(G, -1), finish(G, -1), predecessor(G, nullptr) { }

	void dfs_visit(node u) {
		time += 1;
		discovery[u] = time;

		for (adjEntry adj : u->adjEntries) {
			node v = adj->twinNode();
			if (adj->isSource() && discovery[v] < 0) {
				predecessor[v] = u;
				dfs_visit(v);
			}
		}

		time += 1;
		finish[u] = time;
	}

	void dfs() {
		for (node n : G.nodes) {
			if (discovery[n] < 0) {
				dfs_visit(n);
			}
		}
	}
};

int main() { }
