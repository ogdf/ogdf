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
#pragma once
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/Logger.h>

#include <iosfwd>
#include <string>
#include <vector>

namespace ogdf {
class ClusterGraph;
class GraphAttributes;
} // namespace ogdf

#pragma GCC diagnostic ignored "-Wshadow" // TODO remove

namespace ogdf::sync_plan {

void writeLevelGraph(const Graph& G, const std::vector<std::vector<node>>& emb,
		const NodeArray<int>& pos, std::ostream& os);

void readLevelGraph(Graph& G, std::vector<std::vector<node>>& emb, NodeArray<int>& pos,
		std::istream& is);

void layout(const Graph& G, const NodeArray<int>& lvl, const NodeArray<int>& pos,
		GraphAttributes& GA);

void drawSVG(const Graph& G, const NodeArray<int>& lvl, const NodeArray<int>& pos,
		const string& name);

struct AdjCompLess {
	const NodeArray<int>& pos;

	explicit AdjCompLess(const NodeArray<int>& pos) : pos(pos) { }

	bool operator()(adjEntry a, adjEntry b);
};

void embedPLE(Graph& G, const std::vector<std::vector<node>>& emb, NodeArray<int>& lvl,
		NodeArray<int>& pos);

void checkPLE(const Graph& G, const std::vector<std::vector<node>>& emb, const NodeArray<int>& lvl,
		const NodeArray<int>& pos);

void randomProperMaximalLevelPlaneGraph(Graph& G, std::vector<std::vector<node>>& emb, int N, int K,
		bool radial);

void pruneEdges(Graph& G, int max_edges, int min_deg);

void reduceLevelToCluster(const Graph& LG, const std::vector<std::vector<node>>& emb, Graph& G,
		ClusterGraph& CG);

}
