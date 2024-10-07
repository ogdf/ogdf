/** \file
 * \brief Trivial node coloring which assigns every node a different color.
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

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/basic.h>
#include <ogdf/graphalg/NodeColoringModule.h>
#include <ogdf/graphalg/NodeColoringSimple.h>

namespace ogdf {

using NColor = NodeColoringModule::NodeColor;

NColor NodeColoringSimple::call(const Graph& graph, NodeArray<NColor>& colors, NColor start) {
	NColor numberColorsUsed = NColor(graph.numberOfNodes());
	for (node v : graph.nodes) {
		colors[v] = start++;
	}
	OGDF_ASSERT(checkColoring(graph, colors));
	return numberColorsUsed;
}

}
