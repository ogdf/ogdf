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

#include <array>
#include <functional>

namespace ogdf {
class Color;
class GraphAttributes;

extern const std::array<Color, 63> colors;

void spreadParallels(GraphAttributes& GA, double min_spread = 0.1, double max_spread = 0.6,
		double max_abs = 100);

void fixLoops(Graph& G, const std::function<void(edge, edge)>& cb);

void fixParallels(Graph& G, const std::function<void(edge, edge)>& cb);

void bendEdge(GraphAttributes& GA, edge e, double bend);

}
