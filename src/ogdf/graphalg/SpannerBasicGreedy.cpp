/** \file
 * \brief Implementation of the basic greedy (2k-1)-spanner
 * algorithm of Althöfer et al. 2007.
 *
 * \author Finn Stutzenstein
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

#include <ogdf/graphalg/SpannerBasicGreedy.h>

namespace ogdf {

// We need some rounding for maxLookupDist for the integer case:
// maxLookupDist is a double, but the dijkstra implementation reuires an int, so we have to
// ceil the maxLookupDist to not make mistakes when rounding down.
template<>
double SpannerBasicGreedy<int>::distanceInSpanner(node s, node t, double maxLookupDist) {
	NodeArray<int> distances;
	NodeArray<edge> predecessor;
	Dijkstra<int> dijkstra;
	dijkstra.callBound(*m_spanner, m_spannerWeights, s, predecessor, distances,
			false, // directed
			false, // arcs reversed
			t, ceil(maxLookupDist));
	return distances[t];
}

template<>
double SpannerBasicGreedy<double>::distanceInSpanner(node s, node t, double maxLookupDist) {
	NodeArray<double> distances;
	NodeArray<edge> predecessor;
	Dijkstra<double> dijkstra;
	dijkstra.callBound(*m_spanner, m_spannerWeights, s, predecessor, distances,
			false, // directed
			false, // arcs reversed
			t, maxLookupDist);
	return distances[t];
}

}
