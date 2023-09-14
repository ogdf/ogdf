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

#	include <ogdf/basic/CombinatorialEmbedding.h>
#	include <ogdf/basic/Graph.h>
#	include <ogdf/basic/GraphAttributes.h>
#	include <ogdf/geometric/cr_min/geometry/objects/LineSegment.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Point.h>

#	include <map>

namespace ogdf {
namespace internal {
namespace gcm {
namespace tools {

template<typename kernel>
geometry::Point_t<kernel> extract_point(const GraphAttributes& ga, const node node) {
	return {ga.x(node), ga.y(node)};
}

template<typename kernel>
geometry::LineSegment_t<kernel> extract_segment(const GraphAttributes& ga, const edge& edge) {
	return {extract_point<kernel>(ga, edge->source()), extract_point<kernel>(ga, edge->target())};
}

template<typename Graph, typename Polygon>
void extract_polygon(const Graph& g, const face& face, Polygon& p) {
	adjEntry first = face->firstAdj();
	adjEntry current = first;
	do {
		p.push_back(g.get_point(current->theNode()));
		current = current->faceCycleSucc();
		OGDF_ASSERT(current != nullptr);
	} while (current != first);
}

template<typename kernel, typename Polygon>
void extract_polygon(const GraphAttributes& ga, const face& face,
		NodeArray<unsigned int>& node_to_id, Polygon& p) {
	adjEntry first = face->firstAdj();
	adjEntry current = first;

	do {
		node_to_id[current->theEdge()->source()] = p.size();
		p.push_back(extract_point<kernel>(ga, current->theEdge()->source()));
		current = face->nextFaceEdge(current);
	} while (current != nullptr);
}

}
}
}
}

#endif
