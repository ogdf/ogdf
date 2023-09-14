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

#include <ogdf/basic/CombinatorialEmbedding.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphAttributes.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace tools {

inline bool equal(const node a, const node b) { return a->index() == b->index(); }

inline bool equal(const edge& a, const edge& b) {
	return equal(a->source(), b->source()) && equal(a->target(), b->target());
}

inline bool equal(const face& a, const face& b) { return a->index() == b->index(); }

inline std::vector<node> nodes_of_face(const face& face) {
	std::vector<node> nodes;
	adjEntry first = face->firstAdj();
	adjEntry current = first;

	do {
		nodes.push_back(current->theEdge()->source());
		current = face->nextFaceEdge(current);
	} while (current != first && current != NULL);
	return nodes;
}

}
}
}
}
