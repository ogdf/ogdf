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
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/basic.h>
#include <ogdf/cluster/sync_plan/PMatching.h>
#include <ogdf/cluster/sync_plan/PQPlanarity.h>
#include <ogdf/cluster/sync_plan/PQPlanarityComponents.h>
#include <ogdf/cluster/sync_plan/QPartitioning.h>

#include <functional>
#include <ostream>
#include <string>

std::ostream& operator<<(std::ostream& os, const PQPlanarity& pq) {
	return os << "PQPlanarity Instance with " << pq.G->numberOfNodes() << " nodes, "
			  << pq.G->numberOfEdges() << " edges, " << pq.matchings.getPipeCount() << " pipes, "
			  << pq.partitions.qVertexCount() << " Q-Vertices in " << pq.partitions.partitionCount()
			  << " partitions and " << pq.components.connectedCount() << " connected components";
}

std::ostream& operator<<(std::ostream& os, const PQPlanarity::UndoOperation& undo_op) {
	return undo_op.print(os);
}

std::function<std::ostream&(std::ostream&)> PQPlanarity::fmtPQNode(node n, bool include_comp) const {
	OGDF_ASSERT(n == nullptr || n->graphOf() == G);
	return [n, include_comp, this](std::ostream& ss) -> std::ostream& {
		if (n != nullptr) {
			ss << "["
			   << (matchings.isMatchedPVertex(n) ? "mP" : (partitions.isQVertex(n) ? "Q" : "uP"));
			if (include_comp) {
				ss << (components.isCutVertex(n) ? "C" : "B");
			}
			ss << " #" << n->index() << " °" << n->degree();
			if (include_comp) {
				ss << " @" << components.biconnectedId(n) << "/" << components.connectedId(n);
			}
			if (GA != nullptr) {
				ss << " \"" << GA->label(n) << "\"";
			}
			ss << "]";
		} else {
			ss << "[NULL]";
		}
		return ss;
	};
}
