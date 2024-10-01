/** \file
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
#include <ogdf/basic/List.h>
#include <ogdf/basic/basic.h>
#include <ogdf/cluster/ClusterGraph.h>
#include <ogdf/cluster/sync_plan/utils/Bijection.h>
#include <ogdf/cluster/sync_plan/utils/Logging.h>

#include <functional>
#include <sstream>
#include <string>
#include <utility>


using namespace ogdf::sync_plan::internal;

namespace ogdf::sync_plan::internal {

std::string to_string(const std::function<std::ostream&(std::ostream&)>& func) {
	std::stringstream ss;
	const std::ostream& ret = func(ss);
	OGDF_ASSERT(&ret == &ss);
	return ss.str();
}

std::ostream& operator<<(std::ostream& os, const std::function<std::ostream&(std::ostream&)>& func) {
	return func(os);
}

std::ostream& operator<<(std::ostream& os, const Graph& G) {
	return os << "Graph with " << G.numberOfNodes() << " nodes, " << G.numberOfEdges() << " edges";
}

std::ostream& operator<<(std::ostream& os, const ClusterGraph& CG) {
	return os << "ClusterGraph with " << CG.constGraph().numberOfNodes() << " nodes, "
			  << CG.constGraph().numberOfEdges() << " edges and " << CG.numberOfClusters()
			  << " clusters";
}

template<>
std::ostream& operator<<(std::ostream& os, const printIncidentEdges<PipeBij>& inst) {
	for (const auto& pair : inst.container) {
		ogdf::adjEntry adj = pair.first;
		os << "e" << adj->theEdge()->index() << " (" << (adj->isSource() ? ">" : "<") << "n"
		   << adj->twinNode()->index() << "), ";
	}
	return os;
}

template<>
std::ostream& operator<<(std::ostream& os, const printEdges<PipeBij>& inst) {
	for (const auto& pair : inst.container) {
		ogdf::adjEntry adj = pair.first;
		os << "e" << adj->theEdge()->index() << " (n" << adj->theNode()->index()
		   << (adj->isSource() ? "->" : "<-") << "n" << adj->twinNode()->index() << "), ";
	}
	return os;
}

}
