#include "PQPlanarity.h"
#include "utils/Logging.h"

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
