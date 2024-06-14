#include "utils/Logging.h"

#include <sstream>

using namespace ogdf;

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

std::ostream& operator<<(std::ostream& os, const BCTree::BNodeType& obj) {
	if (obj == BCTree::BNodeType::CComp) {
		os << "cut";
	} else if (obj == BCTree::BNodeType::BComp) {
		os << "bicon";
	} else {
		os << "???";
	}
	return os;
}

std::ostream& operator<<(std::ostream& os, const BCTree::GNodeType& obj) {
	if (obj == BCTree::GNodeType::CutVertex) {
		os << "cut-vertex";
	} else if (obj == BCTree::GNodeType::Normal) {
		os << "block-vertex";
	} else {
		os << "???";
	}
	return os;
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
