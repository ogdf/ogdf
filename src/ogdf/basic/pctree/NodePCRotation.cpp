#include "NodePCRotation.h"

#include <ogdf/basic/Queue.h>
#include <ogdf/basic/STNumbering.h>

#include "utils/FilteringBFS.h"

using namespace ogdf;
using namespace pc_tree;
using namespace Dodecahedron;

using RegisteredEdgeSet = RegisteredElementSet<edge, Graph>;

NodePCRotation::NodePCRotation(const Graph& G, node end, const bool mapBundleEdges)
	: m_G(&G), m_n(end), incidentEdgeForLeaf(*this, nullptr), graphNodeForInnerNode(*this, nullptr) {
	// Compute st-numbering for graph and order the nodes accordingly.
	NodeArray<int> numbering(G);
	int nodecount = computeSTNumbering(G, numbering, nullptr, end);
	OGDF_ASSERT(nodecount > 0);
	Array<node> order(nodecount);
	int maxDegree = 1;
	for (node n : FilteringBFS(G, {end})) {
		order[numbering[n] - 1] = n;
		if (n->degree() > maxDegree) {
			maxDegree = n->degree();
		}
	}

	std::vector<PCNode*> consecutiveLeaves;
	consecutiveLeaves.reserve(maxDegree);
	std::vector<edge> outEdges;
	outEdges.reserve(maxDegree);
	EdgeArray<PCNode*> leafRepresentation(G);
	List<edge> twinPoleCandidateInEdges;
	RegisteredEdgeSet bundleEdges;
	if (mapBundleEdges) {
		bundleEdges.init(G);
		bundleEdgesForLeaf.init(*this);
	}
	PCNode* stEdgeLeaf = nullptr;
	for (node n : order) {
		if (n == order[order.size() - 1]) {
			// Graph is planar

			if (mapBundleEdges && isTrivial()) {
				// Map all incoming edges of the partner node to the st-edge, unless the partner node is the first node
				if (!twinPoleCandidateInEdges.empty()) {
					bundleEdgesForLeaf[stEdgeLeaf].swap(twinPoleCandidateInEdges);
				}
			}

			return;
		}

		consecutiveLeaves.clear();
		outEdges.clear();
		// Separate the incident edges to nodes with a higher number from the edges to nodes with a lower number.
		for (adjEntry adj : n->adjEntries) {
			if (numbering[adj->theEdge()->opposite(n)] > numbering[n]) {
				outEdges.push_back(adj->theEdge());
			} else {
				consecutiveLeaves.push_back(leafRepresentation[adj->theEdge()]);
			}
		}

		PCNode* mergedLeaf;
		bool twinPoleCandidate = false;
		if (mapBundleEdges) {
			twinPoleCandidate = true;
			bundleEdges.clear();
		}
		if (n == order[0]) {
			// Create root node for first node of graph.
			mergedLeaf = newNode(PCNodeType::PNode);
		} else {
			// Make all edges to nodes with a lower number consecutive and merge them to a single leaf.

			if (mapBundleEdges) {
				if (consecutiveLeaves.size() + 1 < getLeafCount()) {
					twinPoleCandidate = false;
				} else {
					twinPoleCandidateInEdges.clear();
				}

				for (PCNode* v : consecutiveLeaves) {
					for (edge e : bundleEdgesForLeaf[v]) {
						bundleEdges.insert(e);
					}

					if (twinPoleCandidate) {
						twinPoleCandidateInEdges.pushBack(incidentEdgeForLeaf[v]);
					}
				}
			}

			mergedLeaf = mergeLeaves(consecutiveLeaves);

			if (mapBundleEdges && twinPoleCandidate) {
				if (getLeaves().front() == mergedLeaf) {
					stEdgeLeaf = getLeaves().back();
				} else {
					stEdgeLeaf = getLeaves().front();
				}
			}
		}

		if (mergedLeaf == nullptr) {
			throw GraphNotPlanarException();
		}
		graphNodeForInnerNode[mergedLeaf] = n;

		if (outEdges.size() > 1) {
			consecutiveLeaves.clear(); // re-use consecutiveLeaves
			std::vector<PCNode*>& addedLeaves = consecutiveLeaves;
			if (mergedLeaf->getNodeType() == PCNodeType::Leaf) {
				if (getLeafCount() <= 2) {
					// see also replaceLeaf
					PCNode* parent = mergedLeaf->getParent();
					OGDF_ASSERT(parent == getRootNode());

					mergedLeaf->detach();
					destroyNode(mergedLeaf);

					mergedLeaf = parent;
					graphNodeForInnerNode[mergedLeaf] = n;
				}
				changeNodeType(mergedLeaf, PCNodeType::PNode);
			}
			insertLeaves(outEdges.size(), mergedLeaf, &addedLeaves);

			for (size_t i = 0; i < outEdges.size(); i++) {
				// For each edge store which leaf represents it.
				edge e = outEdges[i];
				PCNode* newLeaf = addedLeaves[i];
				leafRepresentation[e] = newLeaf;
				incidentEdgeForLeaf[newLeaf] = e;

				if (mapBundleEdges) {
					if (twinPoleCandidate) {
						bundleEdgesForLeaf[newLeaf] = {e};
					} else {
						bundleEdgesForLeaf[newLeaf] = bundleEdges.elements();
					}
				}
			}
		} else if (outEdges.size() == 1) {
			// If there is only one edge to a node with a higher number,
			// declare the mergedLeaf as the leaf that represents it.
			leafRepresentation[outEdges.front()] = mergedLeaf;
			incidentEdgeForLeaf[mergedLeaf] = outEdges.front();

			if (mapBundleEdges) {
				if (twinPoleCandidate) {
					bundleEdgesForLeaf[mergedLeaf] = {outEdges.front()};
				} else {
					bundleEdgesForLeaf[mergedLeaf] = bundleEdges.elements();
				}
			}
		}

		OGDF_HEAVY_ASSERT(checkValid());
	}
	OGDF_ASSERT(false); // should either return from processing last node or throw non-planar
}

node NodePCRotation::getTrivialPartnerPole() const {
	if (m_n->degree() < 3) {
		return nullptr;
	}
	if (!isTrivial()) {
		return nullptr;
	}
	node partner = graphNodeForInnerNode[getRootNode()];
	if (partner->degree() < 3) {
		return nullptr;
	} else {
		return partner;
	}
}

node NodePCRotation::getTrivialPartnerPole(const Graph& G, node pole) {
	NodePCRotation pc(G, pole, false);
	return pc.getTrivialPartnerPole();
}

bool NodePCRotation::isEqual(const NodePCRotation& pc) const {
	OGDF_ASSERT(checkValid());
	OGDF_ASSERT(pc.checkValid());
	if (getGraph() != pc.getGraph()) {
		return false;
	}
	if (getNode() != pc.getNode()) {
		return false;
	}
	if (getLeafCount() != pc.getLeafCount()) {
		return false;
	}
	if (getPNodeCount() != pc.getPNodeCount()) {
		return false;
	}
	if (getCNodeCount() != pc.getCNodeCount()) {
		return false;
	}
	if (possibleOrders() != pc.possibleOrders()) {
		return false;
	}

	EdgeArray<PCNode*> mapping(*m_G, nullptr);
	pc.generateLeafForIncidentEdgeMapping(mapping);
	std::vector<PCNode*> order;
	order.reserve(getLeafCount());
	for (PCNode* l : currentLeafOrder()) {
		order.push_back(mapping[getIncidentEdgeForLeaf(l)]);
	}
	if (!pc.isValidOrder(order)) {
		return false;
	}

#ifdef OGDF_HEAVY_DEBUG
	generateLeafForIncidentEdgeMapping(mapping);
	order.clear();
	for (PCNode* l : pc.currentLeafOrder()) {
		order.push_back(mapping[pc.getIncidentEdgeForLeaf(l)]);
	}
	OGDF_HEAVY_ASSERT(isValidOrder(order));
#endif

	return true;
}

std::function<void(std::ostream& os, pc_tree::PCNode*, int)> NodePCRotation::uidPrinter() const {
	return [&](std::ostream& os, pc_tree::PCNode* n, int i) -> void {
		if (n->isLeaf()) {
			os << getIncidentEdgeForLeaf(n)->index();
		} else if (getGraphNodeForInnerNode(n) != nullptr) {
			os << getGraphNodeForInnerNode(n)->index();
		} else {
			os << n->getNodeType();
		}
	};
}

std::function<bool(pc_tree::PCNode*, pc_tree::PCNode*)> NodePCRotation::uidComparer() const {
	return [&](pc_tree::PCNode* a, pc_tree::PCNode* b) -> bool {
		if (a->isLeaf() && b->isLeaf()) {
			return getIncidentEdgeForLeaf(a)->index() < getIncidentEdgeForLeaf(b)->index();
		} else if (a->isLeaf()) {
			return true;
		} else if (b->isLeaf()) {
			return false;
		}

		node ia = getGraphNodeForInnerNode(a);
		node ib = getGraphNodeForInnerNode(b);
		if (ia != nullptr && ib != nullptr) {
			return ia->index() < ib->index();
		} else if (ia != nullptr) {
			return true;
		} else if (ib != nullptr) {
			return false;
		}

		if (a->getDegree() < b->getDegree()) {
			return true;
		} else if (a->getDegree() > b->getDegree()) {
			return false;
		} else {
			return true; // can't compare C-nodes of same degree easily...
		}
	};
}
