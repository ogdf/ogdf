/** \file
 * \brief Implementation for ogdf::pc_tree::NodePCRotation
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

#include <ogdf/basic/Array.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/GraphSets.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/STNumbering.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/pctree/NodePCRotation.h>
#include <ogdf/basic/pctree/PCEnum.h>
#include <ogdf/basic/pctree/PCNode.h>
#include <ogdf/basic/pctree/PCRegistry.h>
#include <ogdf/basic/pctree/util/FilteringBFS.h>
#include <ogdf/basic/pctree/util/IntrusiveList.h>

#include <cstddef>
#include <functional>
#include <ostream>
#include <vector>

using namespace ogdf::pc_tree;

NodePCRotation::NodePCRotation(const Graph& G, node end, const bool mapBundleEdges)
	: m_G(&G)
	, m_n(end)
	, m_incidentEdgeForLeaf(*this, nullptr)
	, m_graphNodeForInnerNode(*this, nullptr) {
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
	EdgeSet bundleEdges;
	if (mapBundleEdges) {
		bundleEdges.init(G);
		m_bundleEdgesForLeaf.init(*this);
	}
	PCNode* stEdgeLeaf = nullptr;
	for (node n : order) {
		if (n == order[order.size() - 1]) {
			// Graph is planar

			if (mapBundleEdges && isTrivial()) {
				// Map all incoming edges of the partner node to the st-edge, unless the partner node is the first node
				if (!twinPoleCandidateInEdges.empty()) {
					m_bundleEdgesForLeaf[stEdgeLeaf].swap(twinPoleCandidateInEdges);
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
					for (edge e : m_bundleEdgesForLeaf[v]) {
						bundleEdges.insert(e);
					}

					if (twinPoleCandidate) {
						twinPoleCandidateInEdges.pushBack(m_incidentEdgeForLeaf[v]);
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
		m_graphNodeForInnerNode[mergedLeaf] = n;

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
					m_graphNodeForInnerNode[mergedLeaf] = n;
				}
				changeNodeType(mergedLeaf, PCNodeType::PNode);
			}
			insertLeaves(outEdges.size(), mergedLeaf, &addedLeaves);

			for (size_t i = 0; i < outEdges.size(); i++) {
				// For each edge store which leaf represents it.
				edge e = outEdges[i];
				PCNode* newLeaf = addedLeaves[i];
				leafRepresentation[e] = newLeaf;
				m_incidentEdgeForLeaf[newLeaf] = e;

				if (mapBundleEdges) {
					if (twinPoleCandidate) {
						m_bundleEdgesForLeaf[newLeaf] = {e};
					} else {
						m_bundleEdgesForLeaf[newLeaf] = bundleEdges.elements();
					}
				}
			}
		} else if (outEdges.size() == 1) {
			// If there is only one edge to a node with a higher number,
			// declare the mergedLeaf as the leaf that represents it.
			leafRepresentation[outEdges.front()] = mergedLeaf;
			m_incidentEdgeForLeaf[mergedLeaf] = outEdges.front();

			if (mapBundleEdges) {
				if (twinPoleCandidate) {
					m_bundleEdgesForLeaf[mergedLeaf] = {outEdges.front()};
				} else {
					m_bundleEdgesForLeaf[mergedLeaf] = bundleEdges.elements();
				}
			}
		}

		OGDF_HEAVY_ASSERT(checkValid());
	}
	OGDF_ASSERT(false); // should either return from processing last node or throw non-planar
}

ogdf::node NodePCRotation::getTrivialPartnerPole() const {
	if (m_n->degree() < 3) {
		return nullptr;
	}
	if (!isTrivial()) {
		return nullptr;
	}
	node partner = m_graphNodeForInnerNode[getRootNode()];
	if (partner->degree() < 3) {
		return nullptr;
	} else {
		return partner;
	}
}

ogdf::node NodePCRotation::getTrivialPartnerPole(const Graph& G, node pole) {
	NodePCRotation pc(G, pole, false);
	return pc.getTrivialPartnerPole();
}

bool NodePCRotation::isEqual(const NodePCRotation& pc) const {
	OGDF_HEAVY_ASSERT(checkValid());
	OGDF_HEAVY_ASSERT(pc.checkValid());
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
	if (possibleOrders<int>() != pc.possibleOrders<int>()) {
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

std::function<void(std::ostream& os, PCNode*, int)> NodePCRotation::uidPrinter() const {
	return [&](std::ostream& os, PCNode* n, int i) -> void {
		if (n->isLeaf()) {
			os << getIncidentEdgeForLeaf(n)->index();
		} else if (getGraphNodeForInnerNode(n) != nullptr) {
			os << getGraphNodeForInnerNode(n)->index();
		} else {
			os << n->getNodeType();
		}
	};
}

std::function<bool(PCNode*, PCNode*)> NodePCRotation::uidComparer() const {
	return [&](PCNode* a, PCNode* b) -> bool {
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
