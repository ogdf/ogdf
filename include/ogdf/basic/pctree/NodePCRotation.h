#pragma once

#include <ogdf/basic/RegisteredSet.h>
#include <ogdf/basic/pctree/PCTree.h>

namespace pc_tree {
using namespace ogdf;

class NodePCRotation : public PCTree {
protected:
	const Graph* m_G;
	node m_n;

	PCTreeNodeArray<edge> incidentEdgeForLeaf;
	PCTreeNodeArray<node> graphNodeForInnerNode;
	PCTreeNodeArray<List<edge>> bundleEdgesForLeaf;

	NodePCRotation() : m_G(nullptr), m_n(nullptr) { }

public:
	explicit NodePCRotation(const Graph& G, node n, bool mapBundleEdges = true);

	static node getTrivialPartnerPole(const Graph& G, node n);

	node getTrivialPartnerPole() const;

	void generateLeafForIncidentEdgeMapping(EdgeArray<PCNode*>& mapping) const {
		for (auto leaf : getLeaves()) {
			mapping[getIncidentEdgeForLeaf(leaf)] = leaf;
		}
	}

	edge getIncidentEdgeForLeaf(PCNode* n) const { return incidentEdgeForLeaf[n]; }

	/*
         * This is needed so that PQPlanarity::propagatePQ can fix its mapping when propagating into an adjacent cut.
         * TODO replace by better interface
         */
	void setIncidentEdgeForLeaf(PCNode* n, edge e) { incidentEdgeForLeaf[n] = e; }

	void generateInnerNodeForGraphNodeMapping(NodeArray<PCNode*>& mapping) const {
		for (auto leaf : getLeaves()) {
			mapping[getGraphNodeForInnerNode(leaf)] = leaf;
		}
	}

	node getGraphNodeForInnerNode(PCNode* n) const { return graphNodeForInnerNode[n]; }

	void generatePartnerEdgesForIncidentEdge(EdgeArray<const List<edge>*>& mapping) const {
		for (auto leaf : getLeaves()) {
			mapping[getIncidentEdgeForLeaf(leaf)] = &getPartnerEdgesForLeaf(leaf);
		}
	}

	const List<edge>& getPartnerEdgesForLeaf(PCNode* l) const { return bundleEdgesForLeaf[l]; }

	bool knowsPartnerEdges() const {
		if (bundleEdgesForLeaf.registeredAt() == nullptr) {
			return false;
		} else {
			return &bundleEdgesForLeaf.registeredAt()->getForest() == getForest();
		}
	}

	const Graph* getGraph() const { return m_G; }

	node getNode() const { return m_n; }

	bool isEqual(const NodePCRotation& pc) const;

	std::function<void(std::ostream& os, PCNode*, int)> uidPrinter() const;

	std::function<bool(PCNode*, PCNode*)> uidComparer() const;
};

struct GraphNotPlanarException : public std::exception {
public:
	[[nodiscard]] const char* what() const noexcept override { return "Graph is not planar"; }
};
}