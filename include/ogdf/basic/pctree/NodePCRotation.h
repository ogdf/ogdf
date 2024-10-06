/** \file
 * \brief An embedding tree representing, for a single node, all orders of incident edges in all planar embeddings.
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

#pragma once

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/pctree/PCEnum.h>
#include <ogdf/basic/pctree/PCNode.h>
#include <ogdf/basic/pctree/PCRegistry.h>
#include <ogdf/basic/pctree/PCTree.h>
#include <ogdf/basic/pctree/util/IntrusiveList.h>

#include <exception>
#include <functional>
#include <iosfwd>

namespace ogdf::pc_tree {
/**
 * This class represents the embedding tree of a single node in a biconnected component.
 * The embedding tree is a PCTree with one leaf for each incident edge that describes all possible rotations of the
 * incident edges in all planar embeddings of the component.
 * Uses the Booth&Lueker Planarity Test for computing the embedding tree, see https://doi.org/10.15475/cpatp.2024 (Section 9.2).
 */
class OGDF_EXPORT NodePCRotation : public PCTree {
protected:
	const Graph* m_G;
	node m_n;

	PCTreeNodeArray<edge> m_incidentEdgeForLeaf;
	PCTreeNodeArray<node> m_graphNodeForInnerNode;
	PCTreeNodeArray<List<edge>> m_bundleEdgesForLeaf;

	NodePCRotation() : m_G(nullptr), m_n(nullptr) { }

public:
	/**
	 * Calculate the embedding tree of node \p n in graph \p G. The connected component of \p n has to be biconnected.
	 * @param G The graph of n.
	 * @param n The node for which we want an embedding tree.
	 * @param mapBundleEdges By default, the data needed for method getPartnerEdgesForLeaf() will also be generated.
	 * 	Set to \c false to turn this off.
	 * @throws GraphNotPlanarException if the component of \p n is non-planar.
	 */
	explicit NodePCRotation(const Graph& G, node n, bool mapBundleEdges = true);

	static node getTrivialPartnerPole(const Graph& G, node n);

	/**
	 * If this embedding tree is trivial (i.e., consists of a single P-node allowing arbitrary rotations), the graphs
	 * must contain another node called pole that mirrors the rotation of this node.
	 * @return the pole corresponding to the node that this embedding tree represents, or \c nullptr if no pole exists
	 */
	node getTrivialPartnerPole() const;

	void generateLeafForIncidentEdgeMapping(EdgeArray<PCNode*>& mapping) const {
		for (auto leaf : getLeaves()) {
			mapping[getIncidentEdgeForLeaf(leaf)] = leaf;
		}
	}

	/**
	 * Returns which incident edge corresponds to which leaf.
	 * @param n A leaf of this PCtree.
	 * @return An edge incident to getNode() corresponding to \p n.
	 */
	edge getIncidentEdgeForLeaf(PCNode* n) const { return m_incidentEdgeForLeaf[n]; }

	/**
	 * This is only needed so that SyncPlan::propagatePQ can fix its mapping when propagating into an adjacent cut.
	 * This method is thus considered internal.
	 */
	void setIncidentEdgeForLeaf(PCNode* n, edge e) { m_incidentEdgeForLeaf[n] = e; }

	void generateInnerNodeForGraphNodeMapping(NodeArray<PCNode*>& mapping) const {
		for (auto leaf : getLeaves()) {
			mapping[getGraphNodeForInnerNode(leaf)] = leaf;
		}
	}

	/**
	 * Returns which graph node created a P-node during the run of the planarity test that resulted in this embedding tree.
	 * @param n A P-node of this PCtree.
	 * @return A node of \p getGraph() corresponding to \p n.
	 */
	node getGraphNodeForInnerNode(PCNode* n) const { return m_graphNodeForInnerNode[n]; }

	void generatePartnerEdgesForIncidentEdge(EdgeArray<const List<edge>*>& mapping) const {
		for (auto leaf : getLeaves()) {
			mapping[getIncidentEdgeForLeaf(leaf)] = &getPartnerEdgesForLeaf(leaf);
		}
	}

	/**
	 * If this embedding tree is trivial and getNode() thus has a partner pole, returns the set of edges incident to the
	 * pole that correspond to a given leaf \p l. If the pole is also trivial, this is a one-to-one mapping. Otherwise,
	 * multiple edges incident to getTrivialPartnerPole() may correspond to one edge incident to getNode(). Setting
	 * \c mapBundleEdges to \c false disables this functionality.
	 * @param l A leaf of this PC-tree.
	 * @return A list of edges incident to getTrivialPartnerPole() corresponding to \p l.
	 */
	const List<edge>& getPartnerEdgesForLeaf(PCNode* l) const { return m_bundleEdgesForLeaf[l]; }

	/**
	 * @return the value passed to \c mapBundleEdges in the constructor
	 */
	bool knowsPartnerEdges() const {
		if (m_bundleEdgesForLeaf.registeredAt() == nullptr) {
			return false;
		} else {
			return &m_bundleEdgesForLeaf.registeredAt()->getForest() == getForest();
		}
	}

	const Graph* getGraph() const { return m_G; }

	node getNode() const { return m_n; }

	/**
	 * Equality check where to leaves are considered the same if they map to the same edge via getIncidentEdgeForLeaf() in both trees.
	 * @param pc The other tree.
	 * @return \c true iff the leaves of both tree map to the same edges and they represent the same restrictions on these leaves/edges.
	 */
	bool isEqual(const NodePCRotation& pc) const;

	/**
	 * @return A function that prints the ID of the respective graph node / edge for an inner node / leaf.
	 * @sa PCTree::uniqueID()
	 */
	std::function<void(std::ostream& os, PCNode*, int)> uidPrinter() const;

	/**
	 * @return A comparison function that uses the same info as uidPrinter() to make the comparison deterministic if the underlying graph structure is the same.
	 * @sa PCTree::uniqueID()
	 */
	std::function<bool(PCNode*, PCNode*)> uidComparer() const;
};

struct OGDF_EXPORT GraphNotPlanarException : public std::exception {
public:
	[[nodiscard]] const char* what() const noexcept override { return "Graph is not planar"; }
};
}
