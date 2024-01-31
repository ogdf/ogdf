/** \file
 * \brief Implementation of the auxiliary graph as well as the edges
 * and nodes of it for the Blossom V algorithm.
 *
 * \author Joshua Sangmeister
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

#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/NodeArray.h>
#include <ogdf/graphalg/matching_blossom/AlternatingTree.h>
#include <ogdf/graphalg/matching_blossom/BlossomVHelper.h>
#include <ogdf/graphalg/matching_blossom/Cycle.h>
#include <ogdf/graphalg/matching_blossom/PQ.h>
#include <ogdf/graphalg/matching_blossom/Pseudonode.h>
#include <ogdf/graphalg/matching_blossom/utils.h>

#include <unordered_set>
#include <vector>

namespace ogdf {
namespace matching_blossom {

template<class TWeight>
class BlossomVHelper;
template<class TWeight>
class AuxEdge;

template<class TWeight>
class AuxNode {
	using EdgePQ = BlossomPQ<edge, TWeight>;
	using NodePQ = BlossomPQ<node, TWeight>;

	//! The actual node in the auxiliary graph.
	node m_node;

	//! The auxiliary graph this node belongs to.
	BlossomVHelper<TWeight>& m_helper;

	//! The alternating tree this node represents.
	AlternatingTree<TWeight> m_tree;

	//! The cummulated delta of this tree.
	double m_delta = 0;

	//! Eges between even nodes in this tree.
	EdgePQ m_evenEvenEdges;
	//! Edges between an even node in this tree and a free node.
	EdgePQ m_evenFreeEdges;
	//! All odd pseudonodes in this tree.
	NodePQ m_oddPseudonodes;

	//! Structure representing the current edge of this tree. To avoid resetting all edges of all
	//! trees after every iteration, we save the iteration number to see if the edge pointer is invalid.
	struct {
		AuxEdge<TWeight>* edge = nullptr;
		long iteration = -1;
	} m_currentEdge;

public:
	AuxNode(node auxGraphNode, node graphNode, BlossomVHelper<TWeight>& helper)
		: m_node(auxGraphNode), m_helper(helper), m_tree(m_helper, graphNode) {
		OGDF_ASSERT(auxGraphNode->graphOf() != &m_helper.graph());
		OGDF_ASSERT(graphNode->graphOf() == &m_helper.graph());
	}

	/* Getters & setters */

	//! The aux edge pointing to the current aux node.
	AuxEdge<TWeight>* currentEdge() {
		if (m_currentEdge.iteration == m_helper.currentIteration) {
			return m_currentEdge.edge;
		} else {
			return nullptr;
		}
	};

	//! Sets the current edge of this tree to \p edge and update the iteration to the current one.
	void setCurrentEdge(AuxEdge<TWeight>* edge) {
		m_currentEdge.edge = edge;
		m_currentEdge.iteration = m_helper.currentIteration;
	}

	node graphNode() { return m_node; }

	AlternatingTree<TWeight>& tree() { return m_tree; }

	EdgePQ& evenEvenEdges() { return m_evenEvenEdges; }

	EdgePQ& evenFreeEdges() { return m_evenFreeEdges; }

	NodePQ& oddPseudonodes() { return m_oddPseudonodes; }

	double delta() { return m_delta; }

	/* End of getters & setters */

	//! The delta of this tree for the given node.
	double delta(node v) {
		OGDF_ASSERT(m_tree.contains(v));
		if (m_tree.isEven(v)) {
			return m_delta;
		} else {
			return -m_delta;
		}
	}

	//! Add \p delta to this tree. \p delta must be non-negative.
	void addDelta(double delta) {
		OGDF_ASSERT(delta >= 0);
		m_delta += delta;
	}

	//! Add \p v to the list of odd pseudonodes of this tree.
	void addOddPseudonode(node v) { m_oddPseudonodes.push(v, m_helper.y(v)); }

	//! Add \p e to the list of even-even edges of this tree.
	void addEvenEvenEdge(edge e) { m_evenEvenEdges.push(e, m_helper.getReducedWeight(e)); }

	//! Add \p e to the list of even-free edges of this tree.
	void addEvenFreeEdge(edge e) { m_evenFreeEdges.push(e, m_helper.getReducedWeight(e)); }
};

template<class TWeight>
class AuxEdge {
private:
	using EdgePQ = BlossomPQ<edge, TWeight>;

	// the actual edge in the auxiliary graph
	edge m_edge;

	BlossomVHelper<TWeight>& m_helper;

	// edges between an even node in the source tree and an even node in the target tree
	EdgePQ m_evenEvenEdges;
	// edges between an even node in the source tree and an odd node in the target tree
	EdgePQ m_evenOddEdges;
	// edges between an odd node in the source tree and an even node in the target tree
	EdgePQ m_oddEvenEdges;

	//! Helper function to add an edge to a priority queue with its current reduced weight.
	void addEdgeToPQ(edge e, EdgePQ& pq) {
		OGDF_ASSERT(e->graphOf() == &m_helper.graph());
		pq.push(e, m_helper.getReducedWeight(e));
	}

public:
	AuxEdge(edge e, BlossomVHelper<TWeight>& helper) : m_edge(e), m_helper(helper) {
		OGDF_ASSERT(e->graphOf() != &m_helper.graph());
	}

	/* Getters */

	edge graphEdge() { return m_edge; }

	EdgePQ& evenEvenEdges() { return m_evenEvenEdges; }

	EdgePQ& evenOddEdges() { return m_evenOddEdges; }

	EdgePQ& oddEvenEdges() { return m_oddEvenEdges; }

	/* End getters */

	//! Returns evenOddEdges or oddEvenEdges, depending on the perspective of \p auxNode (the even node).
	EdgePQ& evenOddEdgesFromPerspective(AuxNode<TWeight>* auxNode) {
		OGDF_ASSERT(m_edge->isIncident(auxNode->graphNode()));
		if (m_edge->source() == auxNode->graphNode()) {
			return m_evenOddEdges;
		} else {
			return m_oddEvenEdges;
		}
	}

	//! Adds \p e to evenEvenEdges.
	void addEvenEvenEdge(edge e) { addEdgeToPQ(e, m_evenEvenEdges); }

	//! Adds \p e to evenOddEdges or oddEvenEdges depending on the perspective of \p auxNode.
	void addEvenOddEdgeFromPerspective(edge e, AuxNode<TWeight>* auxNode) {
		addEdgeToPQ(e, evenOddEdgesFromPerspective(auxNode));
	}

	//! Whether or not this edge connects two trees via an alternating equality edge.
	bool hasEvenOddEqualityEdge() {
		return (!m_evenOddEdges.empty() && m_helper.isEqualityEdge(m_evenOddEdges.topElement()))
				|| (!m_oddEvenEdges.empty() && m_helper.isEqualityEdge(m_oddEvenEdges.topElement()));
	}
};

template<class TWeight>
class AuxGraph {
	BlossomVHelper<TWeight>& m_helper;

	// the auxiliary graph
	GraphCopySimple m_graph;

	//! maps auxiliary graph nodes to their AuxNode objects
	NodeArray<AuxNode<TWeight>*> m_auxGraphNodeMap;
	//! maps auxiliary graph edges to their AuxEdge objects
	EdgeArray<AuxEdge<TWeight>*> m_auxGraphEdgeMap;
	//! maps normal graph nodes to the AuxNode whose tree they belong to
	NodeArray<AuxNode<TWeight>*> m_nodeAuxNodeMap;

	//! Creates an AuxNode for \p v and stores it in the appropriate maps.
	AuxNode<TWeight>* createAuxNode(node v) {
		node orig = m_graph.original(v);
		auto auxNode = new AuxNode<TWeight>(v, orig, m_helper);
		m_auxGraphNodeMap[v] = auxNode;
		m_nodeAuxNodeMap[orig] = auxNode;
		return auxNode;
	}

	//! Creates an AuxEdge for \p e and stores it in the appropriate maps.
	AuxEdge<TWeight>* createAuxEdge(edge e) {
		auto auxEdge = new AuxEdge<TWeight>(e, m_helper);
		m_auxGraphEdgeMap[e] = auxEdge;
		return auxEdge;
	}

public:
	AuxGraph(BlossomVHelper<TWeight>& helper)
		: m_helper(helper)
		, m_auxGraphNodeMap(m_graph, nullptr)
		, m_auxGraphEdgeMap(m_graph, nullptr)
		, m_nodeAuxNodeMap(m_helper.graph(), nullptr) {
		reset();
	}

	/* Getters */

	GraphCopySimple& graph() { return m_graph; }

	const NodeArray<AuxNode<TWeight>*>& nodes() const { return m_auxGraphNodeMap; }

	const EdgeArray<AuxEdge<TWeight>*>& edges() const { return m_auxGraphEdgeMap; }

	//! Returns the AuxEdge corresponding to \p e, which must be an edge of the auxiliary graph.
	AuxEdge<TWeight>* auxEdge(edge e) { return m_auxGraphEdgeMap[e]; }

	//! Returns the AuxNode corresponding to \p v, which must be a node of the auxiliary graph.
	AuxNode<TWeight>* auxNode(node v) { return m_auxGraphNodeMap[v]; }

	//! Returns the AuxNode of whose tree the node \p v of the actual graph is a part of.
	AuxNode<TWeight>* treeAuxNode(node v) { return m_nodeAuxNodeMap[v]; }

	/* End of getters */

	//! Rebuilds the auxiliary graph from the current graph.
	void reset() {
		m_graph.clear();
		m_graph.setOriginalGraph(m_helper.graph());

		// create all aux nodes
		for (node v : m_helper.graph().nodes) {
			if (!m_helper.matching(v)) {
				createAuxNode(m_graph.newNode(v));
			}
		}
		// fill all edges
		for (node auxGraphNode : m_graph.nodes) {
			node u = m_graph.original(auxGraphNode);
			for (auto adj : u->adjEntries) {
				edge e = adj->theEdge();
				node v = adj->twinNode();
				if (m_graph.copy(v) != nullptr) {
					if (m_graph.copy(e) == nullptr) {
						auto auxEdge = createAuxEdge(m_graph.newEdge(e));
						auxEdge->addEvenEvenEdge(e);
					}
				} else {
					auxNode(auxGraphNode)->addEvenFreeEdge(e);
				}
			}
		}
	}

	//! Returns the AuxNode of the auxiliary graph whose tree contains at least one endpoint of
	//! \p e. \p e must be part of the actual graph.
	AuxNode<TWeight>* auxNodeForEdge(edge e) {
		for (node v : e->nodes()) {
			if (auto auxNode = treeAuxNode(v)) {
				return auxNode;
			}
		}
		return nullptr;
	}

	//! Returns the tree which contains \p v, or nullptr if \p v is free.
	AlternatingTree<TWeight>* treeOf(node v) {
		auto auxNode = treeAuxNode(v);
		return auxNode ? &auxNode->tree() : nullptr;
	}

	//! Creates an AuxEdge between \p auxNode and \p current if it does not exist and sets the
	//! currentEdge of \p auxNode accordingly.
	AuxEdge<TWeight>* assertCurrentEdge(AuxNode<TWeight>* auxNode, AuxNode<TWeight>* current) {
		if (auxNode->currentEdge() == nullptr) {
			edge e = m_graph.newEdge(auxNode->graphNode(), current->graphNode());
			auto auxEdge = createAuxEdge(e);
			auxNode->setCurrentEdge(auxEdge);
			return auxEdge;
		} else {
			return auxNode->currentEdge();
		}
	}

	//! Sets the AuxNode of \p v to \p auxNode.
	void setAuxNode(node v, AuxNode<TWeight>* auxNode) { m_nodeAuxNodeMap[v] = auxNode; }

	//! Removes \p auxNode and all incident edges from the auxiliary graph.
	void deleteNode(AuxNode<TWeight>* auxNode) {
		auto tree = auxNode->tree();
		for (auto treeNodes : {tree.evenNodes, tree.oddNodes}) {
			for (node v : treeNodes) {
				m_helper.y(v) = m_helper.realY(v);
				m_nodeAuxNodeMap[v] = nullptr;
			}
		}
		for (auto adj : auxNode->graphNode()->adjEntries) {
			AuxEdge<TWeight>* ae = auxEdge(adj->theEdge());
			delete ae;
		}
		m_graph.delNode(auxNode->graphNode());
		delete auxNode;
	}

	//! Calculates the connected components of the auxiliary graph and stores them in \p components.
	//! Only tight even-odd/odd-even edges between different trees are taken into account.
	//!
	//! Note: We cannot use connectedComponents from simple_graph_alg.h since we need to iterate
	//! the components one after another and not all edges are taken into account.
	void connectedComponents(std::vector<std::unordered_set<node>>& components) {
		NodeArray<bool> visited(m_graph, false);
		std::vector<node> stack;
		for (node u : m_graph.nodes) {
			if (!visited[u]) {
				visited[u] = true;
				stack.push_back(u);
				std::unordered_set<node> component = {u};
				while (!stack.empty()) {
					node v = stack.back();
					stack.pop_back();
					for (adjEntry adj : v->adjEntries) {
						node w = adj->twinNode();
						if (!visited[w]) {
							auto auxEdge = this->auxEdge(adj->theEdge());
							if (auxEdge->hasEvenOddEqualityEdge()) {
								component.insert(w);
								visited[w] = true;
								stack.push_back(w);
							}
						}
					}
				}
				components.push_back(component);
			}
		}
#ifdef OGDF_DEBUG
		// check that all nodes are in a component
		int count = 0;
		for (auto component : components) {
			count += component.size();
		}
		OGDF_ASSERT(count == m_graph.numberOfNodes());
#endif
	};
};

}
}
