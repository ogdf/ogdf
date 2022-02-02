/** \file
 * \brief Includes declaration of dual graph class.
 *
 * \author Max Ilsen, Michael Schulz
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
#include <ogdf/basic/NodeArray.h>
#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/FaceArray.h>

namespace ogdf {

template<bool isConst> class DualGraphBase;

using DualGraph = DualGraphBase<true>;

using DynamicDualGraph = DualGraphBase<false>;

/**
 * A dual graph including its combinatorial embedding of an embedded graph.
 *
 * Dual edges are rotated counter-clockwise compared to the primal ones.
 *
 * @ingroup graphs
 */
template<bool isConst>
class OGDF_EXPORT DualGraphBase : public CombinatorialEmbedding
{
public:
	using Embedding = typename std::conditional<isConst,
		  const ConstCombinatorialEmbedding, CombinatorialEmbedding>::type;

	//! Constructor; creates dual graph and its combinatorial embedding
	explicit DualGraphBase(Embedding &CE) : m_primalEmbedding(CE)
	{
		const Graph &primalGraph = CE.getGraph();
		init(*(new Graph));
		Graph &dualGraph = getGraph();

		m_dualNode.init(CE);
		m_dualEdge.init(primalGraph);
		m_dualFace.init(primalGraph);
#ifdef OGDF_DEBUG
		m_primalNode.init(*this, nullptr);
#else
		m_primalNode.init(*this);
#endif
		m_primalFace.init(dualGraph);
		m_primalEdge.init(dualGraph);

		// create dual nodes
		for(face f : CE.faces)
		{
			node vDual = dualGraph.newNode();
			m_dualNode[f] = vDual;
			m_primalFace[vDual] = f;
		}

		// create dual edges
		for(edge e : primalGraph.edges)
		{
			adjEntry aE = e->adjSource();
			node vDualSource = m_dualNode[CE.rightFace(aE)];
			node vDualTarget = m_dualNode[CE.leftFace(aE)];
			edge eDual = dualGraph.newEdge(vDualSource, vDualTarget);
			m_primalEdge[eDual] = e;
			m_dualEdge[e] = eDual;
		}

		// sort adjElements of every dual node corresponding to dual embedding
		for(face f : CE.faces)
		{
			node vDual = m_dualNode[f];
			List<adjEntry> newOrder;

			for(adjEntry adj : f->entries) {
				edge e = adj->theEdge();
				edge eDual = m_dualEdge[e];
				bool isSource = adj == e->adjSource();
				adjEntry adjDual = isSource ?
					eDual->adjSource() : eDual->adjTarget();
				newOrder.pushBack(adjDual);
			}

			dualGraph.sort(vDual, newOrder);
		}

		// calculate dual faces and links to corresponding primal nodes
		computeFaces();

		for(node v : primalGraph.nodes)
		{
			edge ePrimal = v->firstAdj()->theEdge();
			edge eDual = m_dualEdge[ePrimal];
			face fDual = rightFace(eDual->adjSource());
			if(ePrimal->source()==v)
				fDual = leftFace(eDual->adjSource());

			OGDF_ASSERT(m_primalNode[fDual] == nullptr);

			m_dualFace[v] = fDual;
			m_primalNode[fDual] = v;
		}
	}

	//! Destructor
	~DualGraphBase<isConst>()
	{
		clear();
		delete m_cpGraph;
	}

	//! Returns a reference to the combinatorial embedding of the primal graph
	Embedding &getPrimalEmbedding() const { return m_primalEmbedding; }

	//! Returns a reference to the primal graph
	const Graph &getPrimalGraph() const { return m_primalEmbedding.getGraph(); }

	//! @name Lookup functions
	//! @{

	//! Returns the node in the primal graph corresponding to \p f.
	/**
	* @param f is a face in the embedding of the dual graph
	* \return the corresponding node in the primal graph
	*/
	const node &primalNode(face f) const { return m_primalNode[f]; }

	//! Returns the edge in the primal graph corresponding to \p e.
	/**
	* @param e is an edge in the dual graph
	* \return the corresponding edge in the primal graph
	*/
	const edge &primalEdge(edge e) const { return m_primalEdge[e]; }

	//! Returns the face in the embedding of the primal graph corresponding to \p v.
	/**
	* @param v is a node in the dual graph
	* \return the corresponding face in the embedding of the primal graph
	*/
	const face &primalFace(node v) const { return m_primalFace[v]; }

	//! Returns the node in the dual graph corresponding to \p f.
	/**
	* @param f is a face in the embedding of the primal graph
	* \return the corresponding node in the dual graph
	*/
	const node &dualNode(face f) const { return m_dualNode[f]; }

	//! Returns the edge in the dual graph corresponding to \p e.
	/**
	* @param e is an edge in the primal graph
	* \return the corresponding edge in the dual graph
	*/
	const edge &dualEdge(edge e) const { return m_dualEdge[e]; }

	//! Returns the face in the embedding of the dual graph corresponding to \p v.
	/**
	* @param v is a node in the primal graph
	* \return the corresponding face in the embedding of the dual graph
	*/
	const face &dualFace(node v) const { return m_dualFace[v]; }

	//! @}
	//! @name Updating the dual graph (also updates primal embedding)
	//! @{

	//! @copydoc CombinatorialEmbedding::split(edge)
	template<bool isConstSFINAE = isConst, typename std::enable_if<!isConstSFINAE, int>::type = 0>
	edge splitPrimal(edge e) {
		edge oldDualEdge {m_dualEdge[e]};

#ifdef OGDF_DEBUG
		node oldPrimalNode {e->target()};
		face oldDualFace {rightFace(oldDualEdge->adjSource())};
#endif

		// Create new edge in the primal graph.
		edge newPrimalEdge {m_primalEmbedding.split(e)};
		node newPrimalNode {newPrimalEdge->source()};

		// Create new edge in the dual graph.
		edge newDualEdge {CombinatorialEmbedding::splitFace(
				oldDualEdge->adjSource(),
				oldDualEdge->adjSource()->faceCycleSucc())};
		face newDualFace {leftFace(newDualEdge->adjSource())};

		// Update node and edge mappings.
		m_dualEdge[newPrimalEdge] = newDualEdge;
		m_primalEdge[newDualEdge] = newPrimalEdge;

		m_dualFace[newPrimalNode] = newDualFace;
		m_primalNode[newDualFace] = newPrimalNode;

		OGDF_ASSERT(m_dualFace[oldPrimalNode] == oldDualFace);
		OGDF_ASSERT(m_primalNode[oldDualFace] == oldPrimalNode);

#ifdef OGDF_HEAVY_DEBUG
		consistencyCheck();
#endif

		return newPrimalEdge;
	}

	//! @copydoc CombinatorialEmbedding::unsplit(edge, edge)
	template<bool isConstSFINAE = isConst, typename std::enable_if<!isConstSFINAE, int>::type = 0>
	void unsplitPrimal(edge eIn, edge eOut) {
		// Join faces in the dual graph, unsplit edge in the primal graph.
		face oldDualFace {CombinatorialEmbedding::joinFaces(m_dualEdge[eOut])};
		m_primalEmbedding.unsplit(eIn, eOut);
		node oldPrimalNode {eIn->target()};

		// Update node and edge mappings.
		m_dualFace[oldPrimalNode] = oldDualFace;
		m_primalNode[oldDualFace] = oldPrimalNode;
		OGDF_ASSERT(m_primalEdge[m_dualEdge[eIn]] == eIn);

#ifdef OGDF_HEAVY_DEBUG
		consistencyCheck();
#endif
	}

	//! @copydoc CombinatorialEmbedding::splitNode(adjEntry, adjEntry)
	template<bool isConstSFINAE = isConst, typename std::enable_if<!isConstSFINAE, int>::type = 0>
	node splitNodePrimal(adjEntry adjStartLeft, adjEntry adjStartRight) {
		node oldPrimalNode {adjStartLeft->theNode()};
		face oldDualFace {m_dualFace[oldPrimalNode]};

		// Split node in the primal graph.
		node newPrimalNode {m_primalEmbedding.splitNode(adjStartLeft, adjStartRight)};
		edge newPrimalEdge {adjStartLeft->cyclicPred()->theEdge()};

		// Split face in the dual graph.
		adjEntry dualAdjLeft {dualAdj(adjStartLeft, true)};
		adjEntry dualAdjRight {dualAdj(adjStartRight, true)};
		OGDF_ASSERT(dualAdjLeft->theNode() ==
			m_dualNode[m_primalEmbedding.leftFace(adjStartLeft)]);
		OGDF_ASSERT(dualAdjRight->theNode() ==
			m_dualNode[m_primalEmbedding.leftFace(adjStartRight)]);
		edge newDualEdge {CombinatorialEmbedding::splitFace(dualAdjLeft, dualAdjRight)};
		face newDualFace {leftFace(newDualEdge->adjSource())};

		// Update node and edge mappings.
		m_dualEdge[newPrimalEdge] = newDualEdge;
		m_primalEdge[newDualEdge] = newPrimalEdge;

		m_dualFace[newPrimalNode] = oldDualFace;
		m_primalNode[oldDualFace] = newPrimalNode;

		m_dualFace[oldPrimalNode] = newDualFace;
		m_primalNode[newDualFace] = oldPrimalNode;

#ifdef OGDF_HEAVY_DEBUG
		consistencyCheck();
#endif

		return newPrimalNode;
	}

	//! @copydoc CombinatorialEmbedding::contract(edge, bool)
	template<bool isConstSFINAE = isConst, typename std::enable_if<!isConstSFINAE, int>::type = 0>
	node contractPrimal(edge e, bool keepSelfLoops = false) {
		edge dualEdge {m_dualEdge[e]};

		// Contract e in in the primal graph, join faces in the dual graph.
		// TODO Joining faces in the dual graph currently always acts as if
		// keepSelfLoops is true.
		node newPrimalNode {m_primalEmbedding.contract(e, keepSelfLoops)};
		face newDualFace {CombinatorialEmbedding::joinFaces(dualEdge)};

		// Update node and edge mappings.
		m_dualFace[newPrimalNode] = newDualFace;
		m_primalNode[newDualFace] = newPrimalNode;

#ifdef OGDF_HEAVY_DEBUG
		consistencyCheck();
#endif

		return newPrimalNode;
	}

	//! @copydoc CombinatorialEmbedding::splitFace(adjEntry, adjEntry, bool)
	template<bool isConstSFINAE = isConst, typename std::enable_if<!isConstSFINAE, int>::type = 0>
	edge splitFacePrimal(adjEntry adjSrc, adjEntry adjTgt, bool sourceAfter = false) {
#ifdef OGDF_DEBUG
		face oldPrimalFace {m_primalEmbedding.rightFace(adjSrc)};
		node oldDualNode {m_dualNode[oldPrimalFace]};
#endif

		// Create new edge in the primal graph.
		edge newPrimalEdge {m_primalEmbedding.splitFace(adjSrc, adjTgt, sourceAfter)};
		face newPrimalFace {m_primalEmbedding.leftFace(newPrimalEdge->adjSource())};

		// Create new edge in the dual graph.
		adjEntry leftAdj {dualAdj(adjTgt)};
		adjEntry rightAdj {dualAdj(adjSrc)};
		OGDF_ASSERT(leftAdj->theNode() == oldDualNode);
		OGDF_ASSERT(rightAdj->theNode() == oldDualNode);

		node newDualNode {CombinatorialEmbedding::splitNode(leftAdj, rightAdj)};
		edge newDualEdge {leftAdj->cyclicPred()->theEdge()};

		// Update node and edge mappings.
		m_dualEdge[newPrimalEdge] = newDualEdge;
		m_primalEdge[newDualEdge] = newPrimalEdge;

		m_dualNode[newPrimalFace] = newDualNode;
		m_primalFace[newDualNode] = newPrimalFace;

		OGDF_ASSERT(m_dualNode[oldPrimalFace] == oldDualNode);
		OGDF_ASSERT(m_primalFace[oldDualNode] == oldPrimalFace);

#ifdef OGDF_HEAVY_DEBUG
		consistencyCheck();
#endif

		return newPrimalEdge;
	}

	//! @copydoc CombinatorialEmbedding::addEdgeToIsolatedNode(node, adjEntry)
	template<bool isConstSFINAE = isConst, typename std::enable_if<!isConstSFINAE, int>::type = 0>
	edge addEdgeToIsolatedNodePrimal(node v, adjEntry adjTgt) {
		return addEdgeToIsolatedNodePrimal(adjTgt, v, false);
	}

	//! @copydoc CombinatorialEmbedding::addEdgeToIsolatedNode(adjEntry, node)
	template<bool isConstSFINAE = isConst, typename std::enable_if<!isConstSFINAE, int>::type = 0>
	edge addEdgeToIsolatedNodePrimal(adjEntry adjSrc, node v) {
		return addEdgeToIsolatedNodePrimal(adjSrc, v, true);
	}

	//! @copydoc CombinatorialEmbedding::joinFaces(edge)
	template<bool isConstSFINAE = isConst, typename std::enable_if<!isConstSFINAE, int>::type = 0>
	face joinFacesPrimal(edge e) {
		edge eDual {m_dualEdge[e]};

		// Join faces in the primal graph, contract edges in the dual graph.
		face oldPrimalFace {m_primalEmbedding.joinFaces(e)};
		node oldDualNode {CombinatorialEmbedding::contract(eDual, true)};

		m_primalFace[oldDualNode] = oldPrimalFace;
		m_dualNode[oldPrimalFace] = oldDualNode;

#ifdef OGDF_HEAVY_DEBUG
		consistencyCheck();
#endif

		return oldPrimalFace;
	}

	//! @copydoc CombinatorialEmbedding::removeDeg1(node)
	template<bool isConstSFINAE = isConst, typename std::enable_if<!isConstSFINAE, int>::type = 0>
	void removeDeg1Primal(node v) {
		OGDF_ASSERT(v->degree() == 1);

		edge primalEdge {v->firstAdj()->theEdge()};
		edge dualEdge {m_dualEdge[primalEdge]};
#ifdef OGDF_DEBUG
		node otherPrimalNode {primalEdge->opposite(v)};
#endif

		m_primalEmbedding.removeDeg1(v);
#ifdef OGDF_DEBUG
		face newDualFace =
#endif
			CombinatorialEmbedding::joinFaces(dualEdge);

		OGDF_ASSERT(m_dualFace[otherPrimalNode] == newDualFace);
		OGDF_ASSERT(m_primalNode[newDualFace] == otherPrimalNode);

#ifdef OGDF_HEAVY_DEBUG
		consistencyCheck();
#endif
	}

	//! @copydoc CombinatorialEmbedding::reverseEdge(edge)
	template<bool isConstSFINAE = isConst, typename std::enable_if<!isConstSFINAE, int>::type = 0>
	void reverseEdgePrimal(edge e) {
		m_primalEmbedding.reverseEdge(e);
		CombinatorialEmbedding::reverseEdge(m_dualEdge[e]);

#ifdef OGDF_HEAVY_DEBUG
		consistencyCheck();
#endif
	}

#ifdef OGDF_DEBUG
	//! Asserts that this embedding is consistent.
	void consistencyCheck() const {
		Embedding::consistencyCheck();

		const Graph &primalGraph {m_primalEmbedding.getGraph()};
		const Graph &dualGraph {getGraph()};
		OGDF_ASSERT(primalGraph.numberOfNodes() == numberOfFaces());
		OGDF_ASSERT(primalGraph.numberOfEdges() == dualGraph.numberOfEdges());
		OGDF_ASSERT(m_primalEmbedding.numberOfFaces() == dualGraph.numberOfNodes());

		for (node vDual : dualGraph.nodes) {
			OGDF_ASSERT(m_dualNode[m_primalFace[vDual]] == vDual);
		}
		for (edge eDual : dualGraph.edges) {
			OGDF_ASSERT(m_dualEdge[m_primalEdge[eDual]] == eDual);

			// A dual edge leads from the right to the left face of its primal edge.
			OGDF_ASSERT(leftFace(eDual->adjSource()) ==
				m_dualFace[m_primalEdge[eDual]->source()]);
			OGDF_ASSERT(rightFace(eDual->adjSource()) ==
				m_dualFace[m_primalEdge[eDual]->target()]);
		}
		for (face fDual : Embedding::faces) {
			OGDF_ASSERT(m_dualFace[m_primalNode[fDual]] == fDual);
		}
	}
#endif

protected:
	Embedding &m_primalEmbedding; //!< The embedding of the primal graph.
	FaceArray<node> m_primalNode; //!< The corresponding node in the primal graph.
	NodeArray<face> m_primalFace; //!< The corresponding facee in the embedding of the primal graph.
	EdgeArray<edge> m_primalEdge; //!< The corresponding edge in the primal graph.
	FaceArray<node> m_dualNode; //!< The corresponding node in the dual graph.
	NodeArray<face> m_dualFace; //!< The corresponding face in embedding of the dual graph.
	EdgeArray<edge> m_dualEdge; //!< The corresponding edge in the dual graph.

private:
	//! @copydoc CombinatorialEmbedding::addEdgeToIsolatedNode(adjEntry, node, bool)
	template<bool isConstSFINAE = isConst, typename std::enable_if<!isConstSFINAE, int>::type = 0>
	edge addEdgeToIsolatedNodePrimal(adjEntry adj, node v, bool adjSrc) {
#ifdef OGDF_DEBUG
		face oldPrimalFace {m_primalEmbedding.rightFace(adj)};
		node oldDualNode {m_dualNode[oldPrimalFace]};
#endif
		node oldPrimalNode {adj->theNode()};
		face oldDualFace {m_dualFace[oldPrimalNode]};

		adjEntry primalAdj {adj->faceCyclePred()};
		edge newPrimalEdge {adjSrc ?
			m_primalEmbedding.addEdgeToIsolatedNode(adj, v) :
			m_primalEmbedding.addEdgeToIsolatedNode(v, adj)
		};
		// If the new primal edge goes from v to adj, the new dual self-loop has to
		// go from its right to its left side. Hence, we pass !adjSrc to splitFace().
		adjEntry dualAdjEntry {dualAdj(primalAdj)};
		OGDF_ASSERT(dualAdjEntry->theNode() == oldDualNode);
		edge newDualEdge {CombinatorialEmbedding::splitFace(dualAdjEntry, dualAdjEntry, !adjSrc)};
		face newDualFace {leftFace(newDualEdge->adjSource())};

		// Update node and edge mappings.
		m_dualEdge[newPrimalEdge] = newDualEdge;
		m_primalEdge[newDualEdge] = newPrimalEdge;

		if (adjSrc) {
			m_primalNode[oldDualFace] = v;
			m_dualFace[v] = oldDualFace;

			m_primalNode[newDualFace] = oldPrimalNode;
			m_dualFace[oldPrimalNode] = newDualFace;
		} else {
			m_primalNode[newDualFace] = v;
			m_dualFace[v] = newDualFace;

			OGDF_ASSERT(m_primalNode[oldDualFace] == oldPrimalNode);
			OGDF_ASSERT(m_dualFace[oldPrimalNode] == oldDualFace);
		}

#ifdef OGDF_HEAVY_DEBUG
		consistencyCheck();
#endif

		return newPrimalEdge;
	}

	//! Returns the corresponding adjEntry of the dual edge of \p primalAdj
	//! (or the opposite adjEntry of the dual edge if \p reverse is set).
	inline adjEntry dualAdj(adjEntry primalAdj, bool reverse = false) {
		return primalAdj->isSource() != reverse ?
			m_dualEdge[primalAdj->theEdge()]->adjSource() :
			m_dualEdge[primalAdj->theEdge()]->adjTarget();
	}
};

}
