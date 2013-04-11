/*
 * $Revision: 3368 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-04 20:07:31 +0200 (Do, 04. Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration of class FixEdgeInserterCore and FixEdgeInserterUMLCore,
 * which are the implementation classes for edge insertion with fixed embedding.
 *
 * \author Carsten Gutwenger
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.txt in the root directory of the OGDF installation for details.
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
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/


#ifdef _MSC_VER
#pragma once
#endif

#ifndef OGDF_FIX_EDGE_INS_CORE_H
#define OGDF_FIX_EDGE_INS_CORE_H


#include <ogdf/basic/Timeouter.h>
#include <ogdf/basic/Module.h>
#include <ogdf/planarity/PlanRepLight.h>
#include <ogdf/planarity/RemoveReinsertType.h>


namespace ogdf {

	class FaceSetSimple;
	template<class E> class QueuePure;


	class OGDF_EXPORT FixEdgeInserterCore : public Timeouter
	{
	public:
		FixEdgeInserterCore(
			PlanRepLight &pr,
			const EdgeArray<int>      *pCostOrig,
			const EdgeArray<bool>     *pForbiddenOrig,
			const EdgeArray<__uint32> *pEdgeSubgraphs)
			: m_pr(pr), m_pCost(pCostOrig), m_pForbidden(pForbiddenOrig), m_pSubgraph(pEdgeSubgraphs) { }

		virtual ~FixEdgeInserterCore() { }

		Module::ReturnType call(
			const Array<edge> &origEdges,
			bool keepEmbedding,
			RemoveReinsertType rrPost,
			double percentMostCrossed);

		int runsPostprocessing() const { return m_runsPostprocessing; }

	protected:
		int getCost(edge e, int stSubGraph) const;
		void findShortestPath(const CombinatorialEmbedding &E, edge eOrig, SList<adjEntry> &crossed);
		void findWeightedShortestPath(const CombinatorialEmbedding &E, edge eOrig, SList<adjEntry> &crossed);

		int costCrossed(edge eOrig) const;
		void insertEdge(CombinatorialEmbedding &E, edge eOrig, const SList<adjEntry> &crossed);
		void removeEdge(CombinatorialEmbedding &E, edge eOrig);

		virtual void storeTypeOfCurrentEdge(edge eOrig) { }
		virtual void init(CombinatorialEmbedding &E);
		virtual void cleanup();
		virtual void constructDual(const CombinatorialEmbedding &E);

		virtual void appendCandidates(QueuePure<edge> &queue, node v);
		virtual void appendCandidates(Array<SListPure<edge> > &nodesAtDist, EdgeArray<int> &costDual, int maxCost, node v, int currentDist);

		virtual void insertEdgesIntoDual(const CombinatorialEmbedding &E, adjEntry adjSrc);
		virtual void insertEdgesIntoDualAfterRemove(const CombinatorialEmbedding &E, face f);

		PlanRepLight	&m_pr;

		const EdgeArray<int>		*m_pCost;
		const EdgeArray<bool>		*m_pForbidden;
		const EdgeArray<__uint32>	*m_pSubgraph;

		Graph m_dual;   //!< (Extended) dual graph, constructed/destructed during call.

		EdgeArray<adjEntry> m_primalAdj;   //!< Adjacency entry in primal graph corresponding to edge in dual.
		FaceArray<node>     m_nodeOf;      //!< The node in dual corresponding to face in primal.

		FaceSetSimple       *m_delFaces;
		FaceSetPure         *m_newFaces;

		node m_vS; //!< The node in extended dual representing s.
		node m_vT; //!< The node in extended dual representing t.

		int m_runsPostprocessing; //!< Runs of remove-reinsert method.
	};


	class FixEdgeInserterUMLCore : public FixEdgeInserterCore
	{
	public:
		FixEdgeInserterUMLCore(
			PlanRepLight &pr,
			const EdgeArray<int>      *pCostOrig,
			const EdgeArray<__uint32> *pEdgeSubgraph) : FixEdgeInserterCore(pr, pCostOrig, 0, pEdgeSubgraph) { }

	protected:
		void storeTypeOfCurrentEdge(edge eOrig) { m_typeOfCurrentEdge = m_pr.typeOrig(eOrig); }
		void init(CombinatorialEmbedding &E);
		void cleanup();
		void constructDual(const CombinatorialEmbedding &E);

		void appendCandidates(QueuePure<edge> &queue, node v);
		void appendCandidates(Array<SListPure<edge> > &nodesAtDist, EdgeArray<int> &costDual, int maxCost, node v, int currentDist);

		void insertEdgesIntoDual(const CombinatorialEmbedding &E, adjEntry adjSrc);
		void insertEdgesIntoDualAfterRemove(const CombinatorialEmbedding &E, face f);

		EdgeArray<bool> m_primalIsGen; //!< true iff corresponding primal edge is a generalization.
		Graph::EdgeType	m_typeOfCurrentEdge;
	};

}

#endif
