/*
 * $Revision: 3368 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-04 20:07:31 +0200 (Do, 04. Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration of class VarEdgeInserterCore and VarEdgeInserterUMLCore,
 * which are the implementation classes for edge insertion with variable embedding.
 *
 * \author Carsten Gutwenger<br>Jan Papenfu&szlig;
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

#ifndef OGDF_VAR_EDGE_INS_DYN_CORE_H
#define OGDF_VAR_EDGE_INS_DYN_CORE_H


#include <ogdf/basic/Timeouter.h>
#include <ogdf/basic/Module.h>
#include <ogdf/planarity/PlanRepLight.h>
#include <ogdf/planarity/RemoveReinsertType.h>


namespace ogdf {


	class OGDF_EXPORT VarEdgeInserterDynCore : public Timeouter
	{
	public:
		VarEdgeInserterDynCore(
			PlanRepLight &pr,
			const EdgeArray<int>      *pCostOrig,
			const EdgeArray<bool>     *pForbiddenOrig,
			const EdgeArray<__uint32> *pEdgeSubgraphs)
			: m_pr(pr), m_pCost(pCostOrig), m_pForbidden(pForbiddenOrig), m_pSubgraph(pEdgeSubgraphs) { }

		virtual ~VarEdgeInserterDynCore() { }

		Module::ReturnType call(
			const Array<edge> &origEdges,
			RemoveReinsertType rrPost,
			double percentMostCrossed);

		int runsPostprocessing() const { return m_runsPostprocessing; }

	protected:
		class BCandSPQRtrees;
		class ExpandedGraph;

		int costCrossed(edge eOrig) const;

		void insert(edge eOrig, SList<adjEntry> &eip);
		void blockInsert(node s, node t, List<adjEntry> &L);

		virtual void storeTypeOfCurrentEdge(edge eOrig) { }
		virtual BCandSPQRtrees *createBCandSPQRtrees();
		virtual ExpandedGraph *createExpandedGraph(BCandSPQRtrees &BC);

		virtual void buildSubpath(node v,
			node vPred,
			node vSucc,
			List<adjEntry> &L,
			ExpandedGraph &Exp,
			node s,
			node t);

		static const int c_bigM = 10000;
		PlanRepLight	&m_pr;

		const EdgeArray<int>		*m_pCost;
		const EdgeArray<bool>		*m_pForbidden;
		const EdgeArray<__uint32>	*m_pSubgraph;

		BCandSPQRtrees *m_pBC;

		int m_runsPostprocessing; //!< Runs of remove-reinsert method.
	};


	class VarEdgeInserterDynUMLCore : public VarEdgeInserterDynCore
	{
	public:
		VarEdgeInserterDynUMLCore(
			PlanRepLight &pr,
			const EdgeArray<int>      *pCostOrig,
			const EdgeArray<__uint32> *pEdgeSubgraph) : VarEdgeInserterDynCore(pr, pCostOrig, 0, pEdgeSubgraph) { }

	protected:
		class BCandSPQRtreesUML;
		class ExpandedGraphUML;

		void storeTypeOfCurrentEdge(edge eOrig) { m_typeOfCurrentEdge = m_pr.typeOrig(eOrig); }
		virtual BCandSPQRtrees *createBCandSPQRtrees();
		virtual ExpandedGraph *createExpandedGraph(BCandSPQRtrees &BC);
		virtual void buildSubpath(node v,
			node vPred,
			node vSucc,
			List<adjEntry> &L,
			ExpandedGraph &Exp,
			node s,
			node t);

		Graph::EdgeType	m_typeOfCurrentEdge;
	};


}

#endif
