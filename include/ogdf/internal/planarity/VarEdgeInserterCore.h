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

#ifndef OGDF_VAR_EDGE_INS_CORE_H
#define OGDF_VAR_EDGE_INS_CORE_H


#include <ogdf/basic/Timeouter.h>
#include <ogdf/basic/Module.h>
#include <ogdf/planarity/PlanRepLight.h>
#include <ogdf/planarity/RemoveReinsertType.h>


namespace ogdf {

	class StaticSPQRTree;


	class OGDF_EXPORT VarEdgeInserterCore : public Timeouter
	{
	public:
		VarEdgeInserterCore(
			PlanRepLight &pr,
			const EdgeArray<int>      *pCostOrig,
			const EdgeArray<bool>     *pForbiddenOrig,
			const EdgeArray<__uint32> *pEdgeSubgraphs)
			: m_pr(pr), m_pCost(pCostOrig), m_pForbidden(pForbiddenOrig), m_pSubgraph(pEdgeSubgraphs) { }

		virtual ~VarEdgeInserterCore() { }

		Module::ReturnType call(
			const Array<edge> &origEdges,
			RemoveReinsertType rrPost,
			double percentMostCrossed);

		Module::ReturnType callPostprocessing(
			const Array<edge> &origEdges,
			RemoveReinsertType rrPost,
			double percentMostCrossed);

		int runsPostprocessing() const { return m_runsPostprocessing; }

	protected:
		class BiconnectedComponent;
		class ExpandedGraph;

		void insert(node s, node t, SList<adjEntry> &eip);
		int costCrossed(edge eOrig) const;

		bool dfsVertex(node v, int parent);
		bool dfsComp(int i, node parent, node &repT);

		void blockInsert(
			const BiconnectedComponent &G,
			node s,
			node t,
			List<adjEntry> &L);
		bool pathSearch(node v, edge parent, List<edge> &path);
		virtual void buildSubpath(node v,
			edge eIn,
			edge eOut,
			List<adjEntry> &L,
			ExpandedGraph &Exp,
			node s,
			node t);

		virtual void storeTypeOfCurrentEdge(edge eOrig) { }
		virtual BiconnectedComponent *createBlock();
		virtual ExpandedGraph *createExpandedGraph(const BiconnectedComponent &BC, const StaticSPQRTree &T);

		static const int c_bigM = 10000;
		PlanRepLight	&m_pr;

		const EdgeArray<int>		*m_pCost;
		const EdgeArray<bool>		*m_pForbidden;
		const EdgeArray<__uint32>	*m_pSubgraph;

		node   m_s, m_t;
		edge   m_st;
		SList<adjEntry> *m_pEip;

		// representation of BC tree
		NodeArray<SList<int> > m_compV;
		Array<SList<node> >    m_nodeB;
		Array<SList<edge> >    m_edgeB;
		NodeArray<node>        m_GtoBC;

		node m_v1, m_v2;

		int m_runsPostprocessing; //!< Runs of remove-reinsert method.
	};


	class VarEdgeInserterUMLCore : public VarEdgeInserterCore
	{
	public:
		VarEdgeInserterUMLCore(
			PlanRepLight &pr,
			const EdgeArray<int>      *pCostOrig,
			const EdgeArray<__uint32> *pEdgeSubgraph) : VarEdgeInserterCore(pr, pCostOrig, 0, pEdgeSubgraph) { }

	protected:
		class BiconnectedComponentUML;
		class ExpandedGraphUML;

		void storeTypeOfCurrentEdge(edge eOrig) { m_typeOfCurrentEdge = m_pr.typeOrig(eOrig); }
		BiconnectedComponent *createBlock();
		ExpandedGraph *createExpandedGraph(const BiconnectedComponent &BC, const StaticSPQRTree &T);
		virtual void buildSubpath(node v,
			edge eIn,
			edge eOut,
			List<adjEntry> &L,
			ExpandedGraph &Exp,
			node s,
			node t);

		Graph::EdgeType	m_typeOfCurrentEdge;
	};

}

#endif
