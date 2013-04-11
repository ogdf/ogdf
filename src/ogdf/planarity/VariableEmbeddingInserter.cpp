/*
 * $Revision: 3368 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-04 20:07:31 +0200 (Do, 04. Apr 2013) $
 ***************************************************************/

/** \file
 * \brief implements class VariableEmbeddingInserter
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


#include <ogdf/planarity/VariableEmbeddingInserter.h>
#include <ogdf/internal/planarity/VarEdgeInserterCore.h>


namespace ogdf {

	//---------------------------------------------------------
	// constructor
	// sets default values for options
	//
	VariableEmbeddingInserter::VariableEmbeddingInserter()
	{
		m_rrOption = rrNone;
		m_percentMostCrossed = 25;
	}


	// copy constructor
	VariableEmbeddingInserter::VariableEmbeddingInserter(const VariableEmbeddingInserter &inserter)
		: EdgeInsertionModule(inserter)
	{
		m_rrOption = inserter.m_rrOption;
		m_percentMostCrossed = inserter.m_percentMostCrossed;
	}


	// clone method
	EdgeInsertionModule *VariableEmbeddingInserter::clone() const
	{
		return new VariableEmbeddingInserter(*this);
	}


	// assignment operator
	VariableEmbeddingInserter &VariableEmbeddingInserter::operator=(const VariableEmbeddingInserter &inserter)
	{
		m_timeLimit = inserter.m_timeLimit;
		m_rrOption = inserter.m_rrOption;
		m_percentMostCrossed = inserter.m_percentMostCrossed;
		return *this;
	}


	// actual call method
	Module::ReturnType VariableEmbeddingInserter::doCall(
		PlanRepLight &pr,
		const Array<edge> &origEdges,
		const EdgeArray<int> *pCostOrig,
		const EdgeArray<bool> *pForbiddenOrig,
		const EdgeArray<__uint32> *pEdgeSubgraph)
	{
		VarEdgeInserterCore core(pr, pCostOrig, pForbiddenOrig, pEdgeSubgraph);
		core.timeLimit(timeLimit());

		ReturnType retVal = core.call(origEdges, m_rrOption, m_percentMostCrossed);
		m_runsPostprocessing = core.runsPostprocessing();
		return retVal;
	}


	// actual call method for postprocessing only
	Module::ReturnType VariableEmbeddingInserter::doCallPostprocessing(
			PlanRepLight              &pr,
			const Array<edge>         &origEdges,
			const EdgeArray<int>      *pCostOrig,
			const EdgeArray<bool>     *pForbiddenOrig,
			const EdgeArray<__uint32> *pEdgeSubgraphs)
	{
		VarEdgeInserterCore core(pr, pCostOrig, pForbiddenOrig, pEdgeSubgraphs);
		core.timeLimit(timeLimit());

		ReturnType retVal = core.callPostprocessing(origEdges, m_rrOption, m_percentMostCrossed);
		m_runsPostprocessing = core.runsPostprocessing();
		return retVal;
	}

}
