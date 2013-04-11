/*
 * $Revision: 3368 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-04 20:07:31 +0200 (Do, 04. Apr 2013) $
 ***************************************************************/

/** \file
 * \brief implementation of FixedEmbeddingInserter class
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


#include <ogdf/planarity/FixedEmbeddingInserter.h>
#include <ogdf/internal/planarity/FixEdgeInserterCore.h>


namespace ogdf {

	//---------------------------------------------------------
	// constructor
	// sets default values for options
	//
	FixedEmbeddingInserter::FixedEmbeddingInserter()
	{
		m_rrOption = rrNone;
		m_percentMostCrossed = 25;
		m_keepEmbedding = false;
	}


	// copy constructor
	FixedEmbeddingInserter::FixedEmbeddingInserter(const FixedEmbeddingInserter &inserter)
		: EdgeInsertionModule(inserter)
	{
		m_rrOption = inserter.m_rrOption;
		m_percentMostCrossed = inserter.m_percentMostCrossed;
		m_keepEmbedding = inserter.m_keepEmbedding;
	}


	// clone method
	EdgeInsertionModule *FixedEmbeddingInserter::clone() const
	{
		return new FixedEmbeddingInserter(*this);
	}


	// assignment operator
	FixedEmbeddingInserter &FixedEmbeddingInserter::operator=(const FixedEmbeddingInserter &inserter)
	{
		m_timeLimit = inserter.m_timeLimit;
		m_rrOption = inserter.m_rrOption;
		m_percentMostCrossed = inserter.m_percentMostCrossed;
		m_keepEmbedding = inserter.m_keepEmbedding;
		return *this;
	}


	// actual call method
	Module::ReturnType FixedEmbeddingInserter::doCall(
		PlanRepLight &pr,
		const Array<edge> &origEdges,
		const EdgeArray<int> *pCostOrig,
		const EdgeArray<bool> *pForbiddenOrig,
		const EdgeArray<__uint32> *pEdgeSubgraphs)
	{
		FixEdgeInserterCore core(pr, pCostOrig, pForbiddenOrig, pEdgeSubgraphs);
		core.timeLimit(timeLimit());

		ReturnType retVal = core.call(origEdges, m_keepEmbedding, m_rrOption, m_percentMostCrossed);
		m_runsPostprocessing = core.runsPostprocessing();
		return retVal;
	}

}
