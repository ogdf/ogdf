/*
 * $Revision: 3417 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-18 10:27:30 +0200 (Do, 18. Apr 2013) $
 ***************************************************************/

/** \file
 * \brief implementation of FixedEmbeddingInserterUML class
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


#include <ogdf/uml/FixedEmbeddingInserterUML.h>
#include <ogdf/internal/planarity/FixEdgeInserterCore.h>


namespace ogdf {

	//---------------------------------------------------------
	// constructor
	// sets default values for options
	//
	FixedEmbeddingInserterUML::FixedEmbeddingInserterUML()
	{
		m_rrOption = rrNone;
		m_percentMostCrossed = 25;
		m_keepEmbedding = false;
	}


	// copy constructor
	FixedEmbeddingInserterUML::FixedEmbeddingInserterUML(const FixedEmbeddingInserterUML &inserter)
		: UMLEdgeInsertionModule()
	{
		m_rrOption = inserter.m_rrOption;
		m_percentMostCrossed = inserter.m_percentMostCrossed;
		m_keepEmbedding = inserter.m_keepEmbedding;
	}


	// clone method
	UMLEdgeInsertionModule *FixedEmbeddingInserterUML::clone() const
	{
		FixedEmbeddingInserterUML *pInserter = new FixedEmbeddingInserterUML;
		pInserter->m_rrOption = m_rrOption;
		pInserter->m_percentMostCrossed = m_percentMostCrossed;
		pInserter->m_keepEmbedding = m_keepEmbedding;

		return pInserter;
	}


	// assignment operator
	FixedEmbeddingInserterUML &FixedEmbeddingInserterUML::operator=(const FixedEmbeddingInserterUML &inserter)
	{
		m_rrOption = inserter.m_rrOption;
		m_percentMostCrossed = inserter.m_percentMostCrossed;
		m_keepEmbedding = inserter.m_keepEmbedding;
		return *this;
	}


	// actual call method
	Module::ReturnType FixedEmbeddingInserterUML::doCall(
		PlanRepLight              &pr,
		const Array<edge>         &origEdges,
		const EdgeArray<int>      *pCostOrig,
		const EdgeArray<__uint32> *pEdgeSubgraph)
	{
		FixEdgeInserterUMLCore core(pr, pCostOrig, pEdgeSubgraph);
		core.timeLimit(timeLimit());

		return core.call(origEdges, m_keepEmbedding, m_rrOption, m_percentMostCrossed);
	}

}
