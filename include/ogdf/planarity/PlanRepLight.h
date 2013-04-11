/*
 * $Revision: 3368 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-04 20:07:31 +0200 (Do, 04. Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration of class PlanRepLight.
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

#ifndef OGDF_PLAN_REP_LIGHT_H
#define OGDF_PLAN_REP_LIGHT_H

#include <ogdf/planarity/PlanRep.h>


namespace ogdf {


	//! Light-weight version of a planarized representation, associated with a PlanRep.
	class PlanRepLight : public GraphCopy
	{
		const CCsInfo &m_ccInfo;
		const PlanRep &m_pr;

		int m_currentCC;
		EdgeArray<edge> m_eAuxCopy;

	public:
		//! Creates a light-weight planarized representation.
		PlanRepLight(const PlanRep &pr);

		//! Returns the number of connected components in the original graph.
		int numberOfCCs() const { return m_ccInfo.numberOfCCs(); }

		//! Returns the index of the current connected component.
		int currentCC() const { return m_currentCC; }

		//! Returns the connected component info structure.
		const CCsInfo &ccInfo() const { return m_ccInfo; }

		//! Returns the original edge with index \a i.
		edge e(int i) const { return m_ccInfo.e(i); }

		//! Returns the original node with index \a i.
		node v(int i) const { return m_ccInfo.v(i); }

		//! Returns the index of the first edge in this connected component.
		int startEdge() const { return m_ccInfo.startEdge(m_currentCC); }

		//! Returns the index of (one past) the last edge in this connected component.
		int stopEdge() const { return m_ccInfo.stopEdge(m_currentCC); }

		EdgeType typeOf(edge e) const {
			edge eOrig = m_eOrig[e];
			return (eOrig != 0) ? typeOrig(eOrig) : Graph::association;
		}

		EdgeType typeOrig(edge eOrig) const { return m_pr.typeOrig(eOrig); }

		//! Initializes the planarized representation for connected component \a cc.
		void initCC(int cc);
	};

} // end namespace ogdf


#endif
