/*
 * $Revision: 3368 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-04 20:07:31 +0200 (Do, 04. Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration of class CrossingStructure.
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


#ifndef OGDF_CROSSING_STRUCTURE_H
#define OGDF_CROSSING_STRUCTURE_H


#include <ogdf/planarity/PlanRepLight.h>


namespace ogdf {

	class CrossingStructure
	{
	public:
		CrossingStructure() : m_numCrossings(0) { }

		void init(PlanRepLight &PG, int weightedCrossingNumber);
		void restore(PlanRep &PG, int cc);

		int numberOfCrossings() const { return m_numCrossings; }
		int weightedCrossingNumber() const { return m_weightedCrossingNumber; }
		const SListPure<int> &crossings(edge e) const { return m_crossings[e]; }

	private:
		int m_numCrossings;
		int m_weightedCrossingNumber;
		EdgeArray<SListPure<int> > m_crossings;
	};

}

#endif
