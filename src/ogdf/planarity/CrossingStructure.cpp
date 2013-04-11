/*
 * $Revision: 3368 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-04 20:07:31 +0200 (Do, 04. Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of class CrossingStructure.
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


#include <ogdf/internal/planarity/CrossingStructure.h>


namespace ogdf {


	void CrossingStructure::init(PlanRepLight &PG, int weightedCrossingNumber)
	{
		m_weightedCrossingNumber = weightedCrossingNumber;
		m_crossings.init(PG.original());

		m_numCrossings = 0;
		NodeArray<int> index(PG,-1);
		node v;
		forall_nodes(v,PG)
			if(PG.isDummy(v))
				index[v] = m_numCrossings++;

		edge ePG;
		forall_edges(ePG,PG)
		{
			if(PG.original(ePG->source()) != 0) {
				edge e = PG.original(ePG);
				ListConstIterator<edge> it = PG.chain(e).begin();
				for(++it; it.valid(); ++it) {
					m_crossings[e].pushBack(index[(*it)->source()]);
				}
			}
		}
	}

	void CrossingStructure::restore(PlanRep &PG, int cc)
	{
		Array<node> id2Node(0,m_numCrossings-1,0);

		SListPure<edge> edges;
		PG.allEdges(edges);

		for(SListConstIterator<edge> itE = edges.begin(); itE.valid(); ++itE)
		{
			edge ePG = *itE;
			edge e = PG.original(ePG);

			SListConstIterator<int> it;
			for(it = m_crossings[e].begin(); it.valid(); ++it)
			{
				node x = id2Node[*it];
				edge ePGOld = ePG;
				ePG = PG.split(ePG);
				node y = ePG->source();

				if(x == 0) {
					id2Node[*it] = y;
				} else {
					PG.moveTarget(ePGOld, x);
					PG.moveSource(ePG, x);
					PG.delNode(y);
				}
			}
		}

	}

}
