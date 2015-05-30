/** \file
 * \brief Declaration of class MFUPS, which implements
 *        the maximal feasible upward planar subgraph computation based on
 * 		  satisfiability (Chimani, Zeranski, 2012+)
 *
 * \author Robert Zeranski, Markus Chimani
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

#ifndef OGDF_MAXIMALFUPS_H
#define OGDF_MAXIMALFUPS_H

#include <ogdf/module/FUPSModule.h>

namespace ogdf {

//---------------------------------------------------------
// MaximalFUPS
//---------------------------------------------------------
class MaximalFUPS : public FUPSModule {
	public:
		//constructor
	MaximalFUPS() : m_timelimit(0) {};
	//		(const Graph &Orig, int timelimit);

	private:
		int m_timelimit;

	protected:
		Module::ReturnType doCall(UpwardPlanRep &UPR, List<edge> &delEdges);
	public:
//		int computeMFUPS(GraphCopy &GC);
		int getTimelimit()                 { return m_timelimit;      }
		void setTimelimit(int timelimit)   { m_timelimit = timelimit; }
};

} // end namespace ogdf

#endif
