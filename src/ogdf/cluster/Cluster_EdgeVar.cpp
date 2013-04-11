/*
 * $Revision: 3005 $
 *
 * last checkin:
 *   $Author: chimani $
 *   $Date: 2012-11-12 14:19:48 +0100 (Mo, 12. Nov 2012) $
 ***************************************************************/

/** \file
 * \brief Implementation of the variable class for the Branch&Cut algorithm
 * for the Maximum C-Planar SubGraph problem
 *
 * \author Mathias Jansen
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

#include <ogdf/basic/basic.h>

#ifdef USE_ABACUS

#include <ogdf/internal/cluster/Cluster_EdgeVar.h>
#include <ogdf/internal/cluster/MaxCPlanar_Master.h>

using namespace ogdf;
using namespace abacus;

EdgeVar::EdgeVar(Master *master, double obj, edgeType eType, node source, node target) :
	Variable (master, 0, false, false, obj, eType==CONNECT ? 0.0 : (((MaxCPlanarMaster*)master)->getCheckCPlanar() ? 1.0 : 0.0), 1.0, eType==CONNECT ? VarType::Binary : (((MaxCPlanarMaster*)master)->getCheckCPlanar() ? VarType::Continuous : VarType::Binary)) // TODO-TESTING
{
	m_eType = eType;
	m_source = source;
	m_target = target;
//	m_objCoeff = obj; // not necc.
//TODO no searchedge!
	if (eType == ORIGINAL) m_edge = ((MaxCPlanarMaster*)master)->getGraph()->searchEdge(source,target);
	else m_edge = NULL;
}

EdgeVar::EdgeVar(Master *master, double obj, node source, node target) :
	Variable (master, 0, false, false, obj, 0.0, 1.0, VarType::Binary)
{
	m_eType = CONNECT;
	m_source = source;
	m_target = target;
	m_edge = NULL;
}

EdgeVar::EdgeVar(Master *master, double obj, double lbound, node source, node target) :
	Variable (master, 0, false, false, obj, lbound, 1.0, VarType::Binary)
{
	m_eType = CONNECT;
	m_source = source;
	m_target = target;
	m_edge = NULL;
}

EdgeVar::~EdgeVar() {}

#endif // USE_ABACUS
