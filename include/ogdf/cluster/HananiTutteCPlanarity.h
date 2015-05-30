/** \file
* \brief Defines class HananiTutteCPlanarity, which represents a
*        c-planarity test based on the Hanani-Tutte theorem.
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


#ifndef OGDF_HANANI_TUTTE_CPLANARITY_H
#define OGDF_HANANI_TUTTE_CPLANARITY_H


#include <ogdf/cluster/ClusterGraph.h>


namespace ogdf {


	//! C-planarity testing via Hanani-Tutte approach.
	/**
	 * @ingroup ga-cplanarity
	 */
	class HananiTutteCPlanarity {

		class CLinearSystem;
		class CGraph;

	public:

		enum class Solver { HananiTutte, HananiTutteVerify, ILP };
		enum class Status { invalid, emptyAfterPreproc, cConnectedAfterPreproc, nonPlanarAfterPreproc, applyHananiTutte, applyILP, timeoutILP, errorILP };
		enum class Verification { cPlanar, cPlanarVerified, nonCPlanarVerified, verificationFailed, timeout };

		enum class Type    : uint16_t { vertex, edge };
		enum class SubType : uint16_t { vertex, cluster, edge, innerCluster, outerCluster, vertexCluster, clusterCluster, crossCluster };

		HananiTutteCPlanarity() {
			m_status = Status::invalid;
		}

		Verification isCPlanar(const ClusterGraph &C, bool doPreproc = true, bool forceSolver = false, Solver solver = Solver::HananiTutte);

		Status status() const { return m_status; }

		void preprocessing(ClusterGraph &C, Graph &G);

		int numNodesPreproc() const { return m_numNodesPreproc; }
		int numEdgesPreproc() const { return m_numEdgesPreproc; }
		int numClustersPreproc() const { return m_numClustersPreproc; }

		int numMatrixRows() const { return m_nRows; }
		int numMatrixCols() const { return m_nCols; }

		int64_t timePrepare() const { return m_tPrepare; }
		int64_t timeCreateSparse() const { return m_tCreateSparse; }
		int64_t timesolve() const { return m_tSolve; }

	private:
		int m_nRows;
		int m_nCols;
		int64_t m_tPrepare;
		int64_t m_tCreateSparse;
		int64_t m_tSolve;

		Status m_status;
		int m_numNodesPreproc;
		int m_numEdgesPreproc;
		int m_numClustersPreproc;
	};


}

#endif
