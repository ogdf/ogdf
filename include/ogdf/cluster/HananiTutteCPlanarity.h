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
 * See README.md in the OGDF root directory for details.
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
 * License along with this program; if not, see
 * http://www.gnu.org/copyleft/gpl.html
 */

#pragma once

#include <ogdf/basic/basic.h>
#include <ogdf/cluster/ClusterPlanarityModule.h>

#include <cstdint>
#include <stdexcept>

namespace ogdf {
class ClusterGraph;
class Graph;

//! C-planarity testing via Hanani-Tutte approach.
/**
 * @ingroup ga-cplanarity
 */
class OGDF_EXPORT HananiTutteCPlanarity : public ClusterPlanarityModule {
	class CGraph;
	class CLinearSystem;

public:
	struct Stats {
		int nRows = 0;
		int nColumns = 0;
		int nConditions = 0;
		int nMoves = 0;

		int64_t tPrepare = 0;
		int64_t tCreateSparse = 0;
		int64_t tSolve = 0;
		int64_t tCheck = 0;
	};

	struct HananiTutteSolver {
		virtual ~HananiTutteSolver() = default;
		virtual bool test(Stats& stats) = 0;
		virtual bool verify(Stats& stats) = 0;
	};

	enum class Solver { HananiTutte, HananiTutteVerify, ILP };
	enum class Status {
		invalid,
		emptyAfterPreproc,
		cConnectedAfterPreproc,
		nonPlanarAfterPreproc,
		applyHananiTutte,
		applyILP,
		timeoutILP,
		errorILP
	};
	enum class Verification {
		cPlanar,
		cPlanarVerified,
		nonCPlanarVerified,
		verificationFailed,
		timeout
	};

	enum class Type : uint16_t { tVertex, tEdge };
	enum class SubType : uint16_t {
		stVertex,
		stCluster,
		stEdge,
		stInnerCluster,
		stOuterCluster,
		stVertexCluster,
		stClusterCluster,
		stCrossCluster
	};

	bool isClusterPlanar(const ClusterGraph& CG) override {
		Verification res = isCPlanar(CG, true, false, Solver::HananiTutteVerify);
		if (res == Verification::cPlanar || res == Verification::cPlanarVerified) {
			return true;
		} else if (res == Verification::nonCPlanarVerified) {
			return false;
		} else {
			throw std::runtime_error("Could not solve instance!");
		}
	}

	bool isClusterPlanarDestructive(ClusterGraph& CG, Graph& G) override {
		return isClusterPlanar(CG);
	}

	Verification isCPlanar(const ClusterGraph& C, bool doPreproc = true, bool forceSolver = false,
			Solver solver = Solver::HananiTutte);

	Status status() const { return m_status; }

	//! @sa ogdf::sync_plan::preprocessClusterGraph()
	static void preprocessing(ClusterGraph& C, Graph& G);

	int numNodesPreproc() const { return m_numNodesPreproc; }

	int numEdgesPreproc() const { return m_numEdgesPreproc; }

	int numClustersPreproc() const { return m_numClustersPreproc; }

	int numMatrixRows() const { return m_stats.nRows; }

	int numMatrixCols() const { return m_stats.nColumns; }

	int64_t timePrepare() const { return m_stats.tPrepare; }

	int64_t timeCreateSparse() const { return m_stats.tCreateSparse; }

	int64_t timesolve() const { return m_stats.tSolve; }

	const Stats& stats() const { return m_stats; }

	static HananiTutteSolver* getSolver(const ClusterGraph& C);

private:
	Stats m_stats;
	Status m_status = Status::invalid;
	int m_numNodesPreproc = 0;
	int m_numEdgesPreproc = 0;
	int m_numClustersPreproc = 0;
};


}
