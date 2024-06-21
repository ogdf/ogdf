/** \file
 * \brief TODO Document
 *
 * \author Simon D. Fink <ogdf@niko.fink.bayern>
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


#include <ogdf/basic/Graph.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/Logger.h>
#include <ogdf/basic/basic.h>
#include <ogdf/cluster/sync_plan/PMatching.h>
#include <ogdf/cluster/sync_plan/QPartitioning.h>
#include <ogdf/cluster/sync_plan/SyncPlanComponents.h>
#include <ogdf/cluster/sync_plan/SyncPlanConsistency.h>
#include <ogdf/cluster/sync_plan/utils/Bijection.h>

#include <cstdint>
#include <functional>
#include <ostream>
#include <string>
#include <tuple>

namespace ogdf::pc_tree {
class NodePCRotation;
class PCTree;
} // namespace ogdf::pc_tree

namespace ogdf {
class ClusterGraph;
class ClusterGraphAttributes;
class GraphAttributes;
} // namespace ogdf

using ogdf::pc_tree::NodePCRotation;

// // Profiling with LIKWID
// #ifdef LIKWID_PERFMON
// 	include <likwid.h>
// 	define SYNCPLAN_PROFILE_START(regionTag) likwid_markerStartRegion(regionTag);
// 	define SYNCPLAN_PROFILE_STOP(regionTag) likwid_markerStopRegion(regionTag);
// #else
// 	define SYNCPLAN_PROFILE_START(regionTag)
// 	define SYNCPLAN_PROFILE_STOP(regionTag)
// #endif

// // Collection of JSON Operation Statistics
// define SYNCPLAN_OPSTATS

namespace ogdf::sync_plan {

namespace internal {
using tpc = std::chrono::high_resolution_clock;
using tp = const std::chrono::time_point<std::chrono::high_resolution_clock>;

inline int64_t dur_ms(const tp::duration& d) {
	return std::chrono::duration_cast<std::chrono::milliseconds>(d).count();
}

inline int64_t dur_ns(const tp::duration& d) {
	return std::chrono::duration_cast<std::chrono::nanoseconds>(d).count();
}

int sumPNodeDegrees(const ogdf::pc_tree::PCTree& pct);

class UndoSimplify;
}

enum class Operation {
	ENCAPSULATE_CONTRACT,
	CONTRACT_BICON,
	PROPAGATE_CUT,
	PROPAGATE_BICON,
	SIMPLIFY_TERMINAL,
	SIMPLIFY_TRANSITIVE,
	SIMPLIFY_TOROIDAL,
	BATCH_SPQR
};

std::ostream& operator<<(std::ostream& os, Operation op);

class SyncPlan {
	friend class SyncPlanConsistency;

	friend class SyncPlanDrawer;

	friend class SyncPlanOptions;

	friend std::ostream& operator<<(std::ostream& os, const SyncPlan& pq);

	friend class UndoContract;

	friend class UndoConvertSmall;

	friend class UndoEncapsulate;

	friend class UndoInitCluster;

	friend class UndoPropagate;

	friend class internal::UndoSimplify;

	friend class UndoSimplifyToroidal;

	friend class UndoMakeWheel;

	friend class EmbeddingTrees;

	/// Inner Classes /////////////////////////////////////////////////////////////////////////////////////////////////

public:
	class UndoOperation {
	public:
#ifdef OGDF_DEBUG
		int consistency_nr = -1;
#endif

		virtual ~UndoOperation() = default;

		virtual void undo(SyncPlan& pq) = 0;

		virtual std::ostream& print(std::ostream& os) const = 0;

		virtual std::string name() const;

		friend std::ostream& operator<<(std::ostream& os, const SyncPlan::UndoOperation& undo_op);
	};

	class VerifyPipeBijections : public UndoOperation {
		List<std::tuple<int, int, FrozenPipeBij>> pipes;

	public:
		explicit VerifyPipeBijections(SyncPlan& pq);

		void undo(SyncPlan& pq) override;

		std::ostream& print(std::ostream& os) const override;
	};

	class ResetIndices : public UndoOperation {
		int max_node, max_edge, count_node, count_edge;

	public:
		explicit ResetIndices(SyncPlan& pq);

		void undo(SyncPlan& pq) override;

		std::ostream& print(std::ostream& os) const override;
	};

	enum class Result { SUCCESS, NOT_APPLICABLE, INVALID_INSTANCE };

	/// Members ///////////////////////////////////////////////////////////////////////////////////////////////////////

public:
	Graph* const G;
	PMatching matchings;
	QPartitioning partitions;
	int simplifyToroidalCycleLength = 0;

#ifdef SYNCPLAN_OPSTATS
	std::ofstream stats_out;
	uint64_t stats_pc_time = 0;
	bool stats_first_in_array = true;
#endif

private:
	List<UndoOperation*> undo_stack;
	SyncPlanComponents components;
	GraphAttributes* GA;
	NodeArray<node> node_reg;
	EdgeArray<edge> edge_reg;
	NodeArray<bool> is_wheel;
	Logger log;
	bool indices_saved = false;
	bool allow_contract_bb_pipe = true;
	bool intersect_trees = true;
	bool batch_spqr = true;

#ifdef OGDF_DEBUG
	SyncPlanConsistency consistency;
#endif

	/// Constructors //////////////////////////////////////////////////////////////////////////////////////////////////

public:
	/**
	 * Usage:
	 *
	 * Graph G;
	 * // init your graph here
	 *
	 * SyncPlan PQ(&G);
	 * // init pipes and partitions here:
	 * PQ.matchings.matchNodes(u, v);
	 * int p = PQ.partitions.makeQVertex(u1);
	 * PQ.partitions.makeQVertex(u2, p);
	 * PQ.partitions.makeQVertex(u3, p);
	 *
	 * // if you made any changes to edges in G after creating the SyncPlan instance, call
	 * PQ.initComponents(); // to update the BC-tree if connectivity changed
	 * PQ.matchings.rebuildHeap(); // to sort the pipes if their degree changed
	 *
	 * if (PQ.makeReduced() && PQ.solveReduced()) {
	 *     PQ.embed();
	 *     OGDF_ASSERT(G.representsCombEmbedding());
	 * }
	 */
	explicit SyncPlan(Graph* g, GraphAttributes* ga = nullptr);

	/**
	 * Usage:
	 *
	 * Graph G;
	 * ClusterGraph CG(G);
	 * // init your graph and clusters here
	 *
	 * SyncPlan PQ(&G, &CG);
	 * // you shouldn't change G or CG after creating a SyncPlan instance from them
	 *
	 * if (PQ.makeReduced() && PQ.solveReduced()) {
	 *     PQ.embed();
	 *     OGDF_ASSERT(CG.representsCombEmbedding());
	 * }
	 */
	explicit SyncPlan(Graph* g, ClusterGraph* cg, ClusterGraphAttributes* ga = nullptr);

	explicit SyncPlan(const Graph* sefe, Graph* work, EdgeArray<uint8_t>& edge_types);

	virtual ~SyncPlan() {
		while (!undo_stack.empty()) {
			delete undo_stack.popBackRet();
		}
	}

	/// Graph Utils ////////////////////////////////////////////////////////////////////////////////////////////////////

private:
	void formatNode(node n) const;

	std::function<std::ostream&(std::ostream&)> fmtPQNode(node n, bool include_comp = true) const;

	node nodeFromIndex(int idx) const;

	edge edgeFromIndex(int idx) const;

	void thawPipeBijection(node u, node v, const FrozenPipeBij& in, PipeBij& out) const;

	bool verifyPipeBijection(node u, node v, const FrozenPipeBij& bij) const;

	void pushUndoOperationAndCheck(UndoOperation* operation) {
#ifdef OGDF_DEBUG
		operation->consistency_nr = consistency.getCheckCounter() - 1;
#endif
		undo_stack.pushBack(operation);
		OGDF_ASSERT(consistency.consistencyCheck());
	}

	void pushUndoOperation(UndoOperation* operation) { undo_stack.pushBack(operation); }

	/// Operations /////////////////////////////////////////////////////////////////////////////////////////////////////

public:
	void makeWheel(node centre, bool update_cuts = true);

	void contractWheel(node centre);

	Result convertSmall(node small);

	Result encapsulate(node g_cut);

	Result contract(node u);

	Result propagatePQ(node u, NodePCRotation* pct, NodePCRotation* pct_v = nullptr);

	Result simplify(node u, const NodePCRotation* pc);

	Result checkPCTree(node u);

	Result batchSPQR();

public:
	void initComponents();

	bool makeReduced(int check_planarity_every = 0);

	bool solveReduced(bool fail_fast = false);

	void embed();

	/// Statistics ////////////////////////////////////////////////////////////////////////////////////////////////////

	int undoOperations() const { return undo_stack.size(); }

	const SyncPlanComponents& getComponents() const { return components; }

#ifdef SYNCPLAN_OPSTATS

	void printOPStatsStart(const Pipe* p, Operation op, const NodePCRotation* pct = nullptr);

	void printOPStatsEnd(bool success, int64_t time_ns);

#endif

	PipeType getPipeType(const Pipe* p);

	bool canContract(const Pipe* p);

	/// Config ////////////////////////////////////////////////////////////////////////////////////////////////////////

	bool isAllowContractBBPipe() const { return allow_contract_bb_pipe; }

	void setAllowContractBBPipe(bool allowContractBbPipe) {
		allow_contract_bb_pipe = allowContractBbPipe;
	}

	bool isIntersectTrees() const { return intersect_trees; }

	void setIntersectTrees(bool intersectTrees) { intersect_trees = intersectTrees; }

	bool isBatchSpqr() const { return batch_spqr; }

	void setBatchSpqr(bool batchSpqr) { batch_spqr = batchSpqr; }
};

}
