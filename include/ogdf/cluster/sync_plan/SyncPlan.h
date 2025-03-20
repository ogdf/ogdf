/** \file
 * \brief The main code for modelling and solving Synchronized Planarity instances.
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
#include <ogdf/basic/GraphSets.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/Logger.h>
#include <ogdf/basic/basic.h>
#include <ogdf/cluster/sync_plan/PMatching.h>
#include <ogdf/cluster/sync_plan/QPartitioning.h>
#include <ogdf/cluster/sync_plan/SyncPlanComponents.h>
#include <ogdf/cluster/sync_plan/SyncPlanConsistency.h>
#include <ogdf/cluster/sync_plan/utils/Bijection.h>

#include <chrono>
#include <cstdint>
#include <functional>
#include <ostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

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

OGDF_EXPORT int sumPNodeDegrees(const ogdf::pc_tree::PCTree& pct);

class UndoSimplify;
}

//! The reduction operations (and their distinct cases) implemented by SyncPlan
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

OGDF_EXPORT std::ostream& operator<<(std::ostream& os, Operation op);

//! A class for modelling and solving Synchronized Planarity instances.
/**
 * This implements the algorithm described in the following paper:
 * \remark Thomas Bläsius, Simon D. Fink, and Ignaz Rutter. 2023. Synchronized Planarity with Applications to Constrained Planarity Problems. ACM Trans. Algorithms 19, 4, Article 34 (October 2023), 23 pages. https://doi.org/10.1145/3607474
 *
 * An evaluation of this implementation can be found in the following paper:
 * \remark Simon D. Fink and Ignaz Rutter. 2024. Constrained Planarity in Practice: Engineering the Synchronized Planarity Algorithm. 2024 Proceedings of the Symposium on Algorithm Engineering and Experiments (ALENEX) https://doi.org/10.1137/1.9781611977929.1
 *
 * For more details, see also (open access):
 * \remark Simon D. Fink. 2024. Constrained Planarity Algorithms in Theory and Practice. Doctoral Thesis, University of Passau. https://doi.org/10.15475/cpatp.2024
 */
class OGDF_EXPORT SyncPlan {
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

	// Inner Classes //////////////////////////////////////////////////////////////////////////////////////////////////

public:
	//! The information needed for undoing the changes a specific operation made to the graph while maintaining its embedding.
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

	//! The result of applying a single operation.
	enum class Result { SUCCESS, NOT_APPLICABLE, INVALID_INSTANCE };

	// Members ////////////////////////////////////////////////////////////////////////////////////////////////////////

public:
	//! The underlying graph
	Graph* const G;
	//! Collection of all matched P-nodes and their pipes
	PMatching matchings;
	//! Collection of all Q-nodes and their partitions
	QPartitioning partitions;

#ifdef SYNCPLAN_OPSTATS
	//! The stream to which json statistics should be written
	std::ofstream stats_out;
	uint64_t stats_pc_time = 0;
	bool stats_first_in_array = true;
#endif

private:
	//! Stack of all operations that should be undone when embedding
	List<UndoOperation*> undo_stack;
	//! Data structure maintaining (bi)connected components information
	SyncPlanComponents components;
	//! Set of all edge objects deleted during the reduction, to be restored when undoing all operations.
	Graph::HiddenEdgeSet deletedEdges;
#ifdef OGDF_DEBUG
	//! Set of all node objects deleted during the reduction. Will remain as isolated nodes within \p G.
	NodeSet deletedNodes;
#endif
	//! If non-null, will be updated with debugging information from applied operations.
	GraphAttributes* GA;
	//! A mapping of node-index to node-object used while undoing operations.
	NodeArray<node> node_reg;
	//! A mapping of edge-index to edge-object used while undoing operations.
	EdgeArray<edge> edge_reg;
	//! Labels nodes whether they are already part of a (Q-node-replacement) wheel
	NodeArray<bool> is_wheel;
	//! The logger to use.
	Logger log;
	//! Whether the constructor already saved the maximum indices present in \p G before the reduction.
	bool indices_saved = false;
	//! Whether to allow contracting block-vertex to block-vertex pipes instead of propagating.
	bool allow_contract_bb_pipe = false;
	//! Whether to intersect trees on block-vertex to block-vertex pipes when propagating.
	bool intersect_trees = true;
	//! Whether to apply embedding-tree based operations in batch using SPQR-tree data
	bool batch_spqr = true;
	//! Keeps track of the longest cycle length encountered in simplify toroidal
	int longestSimplifyToroidalCycle = 0;

#ifdef OGDF_DEBUG
	//! Consistency checking utils
	SyncPlanConsistency consistency;
#endif

	// Constructors ///////////////////////////////////////////////////////////////////////////////////////////////////

public:
	//! Create a new Synchronized Planarity instance with a given underlying graph.
	/**
	 * Usage:
	 * \code
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
	 * \endcode
	 *
	 * @param g the underlying graph.
	 * @param ga optional GraphAttributes in which to store debugging information of applied operations.
	 */
	explicit SyncPlan(Graph* g, GraphAttributes* ga = nullptr);

	//! Create a new Synchronized Planarity instance by applying the reduction from some ClusteredPlanarity instance.
	/**
	 * Usage:
	 * \code
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
	 * \endcode
	 *
	 * @note After calling (the reduction used within) this constructor, the ClusterGraph and underlying Graph will be invalid until SyncPlan::embed() is called.
	 * @param g the underlying graph.
	 * @param cg the corresponding ClusterGraph.
	 * @param augmentation if non-null, will be assigned the edges that need to be inserted to make the graph c-connected c-plane once SyncPlan::embed() was called.
	 * @param ga optional GraphAttributes in which to store debugging information of applied operations.
	 * @sa insertAugmentationEdges()
	 */
	explicit SyncPlan(Graph* g, ClusterGraph* cg,
			std::vector<std::pair<adjEntry, adjEntry>>* augmentation = nullptr,
			ClusterGraphAttributes* ga = nullptr);

	//! Create a new SyncPlan instance using the reduction from a 2-SEFE instance with a connected shared graph.
	/**
	 * @param sefe the 2-SEFE instance that shall be embedded.
	 * @param work an empty graph that is used for the reduction.
	 * @param edge_types the types of edges in \p sefe, 1 or 2 for either of the exclusive graphs, 3 for shared.
	 */
	explicit SyncPlan(Graph* sefe, Graph* work, EdgeArray<uint8_t>& edge_types);

	virtual ~SyncPlan() {
		while (!undo_stack.empty()) {
			delete undo_stack.popBackRet();
		}
	}

	// Graph Utils /////////////////////////////////////////////////////////////////////////////////////////////////////

private:
	void formatNode(node n) const;

	std::function<std::ostream&(std::ostream&)> fmtPQNode(node n, bool include_comp = true) const;

	node nodeFromIndex(int idx) const;

	edge edgeFromIndex(int idx) const;

	void thawPipeBijection(node u, node v, const FrozenPipeBij& in, PipeBij& out) const;

	bool verifyPipeBijection(node u, node v, const FrozenPipeBij& bij) const;

	void pushUndoOperationAndCheck(UndoOperation* operation) {
		// room for improvement: don't generate UndoOps if we are never going to embed
#ifdef OGDF_DEBUG
		operation->consistency_nr = consistency.getCheckCounter() - 1;
#endif
		undo_stack.pushBack(operation);
		OGDF_ASSERT(consistency.consistencyCheck());
	}

	void pushUndoOperation(UndoOperation* operation) { undo_stack.pushBack(operation); }

	// Operations //////////////////////////////////////////////////////////////////////////////////////////////////////

public:
	//! @name Operations used by makeReduced() for removing pipes.
	//! @{
	void makeWheel(node centre, bool update_cuts = true);

	void contractWheel(node centre);

	Result convertSmall(node small);

	Result encapsulate(node g_cut);

	Result contract(node u);

	Result propagatePQ(node u, NodePCRotation* pct, NodePCRotation* pct_v = nullptr);

	Result simplify(node u, const NodePCRotation* pc);

	//! Computes an embedding tree and either applies propagatePQ() or simplify()
	Result checkPCTree(node u);

	Result batchSPQR();
	//! @}

public:
	//! @name Methods for solving instances
	//! @{
	//! Recompute biconnected component information after the graph was changed externally.
	void initComponents();

	//! Apply operations (in the order defined by the current PipeQueue) until no pipes are left.
	bool makeReduced(int check_planarity_every = 0);

	//! Solve a reduced instance, creating a combinatorial embedding of the current graph that respects all Q-constraints.
	bool solveReduced(bool fail_fast = false);

	//! Undo all operations while maintaining the current embedding to obtain an embedding of the initial graph.
	void embed();

	//! @}

	// Statistics /////////////////////////////////////////////////////////////////////////////////////////////////////
	//! @name Statistics
	//! @{

	//! The number of operations that need undoing.
	int undoOperations() const { return undo_stack.size(); }

	//! Return the longest cycle length encountered in simplify toroidal.
	int getLongestSimplifyToroidalCycle() const { return longestSimplifyToroidalCycle; }

	//! The maintained (bi)connected components information.
	const SyncPlanComponents& getComponents() const { return components; }

#ifdef SYNCPLAN_OPSTATS

	void printOPStatsStart(const Pipe* p, Operation op, const NodePCRotation* pct = nullptr);

	void printOPStatsEnd(bool success, int64_t time_ns);

#endif

	PipeType getPipeType(const Pipe* p) const;

	//! Check whether a Pipe can be removed by applying contract().
	//! @sa setAllowContractBBPipe()
	bool canContract(const Pipe* p) const;

	//! @}

	/// Config ////////////////////////////////////////////////////////////////////////////////////////////////////////

	//! @name Configuration
	//! @{

	bool isAllowContractBBPipe() const { return allow_contract_bb_pipe; }

	//! Configure whether block-block pipes can be contracted.
	void setAllowContractBBPipe(bool allowContractBbPipe) {
		allow_contract_bb_pipe = allowContractBbPipe;
	}

	bool isIntersectTrees() const { return intersect_trees; }

	//! Configure whether the embedding trees of block-block pipes should be intersected before propagating.
	void setIntersectTrees(bool intersectTrees) { intersect_trees = intersectTrees; }

	bool isBatchSpqr() const { return batch_spqr; }

	//! Configure whether embedding trees should be computed in batch by deriving them from an SPQR-tree.
	void setBatchSpqr(bool batchSpqr) { batch_spqr = batchSpqr; }

	//! @}
};

}
