#pragma once

#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/Logger.h>
#include <ogdf/cluster/ClusterGraphAttributes.h>
#include <ogdf/cluster/sync_plan/PMatching.h>
#include <ogdf/cluster/sync_plan/PQPlanarityComponents.h>
#include <ogdf/cluster/sync_plan/PQPlanarityConsistency.h>
#include <ogdf/cluster/sync_plan/QPartitioning.h>

#include <chrono>
#include <cstdint>
#include <ostream>

#include <NodePCRotation.h>

#ifdef LIKWID_PERFMON

#	include <likwid.h>

#	define PQ_PROFILE_START(regionTag) likwid_markerStartRegion(regionTag);
#	define PQ_PROFILE_STOP(regionTag) likwid_markerStopRegion(regionTag);

#else

#	define PQ_PROFILE_START(regionTag)
#	define PQ_PROFILE_STOP(regionTag)

#endif

#define PQ_OPSTATS


using tpc = std::chrono::high_resolution_clock;
using tp = const std::chrono::time_point<std::chrono::high_resolution_clock>;

inline int64_t dur_ms(const tp::duration& d) {
	return std::chrono::duration_cast<std::chrono::milliseconds>(d).count();
}

inline int64_t dur_ns(const tp::duration& d) {
	return std::chrono::duration_cast<std::chrono::nanoseconds>(d).count();
}

using pc_tree::NodePCRotation;

int sumPNodeDegrees(const pc_tree::PCTree& pct);

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

ostream& operator<<(ostream& os, Operation op);

class PQPlanarity {
	friend class PQPlanarityConsistency;

	friend class PQPlanarityDrawer;

	friend class PQPlanOptions;

	friend std::ostream& operator<<(std::ostream& os, const PQPlanarity& pq);

	friend class UndoContract;

	friend class UndoConvertSmall;

	friend class UndoEncapsulate;

	friend class UndoInitCluster;

	friend class UndoPropagate;

	friend class UndoSimplify;

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

		virtual void undo(PQPlanarity& pq) = 0;

		virtual std::ostream& print(std::ostream& os) const = 0;

		virtual std::string name() const;
	};

	class VerifyPipeBijections : public UndoOperation {
		List<std::tuple<int, int, FrozenPipeBij>> pipes;

	public:
		explicit VerifyPipeBijections(PQPlanarity& pq);

		void undo(PQPlanarity& pq) override;

		ostream& print(ostream& os) const override;
	};

	class ResetIndices : public UndoOperation {
		int max_node, max_edge, count_node, count_edge;

	public:
		explicit ResetIndices(PQPlanarity& pq);

		void undo(PQPlanarity& pq) override;

		ostream& print(ostream& os) const override;
	};

	enum Result { SUCCESS, NOT_APPLICABLE, INVALID_INSTANCE };

	/// Members ///////////////////////////////////////////////////////////////////////////////////////////////////////

public:
	Graph* const G;
	PMatching matchings;
	QPartitioning partitions;
	int simplifyToroidalCycleLength = 0;

#ifdef PQ_OPSTATS
	std::ofstream stats_out;
	uint64_t stats_pc_time = 0;
	bool stats_first_in_array = true;
#endif

private:
	List<UndoOperation*> undo_stack;
	PQPlanarityComponents components;
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
	PQPlanarityConsistency consistency;
#endif

	/// Constructors //////////////////////////////////////////////////////////////////////////////////////////////////

public:
	/**
     * Usage:
     *
     * Graph G;
     * // init your graph here
     *
     * PQPlanarity PQ(&G);
     * // init pipes and partitions here:
     * PQ.matchings.matchNodes(u, v);
     * int p = PQ.partitions.makeQVertex(u1);
     * PQ.partitions.makeQVertex(u2, p);
     * PQ.partitions.makeQVertex(u3, p);
     *
     * // if you made any changes to edges in G after creating the PQPlanarity instance, call
     * PQ.initComponents(); // to update the BC-tree if connectivity changed
     * PQ.matchings.rebuildHeap(); // to sort the pipes if their degree changed
     *
     * if (PQ.makeReduced() && PQ.solveReduced()) {
     *     PQ.embed();
     *     OGDF_ASSERT(G.representsCombEmbedding());
     * }
     */
	explicit PQPlanarity(Graph* g, GraphAttributes* ga = nullptr);

	/**
     * Usage:
     *
     * Graph G;
     * ClusterGraph CG(G);
     * // init your graph and clusters here
     *
     * PQPlanarity PQ(&G, &CG);
     * // you shouldn't change G or CG after creating a PQPlanarity instance from them
     *
     * if (PQ.makeReduced() && PQ.solveReduced()) {
     *     PQ.embed();
     *     OGDF_ASSERT(G.representsCombEmbedding());
     *     OGDF_ASSERT(isClusterPlanarEmbedding(CG));
     * }
     */
	explicit PQPlanarity(Graph* g, ClusterGraph* cg, ClusterGraphAttributes* ga = nullptr);

	explicit PQPlanarity(const Graph* sefe, Graph* work, EdgeArray<uint8_t>& edge_types);

	virtual ~PQPlanarity() {
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

	const PQPlanarityComponents& getComponents() const { return components; }

#ifdef PQ_OPSTATS

	void printOPStatsStart(const Pipe* p, Operation op, const NodePCRotation* pct = nullptr);

	void printOPStatsEnd(bool success, int64_t time_ns);

#endif

	PipeType getPipeType(const Pipe* p);

	bool canContract(const Pipe* p);

	bool isAllowContractBBPipe() const { return allow_contract_bb_pipe; }

	void setAllowContractBBPipe(bool allowContractBbPipe) {
		allow_contract_bb_pipe = allowContractBbPipe;
	}

	bool isIntersectTrees() const { return intersect_trees; }

	void setIntersectTrees(bool intersectTrees) { intersect_trees = intersectTrees; }

	bool isBatchSpqr() const { return batch_spqr; }

	void setBatchSpqr(bool batchSpqr) { batch_spqr = batchSpqr; }
};

std::ostream& operator<<(std::ostream& os, const PQPlanarity& pq);

std::ostream& operator<<(std::ostream& os, const PQPlanarity::UndoOperation& undo_op);
