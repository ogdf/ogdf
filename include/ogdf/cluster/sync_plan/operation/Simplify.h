#pragma once

#include <ogdf/basic/Graph.h>
#include <ogdf/cluster/sync_plan/PQPlanarity.h>

#include <NodePCRotation.h>

using RegisteredEdgeSet = RegisteredElementSet<edge, Graph>;

struct SimplifyMapping {
	adjEntry u2_adj, u_adj, v2_adj;
	List<adjEntry> v_adj;

	explicit SimplifyMapping(adjEntry u2Adj = nullptr, adjEntry uAdj = nullptr,
			adjEntry vAdj = nullptr, adjEntry v2Adj = nullptr);

	friend ostream& operator<<(ostream& os, const SimplifyMapping& mapping);
};

struct FrozenSimplifyMapping {
	int u2_adj, u_adj, v2_adj;
	List<int> v_adj;

	explicit FrozenSimplifyMapping(int u2Adj = -1, int uAdj = -1, int vAdj = -1, int v2Adj = -1);

	friend ostream& operator<<(ostream& os, const FrozenSimplifyMapping& mapping);
};

class UndoSimplify : public PQPlanarity::UndoOperation {
	int u2_idx, u_idx, v_idx, v2_idx;
	SList<FrozenSimplifyMapping> bij;

public:
	UndoSimplify(const List<SimplifyMapping>& in_bij, node u2, node u, node v, node v2 = nullptr);

	void undo(PQPlanarity& pq) override;

	ostream& print(ostream& os) const override;
};

/**
 * Find all cycles the permutation \p mapping defines on the order of the adjEntries of \p u
 * and put all adjEntries of each cycle into their own list in \p cycles.
 */
int findCycles(Graph& G, node u, const std::function<adjEntry(adjEntry)>& mapping,
		List<List<adjEntry>>& cycles, std::ostream& log);

/**
 * Continue a DFS from the last node of \p dfs_stack until an unvisited edge incident to \p v is found.
 * All edges visited on the way (also those that didn't lead to an edge incident to v
 * and the edge incident to v found in the end) are added to the set \p visited.
 * @param u_adj (for debugging purposes) an adjEntry at the bottom of the stack, whose node we should never visit
 * @return an adjEntry of \p v reachable from \p dfs_stack.back() via previously unvisited edges
 */
adjEntry continueNodeDFS(Graph& G, adjEntry u_adj, node v, RegisteredEdgeSet& visited,
		List<adjEntry>& dfs_stack);

/**
 * Use a DFS to find an edge incident to \p v reachable from \p u_adj without crossing the node u.
 * All edges visited on the way (also those that didn't lead to an edge incident to v
 * and the edge incident to v found in the end) are added to the set \p visited.
 * The path taken is recorded in \p dfs_stack, where dfs_stack.front() == u_adj and dfs_stack.back()
 * is the adjEntry of v which is also returned.
 * If you want to find a further alternative adjEntries that match this criterion, pop the last entry
 * off the stack and call continueNodeDFS (see also exhaustiveNodeDFS).
 */
adjEntry startNodeDFS(Graph& G, adjEntry u_adj, node v, RegisteredEdgeSet& visited,
		List<adjEntry>& dfs_stack);

/**
 * Use a DFS to find all edges incident to \p v reachable from \p u_adj without crossing the node u and add them to \p out.
 * All edges visited on the way (also those that didn't lead to an edge incident to v
 * and the edge incident to v found in the end) are added to the set \p visited.
 * Calls startNodeDFS once and then repeatedly calls continueNodeDFS until no new adjEntry is found.
 */
int exhaustiveNodeDFS(Graph& G, adjEntry u_adj, node v, RegisteredEdgeSet& visited,
		List<adjEntry>& out);

#ifdef OGDF_DEBUG

/**
 * Check that the method exhaustiveNodeDFS would return the same adjEntries as contained in \p found.
 */
bool compareWithExhaustiveNodeDFS(Graph& G, adjEntry u_adj, node v,
		const RegisteredEdgeSet& visited, const List<adjEntry>& found);

/**
 * For two poles u,v of an (SPQR-Tree) parallel, check that their embedding PCTrees correctly represent all parallel bundles between u and v.
 * Note: \p u_pc must be trivial, while \p v_pc may be non-trivial.
 */
void validatePartnerPCTree(const NodePCRotation* u_pc, const NodePCRotation* v_pc);

/**
 * Check that /p bij_list contains all adjEntries of \p v which belong to a parallel bundle leading to \p u.
 */
bool validateCollectedAdjs(node v, node u, List<SimplifyMapping>& bij_list,
		RegisteredEdgeSet& visited, PQPlanarityComponents& components);

#endif
