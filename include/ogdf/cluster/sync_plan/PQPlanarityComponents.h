#pragma once

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/NodeSet.h>
#include <ogdf/decomposition/BCTree.h>

#include "operation/Encapsulate.h"
#include "utils/GraphIterators.h"

using namespace ogdf;

class PQPlanarityComponents {
	friend class PQPlanarity;

	friend class PQPlanarityConsistency;

	friend class PQPlanarityDrawer;

private:
	Graph* G;
	Graph BC;
	int conn_count, conn_next_id;

	NodeArray<BCTree::BNodeType> bc_type;
	NodeArray<int> bc_conn_id;
	NodeArray<int> bc_size;
	NodeArray<node> bc_g;
	NodeArray<node> g_bc;

	mutable int counter = 1;
	mutable NodeArray<std::pair<int, edge>> marker; // used by graphEdgeToBCEdge

public:
	explicit PQPlanarityComponents(Graph* g) : G(g), BC(), marker(BC, std::pair(0, nullptr)) {
		reset();
	}

	int connectedCount() const { return conn_count; }

	int connectedId(node n) const { return bcConnectedId(biconnectedComponent(n)); }

	int biconnectedId(node n) const { return biconnectedComponent(n)->index(); }

	node biconnectedComponent(node n) const;

	BCTree::GNodeType nodeType(node n) const;

	bool isCutVertex(node n) const { return nodeType(n) == BCTree::GNodeType::CutVertex; }

	Graph& graph() const { return *G; }

	const Graph& bcTree() const { return BC; }

	BCTree::BNodeType bcNodeType(node n) const;

	bool isCutComponent(node n) const { return bcNodeType(n) == BCTree::BNodeType::CComp; }

	int bcSize(node n) const;

	int bcConnectedId(node n) const;

	node bcRepr(node bc) const;

	node findOtherRepr(node bc) const;

	std::function<std::ostream&(std::ostream&)> fmtBCNode(node bc) const;

	std::function<std::ostream&(std::ostream&)> fmtBCOf(node n) const {
		return fmtBCNode(biconnectedComponent(n));
	}

	///////////////////////////////////////////////////////////////////////////////////////////////

	std::pair<edge, edge> graphEdgeToBCEdge(node bc_src, node bc_tgt) const;

	node findCommonBiconComp(node bc_cut1, node bc_cut2) const;

	void biconReprNodes(node bicon, SList<node>& nodes) const;

	FilteringBFS nodesInBiconnectedComponent(node bicon) const;

	///////////////////////////////////////////////////////////////////////////////////////////////

private:
	void reset() {
		conn_count = conn_next_id = 0;
		bc_conn_id.init(BC, -1);
		g_bc.init(*G, nullptr);
		bc_g.init(BC, nullptr);
		bc_size.init(BC, 0);
		bc_type.init(BC, BCTree::BNodeType::BComp);
	}

	void makeRepr(node bc, node n);

	void nodeInserted(node g_n, node bc_n);

	void insert(BCTree& tmp_bc);

	void cutReplacedByWheel(node centre, const NodeArray<SList<adjEntry>>& block_neigh);

	void relabelExplodedStar(node center1, node center2, List<node> remnants);

	void preJoin(node keep, node merge);

	void postSplitOffEncapsulatedBlock(node cut, EncapsulatedBlock& block);

	void labelIsolatedNodes();
};

class BiconnectedIsolation {
	PQPlanarityComponents& m_comps;
	node m_bicon;
	NodeSet<true> m_to_restore;
	NodeArray<SListPure<adjEntry>> m_adjEntries; // TODO replace by contiguous vector + indices?
	Graph::HiddenEdgeSet m_hiddenEdges;

public:
	explicit BiconnectedIsolation(PQPlanarityComponents& comps, node bicon);

	~BiconnectedIsolation() { restore(); }

	void restore(bool restore_embedding = true);
};
