#pragma once

#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/pctree/NodePCRotation.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/cluster/sync_plan/basic/OverlappingGraphCopies.h>
#include <ogdf/graphalg/Triconnectivity.h>

using namespace ogdf;
using namespace std;

ostream& operator<<(ostream& os, Triconnectivity::CompType t);

namespace spqr_utils {
inline bool isSNode(const Graph& skel) {
	return skel.numberOfNodes() > 2 && skel.numberOfNodes() == skel.numberOfEdges();
}

inline bool isPNode(const Graph& skel) { return skel.numberOfNodes() == 2; }

inline bool isRNode(const Graph& skel) { return !isSNode(skel) && !isPNode(skel); }

inline adjEntry getAdjInSkel(const OverlappingGraphCopy* skel, adjEntry GC_adj) {
	return skel->copy(GC_adj->theEdge())->getAdj(skel->copy(GC_adj->theNode()));
}

inline adjEntry getAdjInOrig(const OverlappingGraphCopy* skel, adjEntry skel_adj) {
	return skel->original(skel_adj->theEdge())->getAdj(skel->original(skel_adj->theNode()));
}
}

struct SimpleSPQRTree {
	using Comp = Triconnectivity::CompStruct;
	static Logger log;
	OverlappingGraphCopy GC;
	OverlappingGraphCopies GC_skels;
	EdgeArray<SList<edge>> par_replacement;
	EdgeArray<OverlappingGraphCopy*> skels;
	EdgeArray<OverlappingGraphCopy*> twins;
	vector<OverlappingGraphCopy*> skel_array;
	bool planar = true;

	OGDF_NO_COPY(SimpleSPQRTree)

	OGDF_NO_MOVE(SimpleSPQRTree)

	SimpleSPQRTree(OverlappingGraphCopies& OGC_base) : GC(OGC_base), GC_skels(GC) {};

	~SimpleSPQRTree() {
		for (auto ptr : skel_array) {
			if (!ptr) {
				continue;
			}
			ptr->breakLinkForMasterDeconstruction();
			delete ptr;
		}
	}

	void init();

	OverlappingGraphCopy* getNonSSkel(node GC_n) const;

	OverlappingGraphCopy* getTwinSkel(OverlappingGraphCopy* skel, edge skel_e) const;

	OverlappingGraphCopy* getTwinSkel_GC(OverlappingGraphCopy* skel, edge GC_e) const;
};

class NodeSSPQRRotation : public pc_tree::NodePCRotation {
	const SimpleSPQRTree& spqr;

	pc_tree::PCNode* process(adjEntry skel_adj, OverlappingGraphCopy& skel,
			pc_tree::PCNode* parent = nullptr);

	void getIncidentRealEdgesInSubtree(adjEntry skel_adj, OverlappingGraphCopy& skel,
			List<edge>& out);

public:
	OGDF_NO_COPY(NodeSSPQRRotation)

	OGDF_NO_MOVE(NodeSSPQRRotation)

	NodeSSPQRRotation(const SimpleSPQRTree& spqr, node n);

	void mapPartnerEdges();
};
