/** \file
 * \brief Derive embedding trees from Triconnectivity information.
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

#include <ogdf/basic/DisjointSets.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/pctree/NodePCRotation.h>
#include <ogdf/cluster/sync_plan/basic/OverlappingGraphCopies.h>
#include <ogdf/graphalg/Triconnectivity.h>

#include <vector>

namespace ogdf::pc_tree {
class PCNode;
} // namespace ogdf::pc_tree

namespace ogdf {
class Logger;
template<class E>
class List;
} // namespace ogdf

namespace ogdf::sync_plan {

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

//! Wrapper class around Triconnectivity information.
struct OGDF_EXPORT SimpleSPQRTree {
	using Comp = Triconnectivity::CompStruct;
	static Logger log;
	OverlappingGraphCopy GC;
	OverlappingGraphCopies GC_skels;
	EdgeArray<SList<edge>> par_replacement;
	EdgeArray<OverlappingGraphCopy*> skels;
	EdgeArray<OverlappingGraphCopy*> twins;
	std::vector<OverlappingGraphCopy*> skel_array;
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

//! Derive embedding trees from Triconnectivity information.
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

}
