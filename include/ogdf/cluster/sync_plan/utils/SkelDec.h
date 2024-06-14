#pragma once

#include <ogdf/basic/Graph_d.h>
#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/NodeArray.h>

using namespace ogdf;

class SkelDec {
	Graph *G;
	Graph skel;
	NodeArray<node> origV;
	EdgeArray<edge> origE;
	EdgeArray<bool> virtE;

	EdgeArray<edge> r_origE;

	void reset() {
		origV.init(skel, nullptr);
		origE.init(skel, nullptr);
		virtE.init(skel, false);
		r_origE.init(*G, nullptr);
	}

	void consistencyCheck() {
		// V and E are exactly the images of origV and origE
		// origV and twinE+origE are total

		// 1 (bicon) Each skeleton is biconnected.
		// 2 (tree) Graph SPQR is simple, loop-free, connected and acyclic, i.e., a tree.
		// 3 (orig-inj) For each skeleton Gµ , the restriction origV|Vµ is injective.
		// 4 (orig-real) For each real edge uv, the endpoints of origE(uv) are origV(u) and origV(v).
		// 5 (orig-virt) Let uv and u′v′ be two virtual edges with uv = twinE(u′v′).
		//      For their respective skeletons Gµ and G′µ (where µ and µ′ are adjacent in SPQR),
		//      it is origV(Vµ) ∩ origV(Vµ′) = origV({u, v}) = origV({u′, v′}).
		// 6 (subgraph) The allocation skeletons of any vertex of GS form a connected subgraph of SPQR.

		// (SPQR 1) each skeleton is either a polygon, a bond, or triconnected ("rigid")
		// (SPQR 2) two skeletons adjacent in SPQR are never both polygons or both bonds
	}
};