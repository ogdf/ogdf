#include <ogdf/basic/GraphSets.h>
#include <ogdf/basic/graph_generators/clustering.h>
#include <ogdf/cluster/CconnectClusterPlanar.h>
#include <ogdf/cluster/ClusterPlanarizationLayout.h>
#include <ogdf/cluster/sync_plan/ClusterPlanarity.h>
#include <ogdf/fileformats/GraphIO.h>

using namespace ogdf;

int main(void) {
	// create a random cluster graph
	Graph G;
	ClusterGraph CG(G);
	randomClusterPlanarGraph(G, CG, 5, 20, 50);

	// set up the sync plan test
	SyncPlanClusterPlanarityModule sp;
	std::vector<std::pair<adjEntry, adjEntry>> augmentation;
	sp.setStoredAugmentation(&augmentation);

	// run the test + embedder
	if (!sp.clusterPlanarEmbedClusterPlanarGraph(CG, G)) {
		throw std::runtime_error("Not c-planar!");
	}
	OGDF_ASSERT(CG.adjAvailable());
	OGDF_ASSERT(CG.representsCombEmbedding());

	// add the computed augmentation edges
	EdgeSet<> added(G);
	for (auto& pair : augmentation) {
		added.insert(G.newEdge(pair.first, Direction::after, pair.second, Direction::before));
	}
	OGDF_ASSERT(G.representsCombEmbedding());
	OGDF_ASSERT(CconnectClusterPlanarityModule().isClusterPlanar(CG));

	// the augmentation edges make the instance c-connected c-plane, fixing the embedding for ClusterPlanarizationLayout
	ClusterPlanarizationLayout cpl;
	ClusterGraphAttributes CGA(CG, ClusterGraphAttributes::all);
	cpl.call(G, CGA, CG);

	// now hide the augmentation edges from the generated drawing
	Graph::HiddenEdgeSet hes(G);
	for (edge e : added) {
		hes.hide(e);
	}

	GraphIO::write(CGA, "clusters.svg");

	return 0;
}
