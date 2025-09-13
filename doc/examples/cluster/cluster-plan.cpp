#include <ogdf/basic/GraphSets.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/graph_generators/clustering.h>
#include <ogdf/cluster/ClusterGraph.h>
#include <ogdf/cluster/ClusterGraphAttributes.h>
#include <ogdf/cluster/ClusterPlanarizationLayout.h>
#include <ogdf/cluster/sync_plan/ClusterPlanarity.h>
#include <ogdf/fileformats/GraphIO.h>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace ogdf;

int main(void) {
	// create a random cluster graph
	Graph G;
	ClusterGraph CG(G);
	ClusterGraphAttributes CGA(CG, ClusterGraphAttributes::all);
	randomClusterPlanarGraph(G, CG, 10, 20, 50);

	// set up the sync plan test
	SyncPlanClusterPlanarityModule sp;
	std::vector<std::pair<adjEntry, adjEntry>> augmentation;
	sp.setStoreAugmentation(&augmentation);

	// run the test + embedder
	if (!sp.clusterPlanarEmbed(CG, G)) {
		throw std::runtime_error("Not c-planar!");
	}

	// add the computed augmentation edges
	EdgeSet added(G);
	insertAugmentationEdges(CG, G, augmentation, &added);

	// the augmentation edges make the instance c-connected c-plane, fixing the embedding for ClusterPlanarizationLayout
	ClusterPlanarizationLayout cpl;
	cpl.call(G, CGA, CG);

	// now hide the augmentation edges from the generated drawing
	Graph::HiddenEdgeSet hes(G);
	for (edge e : added) {
		hes.hide(e);
	}

	GraphIO::write(CGA, "clusters.svg");
	return 0;
}
