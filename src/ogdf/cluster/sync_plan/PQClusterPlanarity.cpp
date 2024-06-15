#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/cluster/sync_plan/PQPlanarity.h>
#include <ogdf/cluster/sync_plan/basic/GraphUtils.h>
#include <ogdf/cluster/sync_plan/utils/Clusters.h>
#include <ogdf/cluster/sync_plan/utils/Logging.h>

struct FrozenCluster {
	int index = -1, parent = -1, parent_node = -1;
	List<int> nodes;

	FrozenCluster(int index, int parent) : index(index), parent(parent) { }
};

class UndoInitCluster : public PQPlanarity::UndoOperation {
public:
	ClusterGraph* cg;
	List<FrozenCluster> clusters;

	explicit UndoInitCluster(ClusterGraph* cg) : cg(cg) {
		for (cluster c = cg->firstPostOrderCluster(); c != nullptr; c = c->pSucc()) {
			auto it = clusters.emplaceFront(c->index(),
					c->parent() != nullptr ? c->parent()->index() : -1);
			for (node n : c->nodes) {
				(*it).nodes.pushBack(n->index());
			}
		}
	}

	void processCluster(PQPlanarity& pq, cluster c, int parent_node) {
		node n = pq.nodeFromIndex(parent_node);
		node t = pq.matchings.getTwin(n);
		pq.log.lout(Logger::Level::Medium)
				<< "Processing cluster " << c << " with node " << pq.fmtPQNode(n, false)
				<< " matched with " << pq.fmtPQNode(t, false) << " in the parent cluster "
				<< c->parent() << std::endl;
		Logger::Indent _(&pq.log);
		pq.log.lout(Logger::Level::Minor) << c->nodes << std::endl;

		PipeBij bij;
		pq.matchings.getIncidentEdgeBijection(t, bij);
		pq.log.lout(Logger::Level::Minor) << printBijection(bij) << std::endl;
		pq.matchings.removeMatching(n, t);
		join(*pq.G, t, n, bij);
		pq.log.lout(Logger::Level::Minor) << printEdges(bij) << std::endl;

		for (const PipeBijPair& pair : bij) {
			c->adjEntries.pushBack(pair.first->twin());
		}
	}

	void undo(PQPlanarity& pq) override {
		OGDF_ASSERT(cg->numberOfClusters() == 1);
		cg->rootCluster()->adjEntries.clear();
		ClusterArray<cluster> cluster_index(*cg, nullptr);
		for (FrozenCluster& fc : clusters) {
			cluster c;
			if (fc.index == cg->rootCluster()->index()) {
				OGDF_ASSERT(fc.parent == -1);
				c = cg->rootCluster();
			} else {
				OGDF_ASSERT(fc.parent >= 0);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
				OGDF_ASSERT(cluster_index[fc.parent] != nullptr);
				OGDF_ASSERT(cluster_index[fc.parent]->index() == fc.parent);
				c = cg->newCluster(cluster_index[fc.parent], fc.index);
#pragma GCC diagnostic pop
				processCluster(pq, c, fc.parent_node);
			}
			cluster_index[c] = c;
			for (int n : fc.nodes) {
				cg->reassignNode(pq.nodeFromIndex(n), c);
			}
		}
		cg->adjAvailable(true);
#ifdef OGDF_DEBUG
		for (cluster c = cg->firstPostOrderCluster(); c != cg->rootCluster(); c = c->pSucc()) {
			for (adjEntry adj : c->adjEntries) {
				OGDF_ASSERT(pq.edge_reg[adj] == adj->theEdge());
			}
		}
		cg->constGraph().consistencyCheck();
		cg->consistencyCheck();
		// OGDF_ASSERT(cg->representsCombEmbedding());
		OGDF_ASSERT(isClusterPlanarEmbedding(*cg));
#endif
	}

	ostream& print(ostream& os) const override { return os << "UndoInitCluster"; }
};

PQPlanarity::PQPlanarity(Graph* g, ClusterGraph* cg, ClusterGraphAttributes* cga)
	: G(g)
	, GA(cga)
	, matchings(G)
	, partitions(G)
	, components(G)
	, is_wheel(*G, false)
#ifdef OGDF_DEBUG
	, consistency(*this)
#endif
{
	OGDF_ASSERT(cg->getGraph() == g);
	undo_stack.pushBack(new ResetIndices(*this));

	auto* op = new UndoInitCluster(cg);
	log.lout() << "Processing " << cg->clusters.size() << " clusters (max id "
			   << cg->maxClusterIndex() << ") from " << cg->firstPostOrderCluster()->index()
			   << " up to, but excluding root " << cg->rootCluster()->index() << "." << std::endl;
	auto fc_it = op->clusters.rbegin();
	for (cluster c = cg->firstPostOrderCluster(); c != cg->rootCluster(); c = c->pSucc()) {
		Logger::Indent _(&log);
		log.lout(Logger::Level::Medium)
				<< "Rerouting perimeter-crossing-edges of cluster " << c->index() << " with parent "
				<< c->parent()->index() << "." << std::endl;

		node cn = G->newNode();
		node pn = G->newNode();
		cg->reassignNode(cn, c);
		cg->reassignNode(pn, c->parent());

		if (GA != nullptr) {
			std::ostringstream ss;
			ss << "CN " << cn->index() << " [" << c->index() << "<" << c->parent()->index() << "]";
			GA->label(cn) = ss.str();
			GA->x(cn) = cga->x(c) + 10;
			GA->y(cn) = cga->y(c) + 10;
			ss.str("");
			ss << "PN " << pn->index() << " [" << c->parent()->index() << ">" << c->index() << "]";
			GA->label(pn) = ss.str();
			GA->x(pn) = cga->x(c) - 10;
			GA->y(pn) = cga->y(c) - 10;
		}

		OGDF_ASSERT(fc_it != op->clusters.rend());
		OGDF_ASSERT((*fc_it).index == c->index());
		(*fc_it).parent_node = cn->index();
		++fc_it;
		log.lout(Logger::Level::Minor)
				<< "Matched child node " << cn->index() << " in cluster " << c->index()
				<< " with parent node " << pn->index() << " in cluster " << c->parent()->index()
				<< ". Now processing " << (c->nodes.size() - 1) << " nodes in child cluster."
				<< std::endl;

		int count = 0;
		for (node n : c->nodes) {
			if (n == cn) {
				continue;
			}
			Logger::Indent _(&log);

			List<adjEntry> adjEntries;
			for (adjEntry adj : n->adjEntries) {
				if (c != cg->clusterOf(adj->twinNode())) {
					adjEntries.pushBack(adj);
				}
			}
			count += adjEntries.size();

			log.lout(Logger::Level::Minor)
					<< "Processing " << n->adjEntries.size() << " incident edges of node "
					<< n->index() << ", of which " << adjEntries.size()
					<< " are perimeter-crossing." << std::endl;
			for (adjEntry adj : adjEntries) {
				splitEdge(*G, adj->twin(), pn, cn);
			}
		}

		if (count == 0) {
			log.lout(Logger::Level::Minor)
					<< "Cluster " << c->index() << " has no perimeter-crossing edges." << std::endl;
		} else {
			G->reverseAdjEdges(pn);
		}

		matchings.matchNodes(cn, pn);
	}
	OGDF_ASSERT(fc_it != op->clusters.rend());
	OGDF_ASSERT((*fc_it).index == cg->rootCluster()->index());
	++fc_it;
	OGDF_ASSERT(fc_it == op->clusters.rend());

	for (cluster c = cg->firstPostOrderCluster(); c != cg->rootCluster(); c = c->pSucc()) {
		while (!c->nodes.empty()) {
			cg->reassignNode(c->nodes.front(), cg->rootCluster());
		}
		cg->delCluster(c);
	}

	initComponents();
	matchings.rebuildHeap();
	pushUndoOperationAndCheck(op);
}
