#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/cluster/ClusterGraph.h>
#include <ogdf/cluster/sync_plan/utils/Clusters.h>
#include <ogdf/fileformats/GraphIO.h>
// #include <ogdf/planarlayout/FPPLayout.h>

using namespace ogdf;
using namespace std;

using tpc = chrono::high_resolution_clock;
using tp = chrono::time_point<chrono::high_resolution_clock>;


static Logger l;

class Clusterer {
	ClusterGraph& CG;
	RandomClusterConfig& config;

	GraphCopy copy;
	NodeArray<cluster> clusters;
	NodeArray<bool> mark;
	int marked = 0;
	tp stop;

	// FPPLayout layout;
	// GraphAttributes minor_GA;
	// int dump_nr = 0;

public:
	Clusterer(ClusterGraph& CG, RandomClusterConfig& config)
		: CG(CG), config(config), copy(CG.constGraph()), clusters(copy, nullptr), mark(copy, false) {
		OGDF_ASSERT(CG.getGraph()->representsCombEmbedding());
		OGDF_ASSERT(CG.numberOfClusters() == 1);
		CG.adjAvailable(true);
		copy.setOriginalEmbedding();
		if (!isConnected(copy)) {
			makeConnected(copy);
		}
		triangulate(copy);
		OGDF_ASSERT(copy.representsCombEmbedding());

		if (config.timeout != 0) {
			stop = tpc::now() + chrono::seconds(config.timeout);
		} else {
			stop = tp::max();
		}

		// {
		//     minor_GA.init(copy, GraphAttributes::all ^ GraphAttributes::nodeId);
		//     layout.separation(30);
		//     for (node n : copy.nodes) { minor_GA.label(n) = to_string(n->index()); }
		//     for (edge e : copy.edges) {
		//         minor_GA.label(e) = to_string(e->index());
		//         if (!copy.original(e)) { minor_GA.strokeColor(e) = Color(Color::Name::Gray); }
		//     }
		//     layout.callFixEmbed(minor_GA);
		//     dump();
		// }
	}

	bool timedout() const { return tpc::now() > stop; }

	bool hasFreeNodes() const {
		return !timedout() && CG.rootCluster()->nCount() > config.min_root_nodes;
	}

	bool shouldAddCluster() const {
		return hasFreeNodes()
				&& (config.max_clusters == 0 || CG.numberOfClusters() < config.max_clusters)
				&& randomDouble(0, 1) > config.prob_no_further_cluster;
	}

	bool shouldAddNode(cluster c) const {
		return hasFreeNodes()
				&& (config.max_nodes_in_cluster == 0 || c->nCount() < config.max_nodes_in_cluster)
				&& randomDouble(0, 1) > config.prob_no_further_node;
	}

	bool isTriangulated() {
		for (node n : copy.nodes) {
			for (adjEntry adj : n->adjEntries) {
				adjEntry succ2 = adj->faceCycleSucc()->faceCycleSucc();
				if (succ2 != adj && succ2->faceCycleSucc() != adj) {
					return false;
				}
			}
		}
		return true;
	}

	void makeClusters() {
		while (shouldAddCluster()) {
			node n = copy.chooseNode([](node n) -> bool { return n->degree() > 0; });
			if (!n) {
				break;
			}
			cluster c = CG.createEmptyCluster();
			mergeIntoCluster(n, c);
			l.lout() << "Created cluster " << c->index() << " with node " << n->index() << "."
					 << std::endl;
			Logger::Indent _(l);

			OGDF_ASSERT(!mark[n]);
			for (adjEntry adj : n->adjEntries) {
				if (!mark[adj->twinNode()]) {
					mark[adj->twinNode()] = true;
					marked++;
				}
			}
			// minor_GA.fillColor(n) = Color(Color::Name::Red);

			while (shouldAddNode(c) && n->degree() > 0) {
				adjEntry rand_adj = nullptr;
				adjEntry adj = n->adjEntries.head();
				int pos = randomNumber(0, n->degree() - 1);
				if (pos < n->degree() / 2) {
					for (int i = 0; i < pos; ++i) {
						adj = adj->cyclicSucc();
					}
				} else {
					for (int i = 0; i < n->degree() - pos; ++i) {
						adj = adj->cyclicPred();
					}
				}
				for (int i = 0; i < n->degree(); ++i) {
					if (!canContract(adj)) {
						adj = adj->cyclicSucc();
					} else {
						rand_adj = adj;
						break;
					}
				}
				if (!rand_adj) {
					break;
				}

				// minor_GA.strokeColor(rand_adj->theEdge()) = Color(Color::Name::Red);
				// minor_GA.strokeWidth(rand_adj->theEdge()) = 3;
				dump();

				OGDF_ASSERT(canContract(rand_adj));
				contractIntoCluster(rand_adj);

				OGDF_ASSERT(copy.representsCombEmbedding());
				node c1, c2;
				OGDF_ASSERT(isTriconnected(copy, c1, c2));
				OGDF_ASSERT(isTriangulated());
			}

			OGDF_ASSERT(!mark[n]);
			for (adjEntry adj : n->adjEntries) {
				if (mark[adj->twinNode()]) {
					mark[adj->twinNode()] = false;
					marked--;
				}
			}
			OGDF_ASSERT(marked == 0);
			// minor_GA.fillColor(n) = Color(Color::Name::Wheat);

			// if (c->nCount() + c->cCount() <= 1) {
			//     clusters[n] = c->parent();
			//     CG.delCluster(c);
			// } else {
			makeClusterAdjs(n);
			// }
		}
		OGDF_ASSERT(isClusterPlanarEmbedding(CG));
	}

	bool canContract(adjEntry adj) {
		if (config.cconnected && !copy.original(adj->theEdge())) {
			return false;
		}
		node neigh1 = nullptr, neigh2 = nullptr;
		for (adjEntry a : adj->twinNode()->adjEntries) {
			if (!mark[a->twinNode()]) {
				continue;
			}
			if (a->twinNode() == neigh1) {
				continue;
			}
			if (a->twinNode() == neigh2) {
				continue;
			}

			if (neigh1 == nullptr) {
				neigh1 = a->twinNode();
			} else if (neigh2 == nullptr) {
				neigh2 = a->twinNode();
			} else {
				return false;
			}
		}
		return true;
	}

	void contractIntoCluster(adjEntry adj) {
		edge to_contract = adj->theEdge();
		node u = adj->theNode();
		node v = adj->twinNode();
		OGDF_ASSERT(!to_contract->isSelfLoop());

		l.lout() << "Contracting node " << v << " (°" << v->degree() << ") into the cluster of "
				 << u << " (°" << u->degree() << ")" << std::endl; // << minor_GA.label(u)
		// minor_GA.label(u).append(" ").append(minor_GA.label(v));
		mergeIntoCluster(v, clusters[u]);

		OGDF_ASSERT(mark[v] != false);
		mark[v] = false;
		marked--;

		for (adjEntry a : v->adjEntries) {
			if (a->twinNode() != u && !mark[a->twinNode()]) {
				mark[a->twinNode()] = true;
				marked++;
			}
		}
		if (!adj->isSource()) {
			copy.reverseEdge(adj->theEdge());
		}
		node r = copy.contract(adj->theEdge());
		OGDF_ASSERT(r == u);
#ifdef OGDF_DEBUG
		CG.consistencyCheck();
#endif
	}

	void mergeIntoCluster(node n, cluster c) {
		if (clusters[n] != nullptr) {
			if (clusters[n] == c) {
				return;
			}
			OGDF_ASSERT(clusters[n] != nullptr);
			CG.moveCluster(clusters[n], c);
		} else {
			CG.reassignNode(copy.original(n), c);
		}
		clusters[n] = c;
	}

	void makeClusterAdjs(node n) {
		std::ostream& ls = l.lout();
		ls << "Cluster border crossing adjacencies: ";
		cluster c = clusters[n];
		c->adjEntries.clear();
		for (adjEntry adj : n->adjEntries) {
			edge Ge = copy.original(adj->theEdge());
			if (!Ge) {
				continue;
			}
			if (c->isDescendant(CG.clusterOf(Ge->source()), true)) {
				OGDF_ASSERT(!c->isDescendant(CG.clusterOf(Ge->target()), true));
				c->adjEntries.pushBack(Ge->adjSource());
				ls << Ge->source()->index() << " --" << Ge->index() << "-> "
				   << Ge->target()->index() << ", ";
			} else {
				OGDF_ASSERT(c->isDescendant(CG.clusterOf(Ge->target()), true));
				c->adjEntries.pushBack(Ge->adjTarget());
				ls << Ge->target()->index() << " <-" << Ge->index() << "-- "
				   << Ge->source()->index() << ", ";
			}
		}
		ls << std::endl;
		CG.adjAvailable(true);
		OGDF_HEAVY_ASSERT(isClusterPlanarEmbedding(CG));
	}

	void dump() {
		// l.lout() << "Writing dump " << dump_nr << std::endl;
		// // layout.callFixEmbed(minor_GA);
		// {
		//     std::stringstream ss;
		//     ss << "random-cplan" << dump_nr << ".svg";
		//     GraphIO::write(minor_GA, ss.str());
		// }
		// {
		//     std::stringstream ss;
		//     ss << "random-cplan" << dump_nr << ".gml";
		//     GraphIO::write(minor_GA, ss.str());
		// }
		// dump_nr++;
	}
};

bool makeClusters(ClusterGraph& CG, RandomClusterConfig& config) {
	Clusterer c(CG, config);
	c.makeClusters();
	return !c.timedout();
}

bool isClusterPlanarEmbedding(const ClusterGraph& CG) {
	if (!CG.adjAvailable()) {
		return false;
	}
#ifdef OGDF_DEBUG
	CG.consistencyCheck();
#endif
	GraphCopy gcopy(CG.constGraph());
	gcopy.setOriginalEmbedding();
	clusterBorderToEdges(CG, gcopy, nullptr, [&gcopy](edge e) -> edge { return gcopy.copy(e); });
	bool copy_valid = gcopy.representsCombEmbedding();
	//bool cluster_valid = CG.representsCombEmbedding();
	//OGDF_ASSERT(copy_valid == cluster_valid);
	return copy_valid;
}

void clusterBorderToEdges(const ClusterGraph& CG, Graph& G,
		EdgeArray<List<std::pair<adjEntry, cluster>>>* subdivisions,
		const std::function<edge(edge)>& translate) {
	for (cluster c = CG.firstPostOrderCluster(); c != CG.rootCluster(); c = c->pSucc()) {
		adjEntry prev_ray = nullptr, first_ray = nullptr;
		for (adjEntry adj : c->adjEntries) {
			bool reverse = adj->isSource();
			edge the_edge = translate(adj->theEdge());
			if (reverse) {
				G.reverseEdge(the_edge);
			}
			edge new_edge = G.split(the_edge);
			adjEntry spike_to_src =
					the_edge->adjTarget(); // adjacency of the new node toward the source of the_edge
			if (subdivisions != nullptr) {
				if (reverse) {
					(*subdivisions)[adj->theEdge()].emplaceBack(spike_to_src, c);
				} else {
					(*subdivisions)[adj->theEdge()].emplaceFront(spike_to_src, c);
				}
			}
			if (reverse) {
				G.reverseEdge(the_edge);
				G.reverseEdge(new_edge);
			}

			if (prev_ray != nullptr) {
				G.newEdge(prev_ray, spike_to_src->cyclicPred());
			}
			prev_ray = spike_to_src;
			if (first_ray == nullptr) {
				first_ray = spike_to_src;
			}
		}
		if (first_ray != nullptr) {
			G.newEdge(prev_ray, first_ray->cyclicPred());
		}
	}
}

std::ostream& printClusters(cluster c, std::ostream& s) {
	s << c->nCount() << "[";
	for (cluster child : c->children) {
		printClusters(child, s);
	}
	return s << "]";
}

void printCG(const ClusterGraph& CG, const string& type) {
	Logger::slout(Logger::Level::High)
			<< type << "ClusterGraph with " << CG.constGraph().numberOfNodes() << " nodes, "
			<< CG.constGraph().numberOfEdges() << " edges, " << CG.numberOfClusters()
			<< " clusters with max depth " << CG.treeDepth() << ". "
			<< (isClusterPlanarEmbedding(CG)
							   ? "Cluster-"
							   : (CG.constGraph().representsCombEmbedding() ? "" : "Non-"))
			<< "Planar embedding." << std::endl;
	if (CG.treeDepth() < 2) {
		Logger::slout(Logger::Level::Alarm)
				<< "Warning: " << type << "Graph contains no clusters!" << std::endl;
	}
	printClusters(CG.rootCluster(), Logger::slout()) << std::endl;
}
