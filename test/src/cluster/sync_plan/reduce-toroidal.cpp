#include <ogdf/basic/Graph_d.h>
#include <ogdf/basic/Logger.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/external/coin.h>
#include <ogdf/fileformats/GraphIO.h>

#include <memory>
#include <ostream>

#include "PQPlanarity.h"
#include "return.h"
#include "utils/Clusters.h"
#include "utils/Preprocess.h"


static std::string progname;

void usage() {
	Logger::slout(Logger::Level::Alarm)
			<< "Usage: " << progname << " input [output]\n\n"
			<< "-l[0-5]  set log level\n"
			<< "-d  write out intermediate instances (only in debug builds)\n"
			<< "-c  disallow contracting independent bicon-bicon pipes\n"
			<< "-i  do not intersect trees when propagating bicon-bicon pipes\n"
			<< "-r  process pipes in a random order\n"
			<< "-b  prefer operation contract\n"
			<< "-a  sort pipes by increasing, not decreasing degree\n"
			<< "-s  swap preference of contract and other operations\n"
			<< "-p  allow batch SPQR embedding tree generation\n"
			<< std::endl;
}

bool stillMatches(const ClusterGraph& CG, const Graph& G, const PQPlanConf& conf) {
	if (G.numberOfNodes() < 5 || G.numberOfEdges() < 3 || CG.numberOfClusters() < 2) {
		return false;
	}
	if (!isPlanar(G)) {
		return false;
	}
	if (!isParallelFree(G)) {
		return false;
	}
	if (!isLoopFree(G)) {
		return false;
	}

	Logger::Level level = Logger::globalLogLevel();
	Logger::globalLogLevel(Logger::Level::Alarm);

	Graph G2;
	ClusterGraph CG2(CG, G2);
	PQPlanarity pq(&G2, &CG2);
	pq.matchings.setPipeQueue(conf.getOrder(&pq));
	pq.setAllowContractBBPipe(conf.allow_contract);
	pq.setIntersectTrees(conf.intersect_trees);
	pq.setBatchSpqr(conf.batch_spqr);

	// pq.stats_out.open("pqplan_stats.json");
	// if (!pq.stats_out.is_open() || !pq.stats_out.good()) {
	// 	Logger::slout(Logger::Level::Alarm)
	// 			<< "IO Warning: Could not open pqplan_stats.json for writing!" << std::endl;
	// } else {
	// 	pq.stats_out << "[";
	// }

	bool matches = pq.makeReduced() && pq.simplifyToroidalCycleLength > 0;
	Logger::globalLogLevel(level);

	if (matches) {
		return true;
	} else {
		return false;
	}

	// if (!pq.stats_out.is_open() || !pq.stats_out.good()) {
	// 	Logger::slout(Logger::Level::Alarm)
	// 			<< "IO Warning: Could not finish writing to pqplan_stats.json!" << std::endl;
	// } else {
	// 	pq.stats_out << "]" << std::endl;
	// }
}

template<typename C>
bool tryGraph(C callback, unique_ptr<ClusterGraph>& CG, unique_ptr<Graph>& G, const PQPlanConf& conf) {
	ClusterArray<cluster> originalClusterTable(*CG);
	NodeArray<node> originalNodeTable(*G);
	unique_ptr<Graph> G2 = make_unique<Graph>();
	unique_ptr<ClusterGraph> CG2 =
			make_unique<ClusterGraph>(*CG, *G2, originalClusterTable, originalNodeTable);

	callback(*CG2, *G2, originalClusterTable, originalNodeTable);
	makeParallelFreeUndirected(*G2);
	preprocessClusterGraph(*CG2, *G2);

	if (stillMatches(*CG2, *G2, conf)) {
		swap(G, G2);
		swap(CG, CG2);
		return true;
	} else {
		return false;
	}
}

bool tryReduce(unique_ptr<ClusterGraph>& CG, unique_ptr<Graph>& G, const PQPlanConf& conf) {
	for (cluster c : CG->clusters) {
		if (c == CG->rootCluster()) {
			continue;
		}
		cout << "\t\tContract cluster " << c->index() << endl;
		if (tryGraph(
					[&c](ClusterGraph& CG2, Graph& G2, ClusterArray<cluster>& oCT,
							NodeArray<node>& oNT) -> void {
						List<node> l;
						oCT[c]->getClusterNodes(l);
						CG2.collapse(l, G2);
					},
					CG, G, conf)) {
			return true;
		}
	}
	cout << "\tFailed to contract further clusters" << endl;

	for (cluster c : CG->clusters) {
		if (c == CG->rootCluster()) {
			continue;
		}
		cout << "\t\tDelete cluster " << c->index() << endl;
		if (tryGraph([&c](ClusterGraph& CG2, Graph& G2, ClusterArray<cluster>& oCT,
							 NodeArray<node>& oNT) -> void { CG2.delCluster(oCT[c]); },
					CG, G, conf)) {
			return true;
		}
	}
	cout << "\tFailed to delete further clusters" << endl;

	for (node n : G->nodes) {
		cout << "\t\tDelete node " << n->index() << endl;
		if (tryGraph([&n](ClusterGraph& CG2, Graph& G2, ClusterArray<cluster>& oCT,
							 NodeArray<node>& oNT) -> void { G2.delNode(oNT[n]); },
					CG, G, conf)) {
			return true;
		}
	}
	cout << "\tFailed to delete further nodes" << endl;

	for (edge e : G->edges) {
		cout << "\t\tContract edge " << e->index() << " (" << e->source()->index() << ", "
			 << e->target()->index() << ")" << endl;
		if (tryGraph(
					[&e](ClusterGraph& CG2, Graph& G2, ClusterArray<cluster>& oCT,
							NodeArray<node>& oNT) -> void {
						G2.contract(G2.searchEdge(oNT[e->source()], oNT[e->target()]));
					},
					CG, G, conf)) {
			return true;
		}
	}
	cout << "\tFailed to contract further edges" << endl;

	for (edge e : G->edges) {
		cout << "\t\tDelete edge " << e->index() << " (" << e->source()->index() << ", "
			 << e->target()->index() << ")" << endl;
		if (tryGraph(
					[&e](ClusterGraph& CG2, Graph& G2, ClusterArray<cluster>& oCT,
							NodeArray<node>& oNT) -> void {
						G2.delEdge(G2.searchEdge(oNT[e->source()], oNT[e->target()]));
					},
					CG, G, conf)) {
			return true;
		}
	}
	cout << "\tFailed to delete further edges, reached reduced instance!" << endl;

	return false;
}

int main(int argc, char* argv[]) {
	CoinManager::CoinLog.localLogMode(Logger::LogMode::Log);
	CoinManager::CoinLog.localLogLevel(Logger::Level::Alarm);
	progname = argv[0];
	Logger::Level log_level = Logger::Level::Minor;
	PQPlanConf pqconf;

	int opt;
	while ((opt = getopt(argc, argv, "l:dm:t:fcirbasp")) != -1) {
		switch (opt) {
		case 'l':
			log_level = parseLogLevel();
			break;
#ifdef OGDF_DEBUG
		case 'd':
			PQPlanarityConsistency::doWriteOut = true;
			break;
#endif
#define TOGGLE(key, var)          \
	case key:                     \
		pqconf.var = !pqconf.var; \
		break;
			PQPlanConf_KEYS
#undef TOGGLE
					default : /* '?' */
							  usage();
			return ERROR_OPTIONS;
		}
	}
	if (optind >= argc) {
		Logger::slout(Logger::Level::Alarm) << "Invalid Options: Expected infile!" << std::endl;
		usage();
		return ERROR_OPTIONS;
	}
	string infile = argv[optind++];
	string outfile;
	if (optind < argc) {
		outfile = argv[optind++];
	} else {
		outfile = "reduced-" + infile;
	}
	if (optind < argc) {
		Logger::slout(Logger::Level::Alarm)
				<< "Invalid Options: Too many positional arguments!" << std::endl;
		usage();
		return ERROR_OPTIONS;
	}

	Graph G;
	ClusterGraph CG(G);
	int read_ret = readCG(infile, log_level, CG, G);
	if (read_ret != SUCCESS) {
		return read_ret;
	}

	// for (edge e : G.edges) {
	// 	if (e->source()->index() > e->target()->index()) {
	// 		G.reverseEdge(e);
	// 	}
	// }

	Logger::slout() << "PQPlanConf: " << pqconf << std::endl;
	preprocessLog.localLogLevel(ogdf::Logger::Level::Alarm);
	preprocessLog.localLogMode(ogdf::Logger::LogMode::Log);
	if (!stillMatches(CG, G, pqconf)) {
		return FAILURE;
	}

	unique_ptr<Graph> Gr = make_unique<Graph>();
	unique_ptr<ClusterGraph> CGr = make_unique<ClusterGraph>(CG, *Gr);
	for (int i = 0; true; i++) {
		cout << "Iteration " << i << endl;
		if (!tryReduce(CGr, Gr, pqconf)) {
			break;
		}
	}

	if (!GraphIO::write(*CGr, outfile)) {
		Logger::slout(Logger::Level::Alarm) << "IO Error: Couldn't write graph!" << std::endl;
		return ERROR_IO;
	}
	return SUCCESS;
}
