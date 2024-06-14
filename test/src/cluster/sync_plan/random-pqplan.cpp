#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/graphics.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/energybased/FMMMLayout.h>
#include <ogdf/fileformats/GraphIO.h>

#include <filesystem>
#include <memory>
#include <new>

#include <getopt.h>

#include "return.h"
#include "utils/Logging.h"
#include "utils/Random.h"

using namespace ogdf;
using nlohmann::json;
static std::string progname;

void usage() {
	Logger::slout(Logger::Level::Alarm)
			<< "Usage: " << progname << " [-options] input [outfile] [svgfile]\n\n"
			<< "-f generate even if outfile already exists\n"
			<< "-l[0-5] set log level\n"
			<< "-s[int] set seed for random number generator or 'r' to use a true random seed\n"

			<< "-a[int|double] first parameter\n"
			<< "-b[int|double] second parameter\n"
			<< "-m[int] method to use:\n"
			<< "   1: PQPlan (int pipe_count, int min_deg = 3)\n"
			<< "   2: SEFE by shared graph(int edges1, int edges2)\n"
			<< "   3: SEFE by union graph(double frac_shared, double frac_g1)\n\n"
			<< "outfile can use the following templates that will be replaced based on infile and the given options:\n"
			<< "%s | tag (seed+m+a+b) | -s123-m1-a10-b4\n"
			<< "%d | infile directory | graphs/\n"
			<< "%n | infile basename  | graph01\n"
			<< "%e | infile extension | .gml" << std::endl;
}

inline string format(string& in, const string& f, const string& r) {
	const auto idx = in.find(f);
	if (idx != string::npos) {
		in.replace(idx, f.length(), r);
	}
	return in;
}

int main(int argc, char* argv[]) {
	progname = argv[0];
	std::string outfile, svgfile;
	Logger::Level log_level = Logger::Level::Minor;
	bool force = false;
	int seed = 0, mode = 0;
	double a = 0, b = 0;

	int opt;
	while ((opt = getopt(argc, argv, "fl:s:m:a:b:")) != -1) {
		switch (opt) {
		case 'f':
			force = true;
			break;
		case 'l':
			log_level = parseLogLevel();
			break;
		case 's':
			if (string("r") == optarg) {
				std::random_device rd;
				seed = rd();
			} else {
				seed = std::stoi(optarg);
			}
			break;
		case 'm':
			mode = std::stoi(optarg);
			break;
		case 'a':
			a = std::stod(optarg);
			break;
		case 'b':
			b = std::stod(optarg);
			break;
		default: /* '?' */
			usage();
			return ERROR_OPTIONS;
		}
	}

	if (optind >= argc) {
		Logger::slout(Logger::Level::Alarm) << "Invalid Options: Expected infile!" << std::endl;
		usage();
		return ERROR_OPTIONS;
	}
	std::string infile = argv[optind++];
	if (optind < argc) {
		outfile = argv[optind++];
	}
	if (optind < argc) {
		svgfile = argv[optind++];
	}
	if (optind < argc) {
		Logger::slout(Logger::Level::Alarm)
				<< "Invalid Options: Too many positional arguments!" << std::endl;
		usage();
		return ERROR_OPTIONS;
	}
	if (mode < 1 || mode > 3) {
		Logger::slout(Logger::Level::Alarm)
				<< "Invalid Options: Mode must be in range [1, 2, 3]!" << std::endl;
		usage();
		return ERROR_OPTIONS;
	}

	Graph G;
	Logger::globalLogLevel(Logger::Level::Alarm);
	bool read_success = GraphIO::read(G, infile);
	Logger::globalLogLevel(log_level);
	if (!read_success) {
		std::ifstream is(infile);
		if (!is.good()) {
			Logger::slout(Logger::Level::Alarm)
					<< "IO Error: Couldn't read file " << infile << "!" << std::endl;
		} else {
			Logger::slout(Logger::Level::Alarm)
					<< "IO Error: Couldn't read graph from file " << infile << "!" << std::endl;
		}
		return ERROR_IO;
	}

	const std::filesystem::path inpath = infile;
	std::stringstream ss;
	ss << "-s" << seed << "-m" << mode << "-a" << a << "-b" << b;
	const string& tag = ss.str();
	format(outfile, "%s", tag);
	format(outfile, "%d", inpath.parent_path());
	format(outfile, "%n", inpath.stem());
	format(outfile, "%e", inpath.extension());
	if (!force) {
		if (exists(std::filesystem::path(outfile))) {
			Logger::slout(Logger::Level::High) << "Skipping input " << infile << " as outfile "
											   << outfile << " exists" << std::endl;
			return SUCCESS;
		}
		if (exists(std::filesystem::path(outfile + ".ignore"))) {
			Logger::slout(Logger::Level::High) << "Skipping input " << infile << " as outfile "
											   << outfile << ".ignore exists" << std::endl;
			return SUCCESS;
		}
	}

	Logger::slout(Logger::Level::High) << "Random seed is " << seed << std::endl;
	setSeed(seed);
	if (mode != 1 && !isConnected(G)) {
		Logger::slout(Logger::Level::Alarm)
				<< "Error: The given graph is not connected!" << std::endl;
		return NOT_PQPLANAR;
	}
	if (!planarEmbed(G)) {
		Logger::slout(Logger::Level::Alarm) << "Error: The given graph is not planar!" << std::endl;
		return NOT_PQPLANAR;
	}
	makeSimpleUndirected(G);
	Logger::slout(Logger::Level::Medium) << "Input " << G << "." << endl;

	// 1: PQPlan (int pipe_count, int min_deg = 3)
	// 2: SEFE by shared graph(int edges1, int edges2)
	// 3: SEFE by union graph(double frac_shared, double frac_g1)
	unique_ptr<Graph> work;
	unique_ptr<PQPlanarity> pq;
	EdgeArray<uint8_t> edge_types;
	if (mode == 1) {
		Logger::slout(Logger::Level::High) << "Calling randomPQPlanInstance(pipe_count=" << (int)a
										   << ", min_deg=" << (int)b << ")" << std::endl;
		if ((int)a <= 0 || (int)b <= 0) {
			Logger::slout(Logger::Level::Alarm)
					<< "Invalid Options: a and b must be positive integers!" << std::endl;
			usage();
			return ERROR_OPTIONS;
		}
		pq = make_unique<PQPlanarity>(&G);
		randomPQPlanInstance(*pq, (int)a, (int)b);
	} else {
		edge_types.init(G, 0);
		if (mode == 2) {
			Logger::slout(Logger::Level::High)
					<< "Calling randomSEFEInstanceBySharedGraph(edges1=" << (int)a
					<< ", edges2=" << (int)b << ")" << std::endl;
			if ((int)a <= 0 || (int)b <= 0) {
				Logger::slout(Logger::Level::Alarm)
						<< "Invalid Options: a and b must be positive integers!" << std::endl;
				usage();
				return ERROR_OPTIONS;
			}
			randomSEFEInstanceBySharedGraph(&G, edge_types, (int)a, (int)b);
		} else {
			Logger::slout(Logger::Level::High)
					<< "Calling randomSEFEInstanceByUnionGraph(frac_shared=" << a
					<< ", frac_g1=" << b << ")" << std::endl;
			if (a <= 0 || b <= 0 || a + b >= 1) {
				Logger::slout(Logger::Level::Alarm)
						<< "Invalid Options: a and b must be within range (0, 1) such that a + b < 1!"
						<< std::endl;
				usage();
				return ERROR_OPTIONS;
			}
			randomSEFEInstanceByUnionGraph(&G, edge_types, a, b);
		}
		work = make_unique<Graph>();
		pq = make_unique<PQPlanarity>(&G, work.get(), edge_types);
	}
	int max_deg = 0, tot_deg = 0;
	for (auto& p : pq->matchings.getPipes()) {
		tot_deg += p.degree();
		if (p.degree() > max_deg) {
			max_deg = p.degree();
		}
	}
	Logger::slout(Logger::Level::Medium) << "Output " << *pq << "." << endl;

	if (!outfile.empty()) {
		Logger::slout(Logger::Level::High) << "Writing to " << outfile << std::endl;
		GraphAttributes GA(*pq->G, GraphAttributes::nodeLabel | GraphAttributes::edgeLabel);
		for (node n : pq->G->nodes) {
			GA.label(n) = to_string(n->index());
		}
		for (edge e : pq->G->edges) {
			GA.label(e) = to_string(e->index());
		}
		if (!GraphIO::write(GA, outfile)) {
			Logger::slout(Logger::Level::Alarm)
					<< "IO Error: Couldn't write graph to " << outfile << "!" << std::endl;
			return ERROR_IO;
		}
		json config;
		PQPlanOptions::generateConfigJSON<std::function<string(node)>, std::function<string(edge)>>(
				*pq, config, [&GA](node n) -> string { return GA.label(n); },
				[&GA](edge e) -> string { return GA.label(e); });
		std::ofstream out(outfile + ".json");
		out << config;
		if (!out.good()) {
			Logger::slout(Logger::Level::Alarm)
					<< "IO Error: Couldn't write pipes to " << outfile << ".json!" << std::endl;
			return ERROR_IO;
		}
	}

	if (!svgfile.empty()) {
		Logger::slout(Logger::Level::High) << "Drawing to " << svgfile << std::endl;
		if (mode == 1) {
			PQPlanarityDrawer draw(pq.get());
			GraphAttributes& GA = draw.ensureGraphAttributes();
			draw.layout();
			if (!GraphIO::write(GA, svgfile)) {
				Logger::slout(Logger::Level::Alarm)
						<< "IO Error: Couldn't draw graph to " << svgfile << "!" << std::endl;
				return ERROR_IO;
			}
		} else {
			GraphAttributes GA(G, GraphAttributes::all);
			for (edge e : G.edges) {
				switch (edge_types[e]) {
				case 1:
					GA.strokeColor(e) = Color("#0F0");
					break;
				case 2:
					GA.strokeColor(e) = Color("#00F");
					break;
				case 3:
					GA.strokeColor(e) = Color("#000");
					break;
				default:
					GA.strokeColor(e) = Color("#F00");
					break;
				}
			}
			FMMMLayout fmmm;
			fmmm.useHighLevelOptions(true);
			fmmm.unitEdgeLength(50);
			fmmm.newInitialPlacement(true);
			fmmm.qualityVersusSpeed(FMMMOptions::QualityVsSpeed::GorgeousAndEfficient);
			fmmm.call(GA);
			if (!GraphIO::write(GA, svgfile)) {
				Logger::slout(Logger::Level::Alarm)
						<< "IO Error: Couldn't draw graph to " << svgfile << "!" << std::endl;
				return ERROR_IO;
			}
		}
	}

	// if (!(pq->makeReduced() && pq->solveReduced())) {
	//     Logger::slout(Logger::Level::Alarm) << "Instance is not PQ-Planar!" << std::endl;
	//     return ERROR_COMB_EMB;
	// }

	return SUCCESS;
}
