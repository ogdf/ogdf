#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/external/coin.h>
#include <ogdf/fileformats/GraphIO.h>

#include <json.hpp>

#include "PQPlanarity.h"
#include "PQPlanarityOptions.h"
#include "return.h"
#include "utils/Logging.h"

using json = nlohmann::json;

static std::string progname;

void usage() {
	Logger::slout(Logger::Level::Alarm)
			<< "Usage: " << progname << " [-p f,s]+ [-q a,b,c,...]+ [-e n:e,f,g...]+\n"
			<< "[-c conffile] [-l 0-5] [-d] input [outfile]" << std::endl;
}

// TODO same utility for (converting from) sefe
//        #define OGDF_ASSERT_RETURNS(VAL, CHECK) {auto &&tmp = VAL; OGDF_ASSERT(VAL CHECK);}
//        #define OGDF_ASSERT_RETURNS(VAL, CHECK) {VAL;}
//        OGDF_ASSERT_RETURNS(minor.contract(to_contract, true), .degree() == n);

int main(int argc, char* argv[]) {
	CoinManager::CoinLog.localLogMode(Logger::LogMode::Log);
	CoinManager::CoinLog.localLogLevel(Logger::Level::Alarm);
	progname = argv[0];
	PQPlanOptions options;
	char* conffile = nullptr;
	char* outfile = nullptr;
	Logger::Level log_level = Logger::Level::Minor;

	try {
		int opt;
		while ((opt = getopt(argc, argv, "p:q:e:c:l:d")) != -1) {
			switch (opt) {
			case 'p':
				options.parseOptionPipes(optarg);
				break;
			case 'q':
				options.parseOptionPartitions(optarg);
				break;
			case 'e':
				options.parseOptionEmbedding(optarg);
				break;
			case 'c':
				conffile = optarg;
				break;
			case 'l':
				log_level = parseLogLevel();
				break;
#ifdef OGDF_DEBUG
			case 'd':
				PQPlanarityConsistency::doWriteOut = true;
				break;
#endif
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
		char* infile = argv[optind++];
		if (optind < argc) {
			outfile = argv[optind++];
		}
		if (optind < argc) {
			Logger::slout(Logger::Level::Alarm)
					<< "Invalid Options: Too many positional arguments!" << std::endl;
			usage();
			return ERROR_OPTIONS;
		}

		Graph G;
		GraphAttributes GA(G, GraphAttributes::all);
		int read_ret = readCG(infile, log_level, GA, G);
		if (read_ret != SUCCESS) {
			return read_ret;
		}
		Logger::slout(Logger::Level::High) << G << std::endl;
		PQPlanarity pq(&G, &GA);

		Logger::slout(Logger::Level::High) << "Applying PQ-plan options." << std::endl;
		options.apply(G, GA, pq);
		if (conffile != nullptr) {
			try {
				std::ifstream i(conffile);
				json j;
				i >> j;
				PQPlanOptions::applyConfigJSON(G, GA, pq, j);
			} catch (json::parse_error& e) {
				Logger::slout(Logger::Level::Alarm)
						<< "IO Error: Couldn't parse config json " << conffile << ": " << e.what()
						<< std::endl
						<< "exception id: " << e.id << std::endl
						<< "byte position of error: " << e.byte << std::endl;
				return ERROR_IO;
			}
		}

		Logger::slout(Logger::Level::High) << pq << std::endl;
		bool solved = isPlanar(G);
		if (solved) {
			solved = pq.makeReduced();
			Logger::slout(Logger::Level::High) << pq << std::endl;
			if (solved) {
				Logger::slout(Logger::Level::High) << "REDUCED" << std::endl;
				solved = pq.solveReduced();
				Logger::slout(Logger::Level::High) << pq << std::endl;
				if (solved) {
					Logger::slout(Logger::Level::High) << "SOLVED" << std::endl;
					pq.embed();
					Logger::slout(Logger::Level::High) << pq << std::endl;
					Logger::slout(Logger::Level::Alarm) << "EMBEDDED" << std::endl;
					// OGDF_ASSERT(G.representsCombEmbedding());
				} else {
					Logger::slout(Logger::Level::Alarm) << "SOLVING FAILED" << std::endl;
				}
			} else {
				Logger::slout(Logger::Level::Alarm) << "REDUCTION FAILED" << std::endl;
			}
		} else {
			Logger::slout(Logger::Level::Alarm) << "NON-PLANAR" << std::endl;
		}

		if (outfile != nullptr && !GraphIO::write(GA, outfile)) {
			Logger::slout(Logger::Level::Alarm) << "IO Error: Couldn't write graph!" << std::endl;
			return ERROR_IO;
		}

		return solved ? PQPLANAR : NOT_PQPLANAR;
	} catch (invalid_option& e) {
		Logger::slout(Logger::Level::Alarm) << "Invalid Options: " << e.what() << std::endl;
		return ERROR_OPTIONS;
	}
}
