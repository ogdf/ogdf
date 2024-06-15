/** \file
 * \brief TODO Document
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
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/cluster/CconnectClusterPlanar.h>
#include <ogdf/cluster/CconnectClusterPlanarEmbed.h>
#include <ogdf/cluster/ClusterPlanarity.h>
#include <ogdf/cluster/HananiTutteCPlanarity.h>
#include <ogdf/external/coin.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/planarlayout/FPPLayout.h>

#include <json.hpp>
#include <ostream>

#include "PQPlanarity.h"
#include "PipeOrder.h"
#include "likwid_utils.h"
#include "return.h"
#include "utils/Logging.h"

#ifdef OGDF_DEBUG
#	define OGDF_DEBUG_PARAM(x) x,
#else
#	define OGDF_DEBUG_PARAM(x)
#endif

using nlohmann::json;

static std::string progname;

void usage() {
	Logger::slout(Logger::Level::Alarm)
			<< "Usage: " << progname << " input\n\n"
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

json solvePQPlan(PQPlanarity& pq, const PQPlanConf& conf);

int exit_code = ERROR_ASSERT;

int main(int argc, char* argv[]) {
	CoinManager::CoinLog.localLogMode(Logger::LogMode::Log);
	CoinManager::CoinLog.localLogLevel(Logger::Level::Alarm);
	progname = argv[0];
	Logger::Level log_level = Logger::Level::Minor;
	PQPlanConf pqconf;

	int opt;
	while ((opt = getopt(argc, argv, "l:dcirbasp")) != -1) {
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
	std::string infile = argv[optind++];
	if (optind < argc) {
		Logger::slout(Logger::Level::Alarm)
				<< "Invalid Options: Too many positional arguments!" << std::endl;
		usage();
		return ERROR_OPTIONS;
	}

	Graph G;
#ifdef OGDF_DEBUG
	GraphAttributes GA(G, GraphAttributes::all);
#else
	GraphAttributes GA(G,
			GraphAttributes::nodeLabel
					| GraphAttributes::edgeLabel); // only for getting options, not passed to PQ
#endif
	int read_ret = readCG(infile, log_level, GA, G);
	if (read_ret != SUCCESS) {
		return read_ret;
	}
	if (!isPlanar(G)) {
		Logger::slout(Logger::Level::Alarm)
				<< "Assertion Error: Input graph is not planar!" << std::endl;
		return ERROR_ASSERT;
	}

#ifdef OGDF_DEBUG
	PQPlanarity pq(&G, &GA);
#else
	PQPlanarity pq(&G);
#endif
	try {
		std::ifstream i(infile + ".json");
		if (!i.good()) {
			Logger::slout(Logger::Level::Alarm)
					<< "IO Error: Couldn't read config json " << infile << ".json!" << std::endl;
			return ERROR_IO;
		}
		json j;
		i >> j;
		PQPlanOptions::applyConfigJSON(G, GA, pq, j);
	} catch (json::parse_error& e) {
		Logger::slout(Logger::Level::Alarm) << "IO Error: Couldn't parse config json " << infile
											<< ".json: " << e.what() << std::endl
											<< "exception id: " << e.id << std::endl
											<< "byte position of error: " << e.byte << std::endl;
		return ERROR_IO;
	}

	Logger::slout() << pq << std::endl;
	json result = solvePQPlan(pq, pqconf);
	std::cerr << result << std::endl;
	return exit_code;
}

json solvePQPlan(PQPlanarity& pq, const PQPlanConf& conf) {
	Logger::slout() << "PQPlanConf: " << conf << std::endl;

	json ret;
	ret["method"] = "PQPlanarity";
	ret["mode"] = "PQPlanarity" + conf.getID();
	ret["config"] = conf;

	likwid_prepare(ret);

	int64_t dur_init_stats, dur_reduced_stats;
	tp start = tpc::now();
	pq.matchings.setPipeQueue(conf.getOrder(&pq));
	pq.setAllowContractBBPipe(conf.allow_contract);
	pq.setIntersectTrees(conf.intersect_trees);
	pq.setBatchSpqr(conf.batch_spqr);
	pq.stats_out.open("pqplan_stats.json");
	if (!pq.stats_out.is_open() || !pq.stats_out.good()) {
		Logger::slout(Logger::Level::Alarm)
				<< "IO Warning: Could not open pqplan_stats.json for writing!" << std::endl;
	} else {
		pq.stats_out << "[";
	}
	{
		tp t = tpc::now();
		json s_init;
		pqPlanStats(pq, s_init);
		ret["init_stats"] = s_init;
		dur_init_stats = dur_ns(tpc::now() - t);
	}
	tp init = tpc::now();
	ret["time_init_ns"] = dur_ns(init - start) - dur_init_stats;

	bool reduced = pq.makeReduced();
	{
		tp t = tpc::now();
		json s_reduced;
		pqPlanStats(pq, s_reduced);
		ret["reduced_stats"] = s_reduced;
		ret["undo_ops"] = pq.undoOperations();
		dur_reduced_stats = dur_ns(tpc::now() - t);
	}
	tp make_reduced = tpc::now();
	ret["time_make_reduced_ns"] = dur_ns(make_reduced - init) - dur_reduced_stats;

	exit_code = NOT_PQPLANAR;
	bool g_comb_emb = false;
	if (reduced) {
		ret["status"] = "preSolve";
		std::cerr << ret << std::endl;

		reduced = pq.solveReduced();
		tp solve_reduced = tpc::now();
		ret["time_solve_reduced_ns"] = dur_ns(solve_reduced - make_reduced);
		ret["time_ns"] = dur_ns(solve_reduced - start) - dur_init_stats - dur_reduced_stats;
		if (reduced) {
			ret["status"] = "preEmbed";
			std::cerr << ret << std::endl;

			pq.embed();
			tp embed = tpc::now();
			g_comb_emb = ret["g_comb_emb"] = pq.G->representsCombEmbedding();
			tp check = tpc::now();

			exit_code = PQPLANAR;
			ret["result"] = true;
			ret["status"] = "embeddedAndVerified";

			ret["time_embed_ns"] = dur_ns(embed - solve_reduced);
			ret["time_check_ns"] = dur_ns(check - embed);
		} else {
			ret["result"] = false;
			ret["status"] = "solvingFailed";
		}
	} else {
		ret["result"] = false;
		ret["status"] = "reductionFailed";
		ret["time_ns"] = dur_ns(make_reduced - start) - dur_init_stats - dur_reduced_stats;
	}
	if (exit_code == PQPLANAR && !g_comb_emb) {
		exit_code = ERROR_COMB_EMB;
	}
	if (!pq.stats_out.is_open() || !pq.stats_out.good()) {
		Logger::slout(Logger::Level::Alarm)
				<< "IO Warning: Could not finish writing to pqplan_stats.json!" << std::endl;
	} else {
		pq.stats_out << "]" << std::endl;
	}
	likwid_finalize(ret);
	return ret;
}
