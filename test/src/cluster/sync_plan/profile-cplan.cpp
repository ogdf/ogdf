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

#include "SyncPlan.h"
#include "PipeOrder.h"
#include "likwid_utils.h"
#include "return.h"
#include "utils/Clusters.h"

#ifdef OGDF_DEBUG
#	define OGDF_DEBUG_PARAM(x) x,
#else
#	define OGDF_DEBUG_PARAM(x)
#endif


static std::string progname;

void usage() {
	Logger::slout(Logger::Level::Alarm)
			<< "Usage: " << progname << " input\n\n"
			<< "-l[0-5]  set log level\n"
			<< "-m METHOD = {SyncPlan, CConnected, HananiTutte, ILP}\n"
			<< "-t HH:MM:SS  timeout for ILP\n\n"
			<< "-f run Hanani Tutte without verifying\n\n"
			<< "Options for SyncPlan:\n"
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

enum class Method {
	Invalid,
	SyncPlan,
	CConnected,
	HananiTutte,
	ILP,
};

NLOHMANN_JSON_SERIALIZE_ENUM(Method,
		{
				{Method::Invalid, "Invalid"},
				{Method::SyncPlan, "SyncPlan"},
				{Method::CConnected, "CConnected"},
				{Method::HananiTutte, "HananiTutte"},
				{Method::ILP, "ILP"},
		})

json solveSyncPlan(OGDF_DEBUG_PARAM(ClusterGraphAttributes& GA) ClusterGraph& CG, Graph& G,
		const SyncPlanConf& conf);

json solveCConnected(ClusterGraph& CG, Graph& G);

json solveHananiTutte(ClusterGraph& CG, Graph& G, bool fast);

json solveILP(ClusterGraph& CG, Graph& G, string& timeout);

int exit_code = ERROR_ASSERT;

int main(int argc, char* argv[]) {
	CoinManager::CoinLog.localLogMode(Logger::LogMode::Log);
	CoinManager::CoinLog.localLogLevel(Logger::Level::Alarm);
	progname = argv[0];
	Logger::Level log_level = Logger::Level::Minor;
	Method method;
	string timeout;
	bool ht_fast = false;
	SyncPlanConf pqconf;

	int opt;
	while ((opt = getopt(argc, argv, "l:dm:t:fcirbasp")) != -1) {
		switch (opt) {
		case 'l':
			log_level = parseLogLevel();
			break;
#ifdef OGDF_DEBUG
		case 'd':
			SyncPlanConsistency::doWriteOut = true;
			break;
#endif
		case 'm': {
			json in = json::parse(string("\"") + optarg + "\"");
			method = in.get<Method>();
			break;
		}
		case 't':
			timeout = optarg;
			break;
		case 'f':
			ht_fast = true;
			break;
#define TOGGLE(key, var)          \
	case key:                     \
		pqconf.var = !pqconf.var; \
		break;
			SyncPlanConf_KEYS
#undef TOGGLE
					default : /* '?' */
							  usage();
			return ERROR_OPTIONS;
		}
	}
	if (method == Method::Invalid) {
		json methods = json::array(
				{Method::SyncPlan, Method::CConnected, Method::HananiTutte, Method::ILP});
		Logger::slout(Logger::Level::Alarm)
				<< "Invalid Options: Please specify a method from " << methods << "!" << std::endl;
		usage();
		return ERROR_OPTIONS;
	}
	if (method == Method::ILP and timeout.empty()) {
		Logger::slout(Logger::Level::Alarm)
				<< "Invalid Options: Missing timeout parameter!" << std::endl;
		usage();
		return ERROR_OPTIONS;
	}
	if (optind >= argc) {
		Logger::slout(Logger::Level::Alarm) << "Invalid Options: Expected infile!" << std::endl;
		usage();
		return ERROR_OPTIONS;
	}
	char* infile = argv[optind++];
	if (optind < argc) {
		Logger::slout(Logger::Level::Alarm)
				<< "Invalid Options: Too many positional arguments!" << std::endl;
		usage();
		return ERROR_OPTIONS;
	}

	Graph G;
	ClusterGraph CG(G);
	ClusterGraphAttributes GA(CG, ClusterGraphAttributes::all);
	int read_ret = readCG(infile, log_level, GA, CG, G);
	if (read_ret != SUCCESS) {
		return read_ret;
	}
	if (!isPlanar(G)) {
		Logger::slout(Logger::Level::Alarm)
				<< "Assertion Error: Input graph is not planar!" << std::endl;
		return ERROR_ASSERT;
	}
	if (!isParallelFree(G)) {
		Logger::slout(Logger::Level::Alarm)
				<< "Assertion Error: Input graph is not parallel free!" << std::endl;
		return ERROR_ASSERT;
	}
	if (!isLoopFree(G)) {
		Logger::slout(Logger::Level::Alarm)
				<< "Assertion Error: Input graph is not loop free!" << std::endl;
		return ERROR_ASSERT;
	}

	for (edge e : G.edges) {
		if (e->source()->index() > e->target()->index()) {
			G.reverseEdge(e);
		}
	}

	json result;
	switch (method) {
	case Method::SyncPlan:
		result = solveSyncPlan(OGDF_DEBUG_PARAM(GA) CG, G, pqconf);
		break;
	case Method::CConnected:
		result = solveCConnected(CG, G);
		break;
	case Method::HananiTutte:
		result = solveHananiTutte(CG, G, ht_fast);
		break;
	case Method::ILP:
		result = solveILP(CG, G, timeout);
		break;
	default:
		Logger::slout(Logger::Level::Alarm) << "Invalid Options: Invalid method!" << std::endl;
		return ERROR_OPTIONS;
	}

	std::cerr << result << std::endl;
	return exit_code;
}

template<typename EC>
std::string CconnectClusterPlanarErrorCode2Str(EC code, bool cplan) {
	switch (code) {
	case EC::none:
		return cplan ? "CConnectedCPlanar" : "unknownError";
	case EC::nonConnected:
		return "nonConnected";
	case EC::nonCConnected:
		return "nonCConnected";
	case EC::nonPlanar:
		return "nonPlanar";
	case EC::nonCPlanar:
		return "nonCPlanar";
	default:
		return "???";
	}
}

json solveSyncPlan(OGDF_DEBUG_PARAM(ClusterGraphAttributes& GA) ClusterGraph& CG, Graph& G,
		const SyncPlanConf& conf) {
	Logger::slout() << "SyncPlanConf: " << conf << std::endl;

	json ret;
	ret["method"] = Method::SyncPlan;
	ret["mode"] = ret["method"].get<string>() + conf.getID();
	ret["config"] = conf;

	likwid_prepare(ret);

	int64_t dur_init_stats, dur_reduced_stats;
	tp start = tpc::now();
#ifdef OGDF_DEBUG
	SyncPlan pq(&G, &CG, &GA);
#else
	SyncPlan pq(&G, &CG);
#endif
	pq.matchings.setPipeQueue(conf.getOrder(&pq));
	pq.setAllowContractBBPipe(conf.allow_contract);
	pq.setIntersectTrees(conf.intersect_trees);
	pq.setBatchSpqr(conf.batch_spqr);
	pq.stats_out.open("syncplan_stats.json");
	if (!pq.stats_out.is_open() || !pq.stats_out.good()) {
		Logger::slout(Logger::Level::Alarm)
				<< "IO Warning: Could not open syncplan_stats.json for writing!" << std::endl;
	} else {
		pq.stats_out << "[";
	}
	{
		tp t = tpc::now();
		json s_init;
		syncPlanStats(pq, s_init);
		ret["init_stats"] = s_init;
		dur_init_stats = dur_ns(tpc::now() - t);
	}
	tp init = tpc::now();
	ret["time_init_ns"] = dur_ns(init - start) - dur_init_stats;

	bool reduced = pq.makeReduced();
	{
		tp t = tpc::now();
		json s_reduced;
		syncPlanStats(pq, s_reduced);
		ret["reduced_stats"] = s_reduced;
		ret["undo_ops"] = pq.undoOperations();
		dur_reduced_stats = dur_ns(tpc::now() - t);
	}
	tp make_reduced = tpc::now();
	ret["time_make_reduced_ns"] = dur_ns(make_reduced - init) - dur_reduced_stats;

	exit_code = NOT_SYNC_PLAN;
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
			ret["g_comb_emb"] = G.representsCombEmbedding();
			ret["cg_comb_emb"] = isClusterPlanarEmbedding(CG);
			tp check = tpc::now();

			exit_code = SYNC_PLAN;
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
	bool g_comb_emb = ret["g_comb_emb"] = G.representsCombEmbedding();
	bool cg_comb_emb = ret["cg_comb_emb"] = isClusterPlanarEmbedding(CG);
	if (exit_code == SYNC_PLAN && (!g_comb_emb || !cg_comb_emb)) {
		exit_code = ERROR_COMB_EMB;
	}
	if (!pq.stats_out.is_open() || !pq.stats_out.good()) {
		Logger::slout(Logger::Level::Alarm)
				<< "IO Warning: Could not finish writing to syncplan_stats.json!" << std::endl;
	} else {
		pq.stats_out << "]" << std::endl;
	}
	likwid_finalize(ret);
	return ret;
}

json solveCConnected(ClusterGraph& CG, Graph& G) {
	json ret;
	ret["mode"] = ret["method"] = Method::CConnected;

	tp t1 = tpc::now();
	bool ccon = isCConnected(CG);
	tp t2 = tpc::now();
	ret["c_connected"] = ccon;
	ret["time_precheck_ns"] = dur_ns(t2 - t1);
	ret["status"] = "preTest";
	std::cerr << ret << std::endl;

	CconnectClusterPlanar ccPlanarityTest;
	bool cplan = ccPlanarityTest.call(CG);
	tp t3 = tpc::now();
	exit_code = cplan ? SYNC_PLAN : NOT_SYNC_PLAN;
	ret["result"] = cplan;
	ret["status"] = CconnectClusterPlanarErrorCode2Str(ccPlanarityTest.errCode(), cplan);
	ret["time_ns"] = dur_ns(t3 - t2);
	ret["status_embed"] = "preEmbed";
	std::cerr << ret << std::endl;

	CconnectClusterPlanarEmbed ccEmbedder;
	bool cemb = ccEmbedder.embed(CG, G);
	tp t4 = tpc::now();
	ret["result_embed"] = cemb;
	ret["status_embed"] = CconnectClusterPlanarErrorCode2Str(ccEmbedder.errCode(), cemb);
	ret["time_embed_ns"] = dur_ns(t4 - t3);

	bool g_comb_emb = ret["g_comb_emb"] = G.representsCombEmbedding();
	bool cg_comb_emb = ret["cg_comb_emb"] = isClusterPlanarEmbedding(CG);
	if (exit_code == SYNC_PLAN && (!g_comb_emb || !cg_comb_emb)) {
		exit_code = ERROR_COMB_EMB;
	}

	return ret;
}

json solveHananiTutte(ClusterGraph& CG, Graph& G, bool fast) {
	std::unique_ptr<HananiTutteCPlanarity::HananiTutteSolver> solver {
			HananiTutteCPlanarity::getSolver(CG)};
	HananiTutteCPlanarity::Stats stats;
	json ret;
	ret["mode"] = ret["method"] = Method::HananiTutte;
	if (fast) {
		ret["mode"] = "HananiTutte-f";
	}

	tp start = tpc::now();
	bool test = solver->test(stats);
	tp split = tpc::now();
	ret["test_result"] = test;
	ret["time_test_ns"] = dur_ns(split - start);
	ret["status"] = "preVerify";
	std::cerr << ret << std::endl;

	if (test) {
		bool verify = fast || solver->verify(stats);
		tp stop = tpc::now();
		ret["verify_result"] = ret["result"] = verify;
		exit_code = SYNC_PLAN;
		if (fast) {
			ret["status"] = "maybeCPlanar/notVerified";
			ret["verify_result"] = false;
		} else if (verify) {
			ret["status"] = "CPlanarVerified"; // TODO also copy generated embedding to ClusterGraph
		} else {
			ret["status"] = "maybeCPlanar/VerificationFailed";
			exit_code = ERROR_COMB_EMB;
		}
		ret["time_verify_ns"] = dur_ns(stop - split);
		ret["time_ns"] = dur_ns(stop - start);
	} else {
		exit_code = NOT_SYNC_PLAN;
		ret["result"] = false;
		ret["status"] = "nonCPlanarVerified";
		ret["time_ns"] = dur_ns(split - start);
	}

	json s;
	s["nRows"] = stats.nRows;
	s["nColumns"] = stats.nColumns;
	s["nConditions"] = stats.nConditions;
	s["nMoves"] = stats.nMoves;
	s["tPrepare"] = stats.tPrepare;
	s["tCreateSparse"] = stats.tCreateSparse;
	s["tSolve"] = stats.tSolve;
	s["tCheck"] = stats.tCheck;
	ret["statistics"] = s;

	bool g_comb_emb = ret["g_comb_emb"] = G.representsCombEmbedding();
	bool cg_comb_emb = ret["cg_comb_emb"] = isClusterPlanarEmbedding(CG);
	if (exit_code == SYNC_PLAN && (!g_comb_emb || !cg_comb_emb)) {
		// exit_code = ERROR_COMB_EMB; // FIXME Hanani Tutte does embed!
	}

	return ret;
}

json solveILP(ClusterGraph& CG, Graph& G, string& timeout) {
	ClusterPlanarity cPlanarity;
	json ret;
	ret["mode"] = ret["method"] = Method::ILP;
	cPlanarity.setTimeLimit(timeout);
	tp start = tpc::now();
	ClusterPlanarity::NodePairs addedEdges;
	bool result = cPlanarity.isClusterPlanar(CG, addedEdges);
	ret["result"] = result;
	ret["time_ns"] = dur_ns(tpc::now() - start);
	ret["status_val"] = cPlanarity.getOptStatus();
	if (cPlanarity.getOptStatus() == abacus::Master::Optimal) {
		ret["status"] = result ? "Optimal-cplan" : "Optimal-not-cplan";
		exit_code = result ? SYNC_PLAN : NOT_SYNC_PLAN;
	} else {
		ret["status"] = string(abacus::Master::STATUS_[cPlanarity.getOptStatus()]);
		exit_code = ERROR_ABACUS + cPlanarity.getOptStatus();
	}

	json s;
	s["AddedEdges"] = addedEdges.size();
	s["TotalTime"] = cPlanarity.getTotalTime();
	s["HeurTime"] = cPlanarity.getHeurTime();
	s["LPTime"] = cPlanarity.getLPTime();
	s["LPSolverTime"] = cPlanarity.getLPSolverTime();
	s["SeparationTime"] = cPlanarity.getSeparationTime();
	s["TotalWTime"] = cPlanarity.getTotalWTime();
	s["NumCCons"] = cPlanarity.getNumCCons();
	s["NumKCons"] = cPlanarity.getNumKCons();
	s["NumLPs"] = cPlanarity.getNumLPs();
	s["NumBCs"] = cPlanarity.getNumBCs();
	s["NumSubSelected"] = cPlanarity.getNumSubSelected();
	s["NumVars"] = cPlanarity.getNumVars();
	ret["statistics"] = s;

	bool g_comb_emb = ret["g_comb_emb"] = G.representsCombEmbedding();
	bool cg_comb_emb = ret["cg_comb_emb"] = isClusterPlanarEmbedding(CG);
	if (exit_code == SYNC_PLAN && (!g_comb_emb || !cg_comb_emb)) {
		// exit_code = ERROR_COMB_EMB; // ILP does not embed // FIXME?
	}

	return ret;
}
