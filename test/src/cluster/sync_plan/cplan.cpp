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
#include <ogdf/cluster/ClusterOrthoLayout.h>
#include <ogdf/cluster/ClusterPlanarizationLayout.h>
#include <ogdf/external/coin.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/planarlayout/FPPLayout.h>

#include "PQPlanarity.h"
#include "PipeOrder.h"
#include "return.h"
#include "utils/Clusters.h"
#include "utils/Logging.h"

static std::string progname;

void usage() {
	Logger::slout(Logger::Level::Alarm) << "Usage: " << progname << " input [outfile]" << std::endl;
}

int main(int argc, char* argv[]) {
	CoinManager::CoinLog.localLogMode(Logger::LogMode::Log);
	CoinManager::CoinLog.localLogLevel(Logger::Level::Alarm);
	progname = argv[0];
	char* outfile = nullptr;
	char* svgfile = nullptr;
	Logger::Level log_level = Logger::Level::Minor;
	int layout_mode = 1;

	int opt;
	while ((opt = getopt(argc, argv, "l:da:")) != -1) {
		switch (opt) {
		case 'l':
			log_level = parseLogLevel();
			break;
#ifdef OGDF_DEBUG
		case 'd':
			PQPlanarityConsistency::doWriteOut = true;
			break;
#endif
		case 'a':
			try {
				layout_mode = std::stoi(optarg);
				if (1 > layout_mode || layout_mode > 3) {
					throw invalid_option("out of bounds");
				}
			} catch (const std::exception& e) {
				throw invalid_option(string_format(
						"Layout algorithm must be a number between 1 and 3! (%s)", e.what()));
			}
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
	char* infile = argv[optind++];
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

	Graph G;
	ClusterGraph CG(G);
	ClusterGraphAttributes GA(CG, ClusterGraphAttributes::all);
	int read_ret = readCG(infile, log_level, GA, CG, G);
	if (read_ret != SUCCESS) {
		return read_ret;
	}
	printClusters(CG.rootCluster(), Logger::slout(Logger::Level::High) << CG << ": ") << std::endl;
#ifdef OGDF_DEBUG
	bool was_planar = isPlanar(G);
#endif
	PQPlanarity pq(&G, &CG, &GA);
	pq.setAllowContractBBPipe(false);
	pq.setIntersectTrees(false);
	pq.setBatchSpqr(true);
	pq.matchings.setPipeQueue(new PipeQueueByDegree());

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
				printClusters(CG.rootCluster(), Logger::slout(Logger::Level::High) << CG << ": ")
						<< std::endl;
				OGDF_ASSERT(G.representsCombEmbedding());
				// OGDF_ASSERT(CG.representsCombEmbedding());
				OGDF_ASSERT(isClusterPlanarEmbedding(CG));
			} else {
				Logger::slout(Logger::Level::Alarm) << "SOLVING FAILED" << std::endl;
			}
		} else {
			Logger::slout(Logger::Level::Alarm) << "REDUCTION FAILED" << std::endl;
		}
	} else {
		OGDF_ASSERT(!was_planar);
		Logger::slout(Logger::Level::Alarm) << "NON-PLANAR" << std::endl;
	}

	if (outfile != nullptr && !GraphIO::write(GA, outfile)) {
		Logger::slout(Logger::Level::Alarm) << "IO Error: Couldn't write graph!" << std::endl;
		return ERROR_IO;
	}
	if (solved && svgfile != nullptr) {
		adjEntry adjExternal = CG.rootCluster()->children.front()->adjEntries.front();
		for (cluster c : CG.clusters) {
			GA.strokeColor(c) = colors[c->index() % colors.size()];
			std::stringstream ss;
			ss << "cluster " << c->index();
			GA.label(c) = ss.str();
		}

		if (layout_mode == 1) {
			const auto pair = drawClusterGraph(CG, GA);
			if (!GraphIO::write(pair->second, svgfile)) {
				Logger::slout(Logger::Level::Alarm) << "IO Error: Couldn't draw graph!" << std::endl;
				return ERROR_IO;
			}
		} else if (layout_mode == 2) {
			ClusterOrthoLayout col;
			ClusterPlanRep rep(GA, CG);
			Layout drawing(G);
			col.call(rep, adjExternal, drawing);
			for (node n : G.nodes) {
				GA.x(n) = drawing.x(n);
				GA.y(n) = drawing.y(n);
			}
			for (edge e : G.edges) {
				GA.bends(e) = drawing.bends(e);
			}

			if (!GraphIO::write(GA, svgfile)) {
				Logger::slout(Logger::Level::Alarm) << "IO Error: Couldn't draw graph!" << std::endl;
				return ERROR_IO;
			}
		} else if (layout_mode == 3) {
			ClusterPlanarizationLayout cpl;
			cpl.call(G, GA, CG, true);

			GA.scale(0.25);
			for (node n : G.nodes) {
				GA.label(n) = "";
				GA.strokeColor(n) = ogdf::Color();
				GA.fillColor(n) = ogdf::Color();
				GA.width(n) = GA.height(n) = 3;
			}
			for (edge e : G.edges) {
				GA.strokeWidth(e) = 1;
			}
			for (cluster c : CG.clusters) {
				GA.label(c) = "";
				GA.strokeWidth(c) = 2;
			}

			if (!GraphIO::write(GA, svgfile)) {
				Logger::slout(Logger::Level::Alarm) << "IO Error: Couldn't draw graph!" << std::endl;
				return ERROR_IO;
			}
		}
	}

	return solved ? PQPLANAR : NOT_PQPLANAR;
}
