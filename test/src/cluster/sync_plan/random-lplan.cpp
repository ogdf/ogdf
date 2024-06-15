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
#include <getopt.h>

#include "return.h"
#include "utils/Levels.h"

static std::string progname;

void usage() {
	Logger::slout(Logger::Level::Alarm)
			<< "Usage: " << progname << " [-options] [outfile] [clusterfile] [svgfile]\n\n"
			<< "-n[int] set number of nodes\n"
			<< "-m[int] set number of edges\n"
			<< "-l[int] set number of levels\n"
			<< "-d[int] set minimum degree\n"
			<< "-s[int] set seed for random number generator or 'r' to use a true random seed\n"
			<< "-r create a radial instance\n"
			<< "-e add extra edge spanning from highest to lowest level to cluster instance\n"
			<< std::endl;
}

int main(int argc, char* argv[]) {
	progname = argv[0];
	std::string outfile = "-", clusterfile, svgfile;
	int seed = 0, nodes = 0, edges = 0, levels = 0, min_deg = 0;
	bool radial = false, extra_edge = false;

	int opt;
	while ((opt = getopt(argc, argv, "ren:m:l:s:d")) != -1) {
		switch (opt) {
		case 'r':
			radial = true;
			break;
		case 'e':
			extra_edge = true;
			break;
		case 's':
			if (string("r") == optarg) {
				std::random_device rd;
				seed = rd();
			} else {
				seed = std::stoi(optarg);
			}
			break;
		case 'n':
			nodes = std::stoi(optarg);
			break;
		case 'm':
			edges = std::stoi(optarg);
			break;
		case 'l':
			levels = std::stoi(optarg);
			break;
		case 'd':
			min_deg = std::stoi(optarg);
			break;
		default: /* '?' */
			usage();
			return ERROR_OPTIONS;
		}
	}

	if (nodes < 1) {
		Logger::slout(Logger::Level::Alarm)
				<< "Invalid Options: Number of nodes is required!" << std::endl;
		usage();
		return ERROR_OPTIONS;
	}
	if (levels < 1) {
		Logger::slout(Logger::Level::Alarm)
				<< "Invalid Options: Number of levels is required!" << std::endl;
		usage();
		return ERROR_OPTIONS;
	}

	if (optind < argc) {
		outfile = argv[optind++];
	}
	if (optind < argc) {
		clusterfile = argv[optind++];
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

	setSeed(seed);

	Graph G;
	NodeArray<int> lvl(G, -1);
	vector<vector<node>> emb;
	NodeArray<int> pos(G, -1);

	randomProperMaximalLevelPlaneGraph(G, emb, nodes, levels, radial);
	embedPLE(G, emb, lvl, pos);
	if (edges > 0) {
		pruneEdges(G, edges, min_deg);
	}

	if (outfile == "-") {
		writeLevelGraph(G, emb, pos, std::cout);
	} else {
		std::ofstream fout(outfile);
		writeLevelGraph(G, emb, pos, fout);
	}

	if (!svgfile.empty()) {
		drawSVG(G, lvl, pos, svgfile);
	}

	if (!clusterfile.empty()) {
		Graph Gc;
		ClusterGraph CG(Gc);
		reduceLevelToCluster(G, emb, Gc, CG);
		if (extra_edge) {
			node u = nullptr;
			for (cluster c : CG.clusters) {
				node v = Gc.newNode();
				CG.reassignNode(v, c);
				if (u) {
					Gc.newEdge(v, u);
				}
				u = v;
			}
		}
		GraphIO::write(CG, clusterfile);
	}

	// std::ifstream fin("filename.txt");
	// readLevelGraph(G, emb, pos, fin);
	// embedPLE(G, emb, lvl, pos);

	//checkPLE(G, emb, lvl, pos);

	return 0;
}
