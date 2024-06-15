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
#include <ogdf/basic/graph_generators.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/fileformats/GraphIO.h>

#include <filesystem>

#include <getopt.h>

#include "PQPlanarityAttributes.h"
#include "return.h"
#include "utils/Clusters.h"
#include "utils/Preprocess.h"

using namespace ogdf;
static std::string progname;

void usage() {
	Logger::slout(Logger::Level::Alarm)
			<< "Usage: " << progname << " [-options] input [outfile] [svgfile]\n\n"
			<< "-f generate even if outfile already exists\n"
			<< "-a preprocess clusters before writing outfile (but breaks embedding information and disables svg writing)\n"
			<< "-l[0-5] set log level\n"
			<< "-s[int] set seed for random number generator or 'r' to use a true random seed\n"

			<< "-p[double] set probability for not adding anything further nodes to the cluster\n"
			<< "-n[double] set the expected number of nodes per cluster\n"
			<< "-m[int] set maximum number of nodes per cluster\n"

			<< "-u[double] set probability for not adding a further cluster\n"
			<< "-c[double] set the expected number of clusters, based on the number of nodes in the input graph\n"
			<< "-d[int] set maximum number of clusters\n"

			<< "-r[int] set minimum number of nodes in the root cluster\n"
			<< "-o only generate connected clusters\n"
			<< "-t[int] timeout in seconds\n\n"

			<< "outfile and svgfile can use the following templates that will be replaced based on infile and the given options:\n"
			<< "%s | tag (seed + prob) | -s123-p0.01\n"
			<< "%d | infile directory  | clusters/\n"
			<< "%n | infile basename   | graph01\n"
			<< "%e | infile extension  | .gml" << std::endl;
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
	RandomClusterConfig config;
	double expected_clusters = -1;
	int seed = 0, mode = 0;
	bool force = false, preprocess = false;

	int opt;
	while ((opt = getopt(argc, argv, "hfal:j:s:p:n:m:u:c:d:r:ot:")) != -1) {
		switch (opt) {
		case 'f':
			force = true;
			break;
		case 'a':
			preprocess = true;
			break;
		case 'l':
			log_level = parseLogLevel();
			break;
		case 'j':
			mode = std::stoi(optarg);
			break;
		case 's':
			if (string("r") == optarg) {
				std::random_device rd;
				seed = rd();
			} else {
				seed = std::stoi(optarg);
			}
			break;
		case 'p':
			config.prob_no_further_node = std::stod(optarg);
			break;
		case 'n':
			config.expected_nodes(std::stoi(optarg));
			break;
		case 'm':
			config.max_nodes_in_cluster = std::stoi(optarg);
			break;
		case 'u':
			config.prob_no_further_cluster = std::stod(optarg);
			break;
		case 'c':
			expected_clusters = std::stod(optarg);
			break;
		case 'd':
			config.max_clusters = std::stoi(optarg);
			break;
		case 'r':
			config.min_root_nodes = std::stoi(optarg);
			break;
		case 'o':
			config.cconnected = true;
			break;
		case 't':
			config.timeout = std::stoi(optarg);
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

	Graph G;
	ClusterGraph CG(G);
	CG.setUpdateDepth(true);
	int read_ret = readCG(infile, log_level, CG, G);
	if (read_ret != SUCCESS) {
		return read_ret;
	}
	if (expected_clusters > 0) {
		config.expected_nodes((double)G.numberOfNodes() / expected_clusters);
	}

	const std::filesystem::path inpath = infile;
	std::stringstream ss;
	if (mode == 0) {
		ss << "-s" << seed << "-p" << config.prob_no_further_node;
	} else {
		ss << "-s" << seed << "-d" << config.max_clusters;
	}
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
	Logger::slout(Logger::Level::High) << "Config is: " << config << std::endl;
	if (!planarEmbed(G)) {
		Logger::slout(Logger::Level::Alarm) << "Error: The given graph is not planar!" << std::endl;
		return NOT_PQPLANAR;
	}
#ifdef OGDF_DEBUG
	CG.consistencyCheck();
#endif
	if (mode == 0) {
		Logger::slout(Logger::Level::High) << "Calling new (planar) Clusterer" << std::endl;
		if (!makeClusters(CG, config)) {
			Logger::slout(Logger::Level::High) << "Clustering timed out!" << std::endl;
		}
	} else {
		List<edge> added;
		makeConnected(G, added);
		if (mode == 1) {
			Logger::slout(Logger::Level::High) << "Calling old planar Clusterer" << std::endl;
			randomClusterPlanarGraph(CG, G, config.max_clusters);
		} else {
			Logger::slout(Logger::Level::High) << "Calling old non-planar Clusterer" << std::endl;
			randomClusterGraph(CG, G, config.max_clusters);
		}
		for (edge e : added) {
			G.delEdge(e);
		}
	}
#ifdef OGDF_DEBUG
	CG.consistencyCheck();
#endif
	printCG(CG, "Output ");

	bool ignore = false;
	if (preprocess && preprocessClusterGraph(CG, G)) {
		printCG(CG, "Preprocessed ");
		if (!isPlanar(G)) {
			ignore = true;
			Logger::slout(Logger::Level::High)
					<< "Graph became non-planar through preprocessing!" << std::endl;
		} else if (CG.numberOfClusters() < 2) {
			ignore = true;
			Logger::slout(Logger::Level::High)
					<< "ClusterGraph became trivial through preprocessing!" << std::endl;
		} else if (isCConnected(CG)) {
			ignore = true;
			Logger::slout(Logger::Level::High) << "ClusterGraph is C-connected!" << std::endl;
		}
		CG.adjAvailable(false);
	}

	Logger::slout(Logger::Level::High) << "Writing to " << outfile << std::endl;
	int write_res = writeCG(outfile, CG);
	if (write_res != PQPLANAR) {
		return write_res;
	}
	if (ignore) {
		Logger::slout(Logger::Level::High) << "Renaming to " << outfile << ".ignore" << std::endl;
		rename((std::filesystem::path)outfile, (std::filesystem::path)(outfile + ".ignore"));
	}

	if (CG.adjAvailable()) {
		if (!svgfile.empty()) {
			format(svgfile, "%s", tag);
			format(svgfile, "%d", inpath.parent_path());
			format(svgfile, "%n", inpath.stem());
			format(svgfile, "%e", inpath.extension());

			ClusterGraphAttributes CGA(CG, GraphAttributes::all);
			const auto pair = drawClusterGraph(CG, CGA);
			Logger::slout(Logger::Level::High) << "Drawing to " << svgfile << std::endl;
			if (!GraphIO::write(pair->second, svgfile)) {
				Logger::slout(Logger::Level::Alarm) << "IO Error: Couldn't draw graph!" << std::endl;
				return ERROR_IO;
			}
		}

		if (!isClusterPlanarEmbedding(CG)) {
			Logger::slout(Logger::Level::Alarm)
					<< "Assertion Error: Generated embedding is not cluster planar!" << std::endl;
			return NOT_PQPLANAR;
		}
	}

	return PQPLANAR;
}
