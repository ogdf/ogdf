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
#include <ogdf/fileformats/GraphIO.h>

#include <getopt.h>

#include "return.h"
#include "utils/Clusters.h"
#include "utils/Preprocess.h"

using namespace std;
using namespace ogdf;
static std::string progname;

void usage() {
	Logger::slout(Logger::Level::Alarm) << "Usage: " << progname << " [-options] input [outfile]\n\n"
										<< "-l[0-5] set log level\n"
										<< std::endl;
}

int main(int argc, char* argv[]) {
	progname = argv[0];
	std::string outfile;
	Logger::Level log_level = Logger::Level::Minor;

	int opt;
	while ((opt = getopt(argc, argv, "hl:")) != -1) {
		switch (opt) {
		case 'l':
			log_level = parseLogLevel();
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

	printCG(CG, "Input ");
	if (preprocessClusterGraph(CG, G)) {
		CG.adjAvailable(false);
	}
	printCG(CG, "Output ");

	return writeCG(outfile, CG);
}
