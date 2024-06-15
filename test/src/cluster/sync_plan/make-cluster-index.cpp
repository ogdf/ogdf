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

#include <json.hpp>
#include <sstream>

#include <getopt.h>

#include "return.h"
#include "utils/Clusters.h"
#include "utils/Preprocess.h"

using nlohmann::json;

using namespace std;
using namespace ogdf;
static std::string progname;

int clusterCrossings(const ClusterGraph& CG) {
	int sum = 0;
	for (edge e : CG.constGraph().edges) {
		int u = CG.clusterOf(e->source())->depth();
		int v = CG.clusterOf(e->target())->depth();
		int p = CG.commonCluster(e->source(), e->target())->depth();
		sum += u + v - (2 * p);
	}
	return sum;
}

int main(int argc, char* argv[]) {
	progname = argv[0];
	for (int i = 1; i < argc; ++i) {
		Graph G;
		ClusterGraph CG(G);
		CG.setUpdateDepth(true);
		string file = argv[i];
		int read_ret = readCG(file, Logger::Level::Default, CG, G);
		if (read_ret != SUCCESS) {
			return read_ret;
		}

		NodeArray<int> components(G);
		int connected_components = connectedComponents(G, components);

		stringstream cluster_tree;
		printClusters(CG.rootCluster(), cluster_tree);

		json j {
				{"file", file},
				{"nodes", G.numberOfNodes()},
				{"edges", G.numberOfEdges()},
				{"components", connected_components},
				{"clusters", CG.numberOfClusters()},
				{"cluster-tree", cluster_tree.str()},
				{"cluster-depth", CG.treeDepth()},
				{"cluster-crossing", clusterCrossings(CG)},
				{"planar", isPlanar(G)},
				{"can-preprocess", canPreprocessClusterGraph(CG, G)},
				{"cconnected", isCConnected(CG)},
		};
		cout << j << endl;
	}
	return SUCCESS;
}
