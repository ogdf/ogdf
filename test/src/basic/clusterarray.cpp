/** \file
 * \brief Tests for ogdf::ClusterArray
 *
 * \author Matthias Pfretzschner
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

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/RegisteredSet.h>
#include <ogdf/basic/graph_generators/clustering.h>
#include <ogdf/basic/graph_generators/randomized.h>
#include <ogdf/cluster/ClusterGraph.h>
#include <ogdf/cluster/ClusterSet.h>

#include <functional>
#include <list>
#include <string>

#include "array_helper.h"
#include <testing.h>

go_bandit([]() {
	auto allClusters = [](const ClusterGraph& C, List<cluster>& list) { C.allClusters(list); };

	auto chooseCluster = [](const ClusterGraph& C) {
		while (true) {
			cluster c = *chooseIteratorFrom<internal::GraphObjectContainer<ClusterElement>, cluster>(
					const_cast<internal::GraphObjectContainer<ClusterElement>&>(C.clusters));
			if (c != C.rootCluster()) {
				return c;
			}
		}
	};
	auto createCluster = [](ClusterGraph& C) { return C.createEmptyCluster(); };
	auto deleteCluster = [](ClusterGraph& C, cluster c) { return C.delCluster(c); };
	auto clearClusters = [](ClusterGraph& C) { return C.clear(); };

	std::list<Graph> graphs;
	auto init = [&](ClusterGraph& C) {
		graphs.emplace_back();
		randomSimpleConnectedGraph(graphs.back(), 42, 168);
		C.init(graphs.back());
		randomClustering(C, 42);
	};

	runBasicArrayTests<ClusterGraph, ClusterArray, cluster>( //
			"ClusterArray", init, chooseCluster, allClusters, createCluster);

	graphs.clear();
	runBasicSetTests<ClusterGraph, ClusterSet, cluster>("ClusterSet", init, chooseCluster,
			allClusters, createCluster, deleteCluster, clearClusters);
});
