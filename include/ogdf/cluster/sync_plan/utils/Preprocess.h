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
#pragma once

#include <ogdf/basic/Logger.h>
#include <ogdf/cluster/ClusterGraph.h>

extern ogdf::Logger preprocessLog;

ogdf::SList<ogdf::node> findSmallNodes(const ogdf::ClusterGraph& C, const ogdf::Graph& G);

ogdf::SList<ogdf::node> findDeg2Nodes(const ogdf::ClusterGraph& C, const ogdf::Graph& G);

ogdf::SList<ogdf::cluster> findDisconnectedClusters(const ogdf::ClusterGraph& C,
		const ogdf::Graph& G, ogdf::ClusterArray<ogdf::node>* centre = nullptr);

ogdf::SList<ogdf::cluster> findSmallClusters(const ogdf::ClusterGraph& C, const ogdf::Graph& G);

bool removeSmallNodes(const ogdf::ClusterGraph& C, ogdf::Graph& G);

bool unsplitDeg2Nodes(const ogdf::ClusterGraph& C, ogdf::Graph& G);

bool disconnectedClustersToStars(ogdf::ClusterGraph& C, ogdf::Graph& G);

bool removeSmallClusters(ogdf::ClusterGraph& C, ogdf::Graph& G);

//! Preprocessing from HananiTutteCPlanarity::preprocessing (Gutwenger, Mutzel, Schaefer: Practical Experience with Hanani-Tutte for Testing c-Planarity)
bool preprocessClusterGraph(ogdf::ClusterGraph& C, ogdf::Graph& G);

bool canPreprocessClusterGraph(const ogdf::ClusterGraph& C, const ogdf::Graph& G);
