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
