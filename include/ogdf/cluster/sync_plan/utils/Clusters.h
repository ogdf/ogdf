#pragma once

#include <ogdf/cluster/ClusterGraph.h>

#include <ostream>

using namespace ogdf;

struct RandomClusterConfig {
	int max_nodes_in_cluster = 0;
	double prob_no_further_node = 0.1;
	double prob_no_further_cluster = 0.0;
	int max_clusters = 0;
	int min_root_nodes = 0;
	bool cconnected = false;
	int timeout = 0;

	double expected_nodes() const { return 1.0 / prob_no_further_node; }

	void expected_nodes(double n) { prob_no_further_node = 1.0 / n; }

	friend std::ostream& operator<<(std::ostream& os, const RandomClusterConfig& config) {
		os << "max_nodes_in_cluster: " << config.max_nodes_in_cluster
		   << " prob_no_further_node: " << config.prob_no_further_node << " ("
		   << config.expected_nodes() << ")"
		   << " prob_no_further_cluster: " << config.prob_no_further_cluster << " ("
		   << 1.0 / config.prob_no_further_cluster << ")"
		   << " max_clusters: " << config.max_clusters
		   << " min_root_nodes: " << config.min_root_nodes << " timeout: " << config.timeout;
		return os;
	}
};

bool makeClusters(ClusterGraph& CG, RandomClusterConfig& config);

bool isClusterPlanarEmbedding(const ClusterGraph& CG);

void clusterBorderToEdges(
		const ClusterGraph& CG, Graph& G,
		EdgeArray<List<std::pair<adjEntry, cluster>>>* subdivisions = nullptr,
		const std::function<edge(edge)>& translate = [](edge e) -> edge { return e; });

std::ostream& printClusters(cluster c, std::ostream& s);

void printCG(const ClusterGraph& CG, const std::string& type = "");
