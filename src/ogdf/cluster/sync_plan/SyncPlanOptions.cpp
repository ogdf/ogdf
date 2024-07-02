/** \file
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

#if 0

#	include <ogdf/cluster/sync_plan/SyncPlan.h>
#	include <ogdf/cluster/sync_plan/SyncPlanAttributes.h>
#	include <ogdf/cluster/sync_plan/SyncPlanOptions.h>

#	include <unordered_map>

using json = nlohmann::json;

using namespace ogdf::sync_plan::internal;
namespace ogdf::sync_plan {

void SyncPlanOptions::parseOptionPipes(char* optarg) {
	try {
		size_t delim;
		int f = std::stoi(optarg, &delim);
		if (optarg[delim] != ',') {
			throw std::invalid_argument("values must be delimited by comma");
		}
		size_t offset = delim + 1;
		int s = std::stoi(optarg + offset, &delim);
		if (optarg[delim + offset] != 0) {
			throw std::invalid_argument("too many values");
		}
		pipes.emplaceBack(f, s);
	} catch (const std::exception& e) {
		throw invalid_option(string_format(
				"Pipe specification must have the form first_index,second_index! (%s)", e.what()));
	}
}

void SyncPlanOptions::parseOptionPartitions(char* optarg) {
	try {
		size_t offset = 0, delim = 0;
		auto list = partitions.emplaceBack();
		while (true) {
			(*list).pushBack(std::stoi(optarg + offset, &delim));
			char next = optarg[delim + offset];
			if (next == 0) {
				break;
			} else if (next != ',') {
				throw std::invalid_argument("values must be delimited by comma");
			}
			offset += delim + 1;
		}
	} catch (const std::exception& e) {
		throw invalid_option(string_format("Partition specification must have the form "
										   "first_index,second_index,third_index,...! (%s)",
				e.what()));
	}
}

void SyncPlanOptions::parseOptionEmbedding(char* optarg) {
	try {
		size_t delim;
		auto order = orders.emplaceBack(std::stoi(optarg, &delim));
		if (optarg[delim] != ':') {
			throw std::invalid_argument("node id must be delimited by colon");
		}
		size_t offset = delim + 1;
		while (true) {
			size_t len = std::strcspn(optarg + offset, ",");
			if (len < 1) {
				throw std::invalid_argument("edge name missing");
			}
			(*order).edges.emplaceBack(optarg + offset, len);
			char next = optarg[offset + len];
			if (next == 0) {
				break;
			} else if (next != ',') {
				throw std::invalid_argument("values must be delimited by comma");
			}
			offset += len + 1;
		}
	} catch (const std::exception& e) {
		throw invalid_option(string_format("Embedding specification must have the form "
										   "node_index:edge1_name,edge2_name,...! (%s)",
				e.what()));
	}
}

void SyncPlanOptions::apply(Graph& G, GraphAttributes& GA, SyncPlan& pq) {
	for (auto pipe : pipes) {
		node first = getNode(G, GA, (json)pipe.first,
				"Couldn't find first node %s in pipe (%d, %d)!", pipe.first, pipe.second);
		node second = getNode(G, GA, pipe.second, "Couldn't find second node %s in pipe (%d, %d)!",
				pipe.first, pipe.second);
		pq.matchings.matchNodes(first, second);
	}

	for (auto& partition : partitions) {
		partition.quicksort();
		int part = NO_PARTITION;
		for (node n : G.nodes) {
			if (n->index() != partition.front()) {
				continue;
			}
			if (part == NO_PARTITION) {
				part = pq.partitions.makeQVertex(n);
			} else {
				pq.partitions.makeQVertex(n, part);
			}
			partition.popFront();
			if (partition.empty()) {
				break;
			}
		}
		if (!partition.empty()) {
			throw invalid_option(string_format("Couldn't find node %d for partition %d!",
					partition.front(), part));
		}
	}

	for (auto& order : orders) {
		node n = getNode(G, GA, order.node, "Couldn't find embedded node %s!");
		List<adjEntry> adjs;
		for (auto edge_name : order.edges) {
			adjs.pushBack(getEdge(G, GA, edge_name, "Couldn't find embedded edge %s!")->getAdj(n));
		}
		G.sort(n, adjs);
	}

	OGDF_ASSERT(pq.consistency.consistencyCheck());
}

void SyncPlanOptions::applyConfigJSON(Graph& G, GraphAttributes& GA, SyncPlan& pq, json& j) {
	unordered_map<json, node> nodes;
	unordered_map<json, edge> edges;
	nodes.reserve(G.numberOfNodes() * 2);
	edges.reserve(G.numberOfEdges() * 2);
	for (node n : G.nodes) {
		nodes[(json)n->index()] = n;
		nodes[(json)GA.label(n)] = n;
	}
	for (edge e : G.edges) {
		edges[(json)e->index()] = e;
		edges[(json)GA.label(e)] = e;
	}
	for (json& pipe : j["pipes"]) {
		node first =
				nodes.at(pipe[0]); // "Couldn't find first node %s in pipe %s!", pipe.dump().c_str());
		node second = nodes.at(
				pipe[1]); // "Couldn't find second node %s in pipe %s!", pipe.dump().c_str());
		pq.matchings.matchNodes(first, second);
		if (pipe.size() == 2) {
			continue;
		}
		List<adjEntry> first_adj;
		List<adjEntry> second_adj;
		if (pipe.size() == 3) {
			for (auto& k_v : pipe[2].items()) {
				// note: key is always string
				first_adj.pushBack(edges.at(k_v.key())->getAdj(first));
				// "Couldn't find first embedded edge %s for pipe %s!", pipe.dump().c_str())
				second_adj.pushBack(edges.at(k_v.value())->getAdj(second));
				// "Couldn't find second embedded edge %s for pipe %s!", pipe.dump().c_str())
			}
		} else if (pipe.size() == 4) {
			for (json& e_id : pipe[2]) {
				first_adj.pushBack(edges.at(e_id)->getAdj(first));
				// "Couldn't find first embedded edge %s for pipe %s!", pipe.dump().c_str())
			}
			for (json& e_id : pipe[3]) {
				second_adj.pushBack(edges.at(e_id)->getAdj(second));
				// "Couldn't find second embedded edge %s for pipe %s!", pipe.dump().c_str())
			}
		}
		G.sort(first, first_adj);
		G.sort(second, second_adj);
		G.reverseAdjEdges(second); // json contains pipes in bijection order, which has one side reversed
	}

	for (json& partition : j["partitions"]) {
		int part = NO_PARTITION;
		if (partition.is_object()) {
			for (auto& k_v : partition.items()) {
				node n = nodes.at(k_v.key()); // "Couldn't find node %s for partition %d!", part);
				if (part == NO_PARTITION) {
					part = pq.partitions.makeQVertex(n);
				} else {
					pq.partitions.makeQVertex(n, part);
				}

				List<adjEntry> adjs;
				for (json& e_id : k_v.value()) {
					adjs.pushBack(edges.at(e_id)->getAdj(n));
					// "Couldn't find incident edge %s of node %s in partition %d!", k_v.key().c_str(), part)
				}
				G.sort(n, adjs);
			}
		} else {
			for (json& n_id : partition) {
				node n = nodes.at(n_id); // "Couldn't find node %s for partition %d!", part);
				if (part == NO_PARTITION) {
					part = pq.partitions.makeQVertex(n);
				} else {
					pq.partitions.makeQVertex(n, part);
				}
			}
		}
	}

	for (json& embedding : j["embeddings"]) {
		if (embedding.is_object()) {
			for (auto& k_v : embedding.items()) {
				node n = nodes.at(k_v.key()); // "Couldn't find node %s for embedding!");

				List<adjEntry> adjs;
				for (json& e_id : k_v.value()) {
					adjs.pushBack(edges.at(e_id)->getAdj(
							n)); // "Couldn't find incident edge %s of node %s!", k_v.key().c_str())
				}
				G.sort(n, adjs);
			}
		} else {
			for (json& n_id_e_ids : embedding) {
				node n = nodes.at(n_id_e_ids[0]); // "Couldn't find node %s for embedding!");

				List<adjEntry> adjs;
				for (json& e_id : n_id_e_ids[1]) {
					adjs.pushBack(edges.at(e_id)->getAdj(n));
					// "Couldn't find incident edge %s of node %s!", n_id_e_ids[0].dump().c_str())
				}
				G.sort(n, adjs);
			}
		}
	}
}

}

#endif
