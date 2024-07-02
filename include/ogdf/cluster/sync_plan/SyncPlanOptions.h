/** \file
 * \brief Class for (de)serializing a SyncPlan instance from/to JSON, enable if you have nlohmann json.hpp available.
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

// enable if you have nlohmann json.hpp available
#if 0

#	include <ogdf/basic/GraphAttributes.h>
#	include <ogdf/cluster/sync_plan/SyncPlan.h>
#	include <ogdf/cluster/sync_plan/utils/Bijection.h>

#	include <json.hpp>

namespace ogdf::sync_plan {

namespace internal {
template<typename... Args>
std::string string_format(const std::string& format, const Args... args) {
	size_t size = snprintf(nullptr, 0, format.c_str(), args...) + 1; // Extra space for '\0'
	if (size <= 0) {
		throw std::runtime_error("Error during formatting.");
	}
	std::unique_ptr<char[]> buf(new char[size]);
	snprintf(buf.get(), size, format.c_str(), args...);
	return std::string(buf.get(), buf.get() + size - 1); // We don't want the '\0' inside
}
}

//! Allows (de)serializing a SyncPlan instance from/to JSON.
class OGDF_EXPORT SyncPlanOptions {
	struct EdgeOrder {
		int node;
		List<string> edges;

		EdgeOrder(int node) : node(node) { }
	};

	FrozenPipeBij pipes;
	List<List<int>> partitions;
	List<EdgeOrder> orders;

public:
	void parseOptionPipes(char* optarg);

	void parseOptionPartitions(char* optarg);

	void parseOptionEmbedding(char* optarg);

	void apply(Graph& G, GraphAttributes& GA, SyncPlan& pq);

	static void applyConfigJSON(Graph& G, GraphAttributes& GA, SyncPlan& pq, nlohmann::json& j);

	template<typename NodeLabeler = std::function<int(node)>,
			typename EdgeLabeler = std::function<int(edge)>>
	static void generateConfigJSON(
			SyncPlan& pq, nlohmann::json& j,
			const NodeLabeler& nl = [](node n) -> int { return n->index(); },
			const EdgeLabeler& el = [](edge e) -> int { return e->index(); }) {
		using nlohmann::json;
		json pipes = json::array();
		for (const auto& pipe : pq.matchings) {
			json adj1 = json::array(), adj2 = json::array();
			for (auto edges : pq.matchings.getIncidentEdgeBijection(pipe.node1)) {
				adj1 += el(edges.first->theEdge());
				adj2 += el(edges.second->theEdge());
			}
			pipes.push_back(json::array({nl(pipe.node1), nl(pipe.node2), adj1, adj2}));
		}
		j["pipes"] = pipes;

		json partitions = json::array();
		json embeddings = json::array();
		for (int p = 0; p < pq.partitions.partitionCount(); p++) {
			json ids = json::array();
			for (node n : pq.partitions.nodesInPartition(p)) {
				ids += nl(n);

				json adj = json::array();
				for (adjEntry a : n->adjEntries) {
					adj += el(a->theEdge());
				}
				embeddings.push_back(json::array({nl(n), adj}));
			}
			partitions.push_back(ids);
		}
		j["partitions"] = partitions;
		j["embeddings"] = embeddings;
	}
};

}

#endif
