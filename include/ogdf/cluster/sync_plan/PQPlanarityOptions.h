#pragma once

#include <ogdf/basic/GraphAttributes.h>

#include <json.hpp>
#include <stdexcept>

#include "PQPlanarity.h"
#include "utils/Bijection.h"

using namespace ogdf;

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

class invalid_option : public std::runtime_error {
public:
	explicit invalid_option(const string& message) : runtime_error(message) { }
};

template<typename... Args>
node getNode(const Graph& G, const GraphAttributes& GA, const nlohmann::json& id,
		const std::string& format, const Args... args) {
	throw new runtime_error(
			"not implemented"); // too inefficient to be usable, see PQPlanOptions::applyConfigJSON for a better approach
	// node n = nullptr;
	// if (id.is_number_integer()) {
	//     n = G.chooseNode([&id](node n) -> bool { return n->index() == id; });
	// } else if (id.is_string()) {
	//     n = G.chooseNode([&GA, &id](node n) -> bool { return GA.label(n) == id; });
	// }
	// if (n == nullptr) {
	//     throw invalid_option(string_format(format, id.dump().c_str(), args...));
	// }
	// return n;
}

template<typename... Args>
edge getEdge(const Graph& G, const GraphAttributes& GA, const nlohmann::json& id,
		const std::string& format, const Args... args) {
	throw new runtime_error("not implemented");
	// edge e = nullptr;
	// if (id.is_number_integer()) {
	//     e = G.chooseEdge([&id](edge e) -> bool { return e->index() == id; });
	// } else if (id.is_string()) {
	//     e = G.chooseEdge([&GA, &id](edge e) -> bool { return GA.label(e) == id; });
	// }
	// if (e == nullptr) {
	//     throw invalid_option(string_format(format, id.dump().c_str(), args...));
	// }
	// return e;
}

class PQPlanOptions {
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

	void apply(Graph& G, GraphAttributes& GA, PQPlanarity& pq);

	static void applyConfigJSON(Graph& G, GraphAttributes& GA, PQPlanarity& pq, nlohmann::json& j);

	template<typename NodeLabeler = std::function<int(node)>,
			typename EdgeLabeler = std::function<int(edge)>>
	static void generateConfigJSON(
			PQPlanarity& pq, nlohmann::json& j,
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
