#pragma once

#include <ogdf/cluster/ClusterGraphAttributes.h>
#include <ogdf/fileformats/GraphIO.h>

#include <chrono>
#include <json.hpp>

#include <getopt.h>

#include "PQPlanarity.h"
#include "PQPlanarityOptions.h"
#include "PipeOrder.h"

using nlohmann::json;

static const int PQPLANAR = 0x00;
static const int SUCCESS = 0x00;
static const int NOT_PQPLANAR = 0x09;
static const int FAILURE = 0x09;
static const int ERROR_ABACUS = 0x20;
static const int ERROR_OPTIONS = 0x40;
static const int ERROR_IO = 0x41;
static const int ERROR_ASSERT = 0x42;
static const int ERROR_COMB_EMB = 0x43;

Logger::Level parseLogLevel() {
	try {
		switch (std::stoi(optarg)) {
		case 0:
		default:
			return ogdf::Logger::Level::Minor;
		case 1:
			return ogdf::Logger::Level::Medium;
		case 2:
			return ogdf::Logger::Level::Default;
		case 3:
			return ogdf::Logger::Level::High;
		case 4:
			return ogdf::Logger::Level::Alarm;
		case 5:
			return ogdf::Logger::Level::Force;
		}
	} catch (const std::exception& e) {
		throw invalid_option(string_format("Log level must be an integer! (%s)", e.what()));
	}
}

template<class... Types>
int readCG(const string& infile, Logger::Level log_level, Types&... args) {
	Logger::globalLogLevel(Logger::Level::Alarm);
	bool read_success;
	if (infile[0] == '-') {
		read_success = GraphIO::read(args..., std::cin);
	} else {
		read_success = GraphIO::read(args..., infile);
	}
	Logger::globalLogLevel(log_level);
	if (!read_success) {
		if (infile[0] == '-') {
			Logger::slout(Logger::Level::Alarm)
					<< "IO Error: Couldn't read graph from stdin!" << std::endl;
		} else {
			std::ifstream is(infile);
			if (!is.good()) {
				Logger::slout(Logger::Level::Alarm)
						<< "IO Error: Couldn't read file " << infile << "!" << std::endl;
			} else {
				Logger::slout(Logger::Level::Alarm)
						<< "IO Error: Couldn't read graph from file " << infile << "!" << std::endl;
			}
		}
		return ERROR_IO;
	}
	return SUCCESS;
}

int writeCG(const string& outfile, const ClusterGraphAttributes& CGA) {
	if (!outfile.empty()) {
		bool write_success;
#ifdef OGDF_DEBUG
		CGA.constClusterGraph().consistencyCheck();
#endif
		if (outfile[0] == '-' && (outfile.length() == 1 || outfile[1] == '.')) {
			const GraphIO::FileType* t = GraphIO::getFileType(outfile);
			if (t == nullptr) {
				write_success = GraphIO::writeGML(CGA, std::cout);
			} else if (t->cluster_attr_writer_func != nullptr) {
				write_success = t->cluster_attr_writer_func(CGA, std::cout);
			} else if (t->cluster_writer_func != nullptr) {
				write_success = t->cluster_writer_func(CGA.constClusterGraph(), std::cout);
			} else {
				Logger::slout(Logger::Level::Alarm)
						<< "Invalid Options: File format for " << outfile
						<< " does not allow writing clusters!" << std::endl;
				return ERROR_OPTIONS;
			}
		} else {
			write_success = GraphIO::write(CGA, outfile);
		}
		if (!write_success) {
			Logger::slout(Logger::Level::Alarm)
					<< "IO Error: Couldn't write graph to " << outfile << "!" << std::endl;
			return ERROR_IO;
		}
	}
	return SUCCESS;
}

int writeCG(const string& outfile, const ClusterGraph& CG) {
	if (!outfile.empty()) {
		bool write_success;
#ifdef OGDF_DEBUG
		CG.consistencyCheck();
#endif
		if (outfile[0] == '-' && (outfile.length() == 1 || outfile[1] == '.')) {
			const GraphIO::FileType* t = GraphIO::getFileType(outfile);
			if (t == nullptr) {
				write_success = GraphIO::writeGML(CG, std::cout);
			} else if (t->cluster_writer_func != nullptr) {
				write_success = t->cluster_writer_func(CG, std::cout);
			} else {
				Logger::slout(Logger::Level::Alarm)
						<< "Invalid Options: File format for " << outfile
						<< " does not allow writing clusters!" << std::endl;
				return ERROR_OPTIONS;
			}
		} else {
			write_success = GraphIO::write(CG, outfile);
		}
		if (!write_success) {
			Logger::slout(Logger::Level::Alarm)
					<< "IO Error: Couldn't write graph to " << outfile << "!" << std::endl;
			return ERROR_IO;
		}
	}
	return SUCCESS;
}

struct PQPlanConf {
	bool allow_contract = true, intersect_trees = true, random_order = false,
		 contract_first = false, invert_degree = false, invert_contract = false, batch_spqr = false;

#define PQPlanConf_KEYS          \
	TOGGLE('c', allow_contract)  \
	TOGGLE('i', intersect_trees) \
	TOGGLE('r', random_order)    \
	TOGGLE('a', invert_degree)   \
	TOGGLE('b', contract_first)  \
	TOGGLE('s', invert_contract) \
	TOGGLE('p', batch_spqr)

	PipeQueue* getOrder(PQPlanarity* PQ) const {
		if (random_order) {
			return new PipeQueueRandom();
		} else if (contract_first) {
			return new PipeQueueByDegreePreferContract(PQ, invert_degree, invert_contract);
		} else {
			return new PipeQueueByDegree(invert_degree);
		}
	}

	friend std::ostream& operator<<(std::ostream& os, const PQPlanConf& conf) {
		os << "allow_contract: " << conf.allow_contract
		   << " intersect_trees: " << conf.intersect_trees << " random_order: " << conf.random_order
		   << " contract_first: " << conf.contract_first << " invert_degree: " << conf.invert_degree
		   << " invert_contract: " << conf.invert_contract << " batch_spqr: " << conf.batch_spqr;
		return os;
	}

	string getID() const {
		const PQPlanConf def;
		string ret;
#define TOGGLE(key, var)  \
	if (var != def.var) { \
		ret += "-";       \
		ret += key;       \
	}
		PQPlanConf_KEYS
#undef TOGGLE
				return ret;
	}

	NLOHMANN_DEFINE_TYPE_INTRUSIVE(PQPlanConf, allow_contract, intersect_trees, random_order,
			contract_first, invert_degree, invert_contract, batch_spqr)
};

string typeOfPipe(const PQPlanarity& pq, const Pipe& p) {
	if (p.degree() <= 3) {
		return "small";
	} else if (pq.getComponents().isCutVertex(p.node1)) {
		if (pq.getComponents().isCutVertex(p.node2)) {
			return "cut_cut";
		} else {
			return "bicon_cut";
		}
	} else {
		if (pq.getComponents().isCutVertex(p.node2)) {
			return "bicon_cut";
		} else {
			if (pq.getComponents().connectedId(p.node1) == pq.getComponents().connectedId(p.node2)) {
				return "bicon_bicon_same_cc";
			} else {
				return "bicon_bicon_diff_cc";
			}
		}
	}
}

void pqPlanStats(const PQPlanarity& pq, json& s) {
	s["nodes"] = pq.G->numberOfNodes();
	s["edges"] = pq.G->numberOfEdges();
	s["partitions"] = pq.partitions.partitionCount();
	s["q_vertices"] = pq.partitions.qVertexCount();
	s["pipes"] = pq.matchings.getPipeCount();
	s["pipes_degrees"] = 0;
	if (pq.matchings.getPipeCount() > 0) {
		s["top_pipe"] = pq.matchings.getTopPipe().degree();
		json degrees_types = {
				{"<=3", {}},
				{"=4", {}},
				{"=5", {}},
				{"<10", {}},
				{"<20", {}},
				{"<50", {}},
				{"<100", {}},
				{">=100", {}},
		};
		json types, degrees;
		int max_deg = 0;
		string max_deg_type = "";
		for (const Pipe& p : pq.matchings) {
			s["pipes_degrees"] = int(s["pipes_degrees"]) + p.degree();
			string type = typeOfPipe(pq, p);
			if (p.degree() > max_deg) {
				max_deg = p.degree();
				max_deg_type = type;
			}
			if (!types.contains(type)) {
				types[type] = 1;
				degrees[type] = p.degree();
			} else {
				types[type] = int(types[type]) + 1;
				degrees[type] = int(degrees[type]) + p.degree();
			}
			string pdeg;
			if (p.degree() <= 3) {
				pdeg = "<=3";
			} else if (p.degree() == 4) {
				pdeg = "=4";
			} else if (p.degree() == 5) {
				pdeg = "=5";
			} else if (p.degree() < 10) {
				pdeg = "<10";
			} else if (p.degree() < 20) {
				pdeg = "<20";
			} else if (p.degree() < 50) {
				pdeg = "<50";
			} else if (p.degree() < 100) {
				pdeg = "<100";
			} else {
				pdeg = ">=100";
			}
			if (!degrees_types[pdeg].contains(type)) {
				degrees_types[pdeg][type] = 1;
			} else {
				degrees_types[pdeg][type] = int(degrees_types[pdeg][type]) + 1;
			}
		}
		s["pipe_types"] = types;
		s["pipe_types_degrees"] = degrees;
		s["pipe_degrees_types"] = degrees_types;
		s["max_pipe_deg"] = max_deg;
		s["max_pipe_deg_type"] = max_deg_type;
	}
	s["connected"] = pq.getComponents().connectedCount();
	s["bc_edges"] = pq.getComponents().bcTree().numberOfEdges();
	s["bicon"] = s["cut"] = 0;
	s["bicon_degrees"] = {
			{"0", 0},
			{"1", 0},
			{"2", 0},
			{"+", 0},
	};
	s["bicon_sizes"] = {
			{"<10", 0},
			{"<50", 0},
			{"<100", 0},
			{"<500", 0},
			{"<1000", 0},
			{">=1000", 0},
	};
	int bicon_sum_size = 0, bicon_max_size = 0;
	for (node tn : pq.getComponents().bcTree().nodes) {
		if (pq.getComponents().isCutComponent(tn)) {
			s["cut"] = s["cut"].get<int>() + 1;
		} else {
			s["bicon"] = s["bicon"].get<int>() + 1;

			int size = pq.getComponents().bcSize(tn);
			bicon_sum_size += size;
			bicon_max_size = max(bicon_max_size, size);
			string ssize;
			if (size < 10) {
				ssize = "<10";
			} else if (size < 50) {
				ssize = "<50";
			} else if (size < 100) {
				ssize = "<100";
			} else if (size < 500) {
				ssize = "<500";
			} else if (size < 1000) {
				ssize = "<1000";
			} else {
				ssize = ">=1000";
			}
			s["bicon_sizes"][ssize] = s["bicon_sizes"][ssize].get<int>() + 1;

			string deg;
			switch (tn->degree()) {
			case 0:
				deg = "0";
				break;
			case 1:
				deg = "1";
				break;
			case 2:
				deg = "2";
				break;
			default:
				deg = "+";
				break;
			}
			s["bicon_degrees"][deg] = s["bicon_degrees"][deg].get<int>() + 1;
		}
	}
	s["bicon_sum_size"] = bicon_sum_size;
	s["bicon_max_size"] = bicon_max_size;
}
