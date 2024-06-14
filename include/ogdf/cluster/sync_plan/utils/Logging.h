#pragma once

#include <ogdf/basic/Graph.h>
#include <ogdf/cluster/ClusterGraph.h>
#include <ogdf/decomposition/BCTree.h>

#include <ostream>

#include "Bijection.h"

std::string to_string(const std::function<std::ostream&(std::ostream&)>& func);

std::ostream& operator<<(std::ostream& os, const std::function<std::ostream&(std::ostream&)>& func);

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const std::pair<T1, T2>& pair) {
	return os << "(" << pair.first << ", " << pair.second << ")";
}

std::ostream& operator<<(std::ostream& os, const ogdf::Graph& G);

std::ostream& operator<<(std::ostream& os, const ogdf::ClusterGraph& CG);

std::ostream& operator<<(std::ostream& os, const ogdf::BCTree::BNodeType& obj);

std::ostream& operator<<(std::ostream& os, const ogdf::BCTree::GNodeType& obj);

#define CONTAINER_PRINTER(NAME)                                                          \
	template<typename Container>                                                         \
	struct NAME {                                                                        \
		const Container& container;                                                      \
		explicit NAME(const Container& container) : container(container) { }             \
		template<typename ContainerT>                                                    \
		friend std::ostream& operator<<(std::ostream& os, const NAME<ContainerT>& inst); \
	}

CONTAINER_PRINTER(printIncidentEdges);

CONTAINER_PRINTER(printEdges);

CONTAINER_PRINTER(printBijection);

CONTAINER_PRINTER(printFrozenBijection);

template<typename Container>
std::ostream& operator<<(std::ostream& os, const printIncidentEdges<Container>& inst) {
	for (ogdf::adjEntry adj : inst.container) {
		os << "e" << adj->theEdge()->index() << " (" << (adj->isSource() ? ">" : "<") << "n"
		   << adj->twinNode()->index() << "), ";
	}
	return os;
}

template<>
std::ostream& operator<<(std::ostream& os, const printIncidentEdges<PipeBij>& inst);

template<typename Container>
std::ostream& operator<<(std::ostream& os, const printEdges<Container>& inst) {
	for (ogdf::adjEntry adj : inst.container) {
		os << "e" << adj->theEdge()->index() << " (n" << adj->theNode()->index()
		   << (adj->isSource() ? "->" : "<-") << "n" << adj->twinNode()->index() << "), ";
	}
	return os;
};

template<>
std::ostream& operator<<(std::ostream& os, const printEdges<PipeBij>& inst);

template<typename Container>
std::ostream& operator<<(std::ostream& os, const printBijection<Container>& inst) {
	bool first = true;
	for (const std::pair<ogdf::adjEntry, ogdf::adjEntry>& pair : inst.container) {
		if (first) {
			first = false;
		} else {
			os << " ";
		}
		os << "(";
		if (pair.first == nullptr) {
			os << "NULL";
		} else {
			os << "n" << pair.first->twinNode()->index() << (pair.first->isSource() ? "<" : ">")
			   << " e" << pair.first->theEdge()->index();
		}
		os << " = ";
		if (pair.second == nullptr) {
			os << "NULL";
		} else {
			os << "e" << pair.second->theEdge()->index() << " "
			   << (pair.second->isSource() ? ">" : "<") << "n" << pair.second->twinNode()->index();
		}
		os << ")";
	}
	return os;
};

template<typename Container>
std::ostream& operator<<(std::ostream& os, const printFrozenBijection<Container>& inst) {
	bool first = true;
	for (const auto& pair : inst.container) {
		os << (first ? "" : " ") << "(e" << pair.first << " = e" << pair.second << ")";
		first = false;
	}
	return os;
};
