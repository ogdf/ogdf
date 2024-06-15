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

#include <ogdf/basic/Graph.h>
#include <ogdf/cluster/ClusterGraph.h>
#include <ogdf/cluster/sync_plan/utils/Bijection.h>
#include <ogdf/decomposition/BCTree.h>

#include <ostream>

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
