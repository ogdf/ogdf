/*
 * $Revision: 3837 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 15:19:30 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief TLP format parser utility declaration.
 *
 * \author Łukasz Hanuszczak
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.txt in the root directory of the OGDF installation for details.
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
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/

#ifdef _MSC_VER
#pragma once
#endif


#ifndef OGDF_TLP_PARSER_H
#define OGDF_TLP_PARSER_H


#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/cluster/ClusterGraph.h>
#include <ogdf/cluster/ClusterGraphAttributes.h>
#include <ogdf/fileformats/Tlp.h>
#include <ogdf/fileformats/TlpLexer.h>

#include <iostream>
#include <string>
#include <sstream>
#include <map>


namespace ogdf {

namespace tlp {


class Parser {
private:
	typedef std::vector<Token>::const_iterator Iterator;
	std::map<int, node> m_idNode;
	std::map<int, edge> m_idEdge;

	std::istream &m_istream;
	Iterator m_begin, m_end;

	bool readEdge(Graph &G);
	bool readNodes(Graph &G, ClusterGraph *C, cluster c);
	bool readCluster(Graph &G, ClusterGraph *C, cluster c);
	bool readProperty(Graph &G, GraphAttributes *GA);

	bool readClusterStatement(
		Graph &G, ClusterGraph *C, cluster c);
	bool readStatement(
		Graph &G, GraphAttributes *GA, ClusterGraph *C);
	bool readPropertyStatement(
		GraphAttributes *GA, const Attribute &attr,
		NodeArray<bool> &nodeDone, std::string &nodeDefault,
		EdgeArray<bool> &edgeDone, std::string &edgeDefault);

	bool readGraph(Graph &G, GraphAttributes *GA, ClusterGraph *C);

	inline bool applyNodes(
		Graph &G, ClusterGraph *C, cluster c,
		const std::string &str);
	inline void tokenError(const char *str, bool got = true);
	inline void tokenError(const std::string &str, bool got = true);

public:
	Parser(std::istream &is);

	bool read(Graph &G) {
		return readGraph(G, NULL, NULL);
	}

	bool read(Graph &G, GraphAttributes &GA) {
		return readGraph(G, &GA, NULL);
	}

	bool read(Graph &G, ClusterGraph &C) {
		return readGraph(G, NULL, &C);
	}

	bool read(Graph &G, ClusterGraph &C, ClusterGraphAttributes &CA) {
		return readGraph(G, &CA, &C);
	}
};


} // end namespace tlp

} // end namespace ogdf


#endif
