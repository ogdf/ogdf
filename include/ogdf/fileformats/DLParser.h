/*
 * $Revision: 3837 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 15:19:30 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration of UCINET DL format parser class.
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

#ifndef OGDF_DL_PARSER_H
#define OGDF_DL_PARSER_H


#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphAttributes.h>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>


namespace ogdf {


class DLParser {
private:
	std::istream &m_istream;
	bool m_initialized;

	int m_nodes;
	enum { fullmatrix, edgelist, nodelist } m_format;
	bool m_embedded;

	std::vector<node> m_nodeId; // For constant-time index to node mapping.
	std::map<std::string, node> m_nodeLabel; // For embedded label mode.

	static inline void toUpper(std::string &str) {
		std::transform(str.begin(), str.end(), str.begin(), toupper);
	}

	static inline void toLower(std::string &str) {
		std::transform(str.begin(), str.end(), str.begin(), tolower);
	}

	inline bool fineId(int vid) {
		return 0 < vid && vid < static_cast<int>(m_nodeId.size());
	}

	inline node requestLabel(
		GraphAttributes *GA, node &nextFree,
		const std::string &label);

	void init();
	bool initGraph(Graph &G);

	bool readMatrix(Graph &G, GraphAttributes *GA);
	bool readEdgeList(Graph &G, GraphAttributes *GA);
	bool readNodeList(Graph &G, GraphAttributes *GA);
	bool readEmbeddedMatrix(Graph &G, GraphAttributes *GA);
	bool readEmbeddedEdgeList(Graph &G, GraphAttributes *GA);
	bool readEmbeddedNodeList(Graph &G, GraphAttributes *GA);

	bool readAssignment(
		Graph &G,
		const std::string &lhs, const std::string &rhs);

	bool readData(Graph &G, GraphAttributes *GA);
	bool readWithLabels(Graph &G, GraphAttributes *GA);
	bool readStatements(Graph &G, GraphAttributes *GA);
	bool readGraph(Graph &G, GraphAttributes *GA);

public:
	DLParser(std::istream &is);

	bool read(Graph &G) {
		return readGraph(G, NULL);
	}

	bool read(Graph &G, GraphAttributes &GA) {
		return readGraph(G, &GA);
	}
};


} // end namespace ogdf


#endif
