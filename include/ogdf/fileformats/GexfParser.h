/*
 * $Revision: 3837 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 15:19:30 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration of GEXF format reading utilities.
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

#ifndef OGDF_GEXF_PARSER_H
#define OGDF_GEXF_PARSER_H


#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/cluster/ClusterGraph.h>
#include <ogdf/cluster/ClusterGraphAttributes.h>
#include <ogdf/basic/HashArray.h>
#include <ogdf/fileformats/XmlParser.h>

#include <iostream>
#include <sstream>

namespace ogdf {

namespace gexf {


class Parser {
private:
	XmlParser m_xml;
	XmlTagObject *m_graphTag, *m_nodesTag, *m_edgesTag;

	HashArray<std::string, node> m_nodeId;
	HashArray<std::string, cluster> m_clusterId;

	HashArray<std::string, std::string> m_nodeAttr, m_edgeAttr;

	bool init();
	bool readNodes(Graph &G, GraphAttributes *GA);
	bool readEdges(Graph &G, ClusterGraph *C, GraphAttributes *GA);
	bool readCluster(
		Graph &G, ClusterGraph &C, ClusterGraphAttributes *CA,
		cluster rootCluster,
		const XmlTagObject &rootTag);
	bool readAttributes(
		GraphAttributes &GA, node v,
		const XmlTagObject &nodeTag);
	bool readAttributes(
		GraphAttributes &GA, edge e,
		const XmlTagObject &edgeTag);

	static void error(const XmlTagObject &tag, const std::string msg);

public:
	Parser(std::istream &is);

	bool read(Graph &G);
	bool read(Graph &G, GraphAttributes &GA);
	bool read(Graph &G, ClusterGraph &C);
	bool read(Graph &G, ClusterGraph &C, ClusterGraphAttributes &CA);
};


}

} // end namespace ogdf


#endif
