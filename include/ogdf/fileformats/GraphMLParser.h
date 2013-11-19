/*
 * $Revision: 3837 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 15:19:30 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Declarations for GraphML Parser
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

#ifndef OGDF_GRAPHML_PARSER_H
#define OGDF_GRAPHML_PARSER_H


#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/cluster/ClusterGraph.h>
#include <ogdf/cluster/ClusterGraphAttributes.h>

#include <ogdf/basic/HashArray.h>
#include <ogdf/basic/List.h>
#include <ogdf/fileformats/XmlParser.h>

#include <sstream>


namespace ogdf {

class GraphMLParser {
private:
	XmlParser m_xml;
	XmlTagObject *m_graphTag; // "Almost root" tag.

	HashArray<string, node> m_nodeId; // Maps GraphML node id to Graph node.
	HashArray<string, string> m_attrName; // Maps attribute id to its name.

	bool readData(
		GraphAttributes &GA,
		const node &v, const XmlTagObject &nodeData);
	bool readData(
		GraphAttributes &GA,
		const edge &e, const XmlTagObject &edgeData);
	bool readData(
		ClusterGraphAttributes &CA,
		const cluster &c, const XmlTagObject &clusterData);

	// Finds all data-keys for given element and calls appropiate "readData".
	template <typename A, typename T>
	bool readAttributes(A &GA, const T &elem, const XmlTagObject &elemTag) {
		List<XmlTagObject *> dataTags;
		elemTag.findSonXmlTagObjectByName("data", dataTags);

		forall_listiterators(XmlTagObject *, it, dataTags) {
			const bool result = readData(GA, elem, **it);

			if(!result) {
				return false;
			}
		}

		return true;
	}

	bool readNodes(
		Graph &G, GraphAttributes *GA,
		const XmlTagObject &rootTag);
	bool readEdges(
		Graph &G, GraphAttributes *GA,
		const XmlTagObject &rootTag);
	bool readClusters(
		Graph &G, ClusterGraph &C, ClusterGraphAttributes *CA,
		const cluster &rootCluster, const XmlTagObject &clusterRoot);

	bool m_error;

public:
	GraphMLParser(istream &in);
	~GraphMLParser();

	bool read(Graph &G);
	bool read(Graph &G, GraphAttributes &GA);
	bool read(Graph &G, ClusterGraph &C);
	bool read(Graph &G, ClusterGraph &C, ClusterGraphAttributes &CA);
};


} // end namespace ogdf


#endif
