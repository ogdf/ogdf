/** \file
 * \brief Implementation of the TsplibXmlParser.
 *
 * \author Finn Stutzenstein
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

#include <ogdf/fileformats/TsplibXmlParser.h>

namespace ogdf {

TsplibXmlParser::TsplibXmlParser(std::istream& in) {
	std::string error;
	m_hasError = !load(in, error);
	if (m_hasError) {
		GraphIO::logger.lout() << error << std::endl;
	}
}

bool TsplibXmlParser::load(std::istream& in, std::string& error) {
	pugi::xml_parse_result result = m_xml.load(in);

	if (!result) {
		error = "XML parser error: " + std::string(result.description());
		return false;
	}

	pugi::xml_node tspTag = m_xml.child("travellingSalesmanProblemInstance");
	if (!tspTag) {
		error = "File root tag is not a <travellingSalesmanProblemInstance>.";
		return false;
	}

	m_graphTag = tspTag.child("graph");
	if (!m_graphTag) {
		error = "<graph> tag not found.";
		return false;
	}

	return true;
}

bool TsplibXmlParser::read(Graph& G, GraphAttributes* GA) {
	if (m_hasError) {
		return false;
	}

	// count "vertex" tags and create nodes
	int n = std::distance(m_graphTag.children("vertex").begin(), m_graphTag.children("vertex").end());
	G.clear();
	Array<node> indexToNode(n);
	for (int i = 0; i < n; i++) {
		indexToNode[i] = G.newNode(i);
	}

	// 2D double array for costs.
	std::vector<std::vector<double>> costs(n, std::vector<double>(n, std::numeric_limits<double>::quiet_NaN()));

	// read costs
	int sourceVertexIndex = 0;
	for (pugi::xml_node& vertex : m_graphTag.children("vertex")) {
		// <edge cost="0.500000000000000e+00">1</edge>
		for (pugi::xml_node& edgeNode : vertex.children("edge")) {
			int targetVertexIndex;
			try {
				targetVertexIndex = std::stoi(edgeNode.child_value());
			} catch(...) {
				GraphIO::logger.lout() << "invalid vertex index encountered: " << edgeNode.child_value() << std::endl;
				return false;
			}

			if (sourceVertexIndex == targetVertexIndex) {
				continue; // skip self loops
			}

			if (targetVertexIndex >= n || targetVertexIndex < 0) {
				GraphIO::logger.lout() << "invalid vertex index encountered: " << targetVertexIndex << std::endl;
				return false;
			}

			pugi::xml_attribute costAttr = edgeNode.first_attribute();
			if (costAttr == nullptr) {
				GraphIO::logger.lout() << "no cost attribute found for source node " << sourceVertexIndex
					<< " and target node " << targetVertexIndex << std::endl;
				return false;
			}
			if (std::string(costAttr.name()) != "cost") {
				GraphIO::logger.lout() << "invalid attribute encountered: " << costAttr.name() << std::endl;
				return false;
			}

			costs[sourceVertexIndex][targetVertexIndex] = costAttr.as_double();
		}
		sourceVertexIndex++;
	}

	// the graph is directed if any costs for arcs (i, j) and (j, i) are different or not set (NaN).
	bool directed = false;
	for (int i = 0; i < n; i++) {
		for (int j = i+1; j < n; j++) {
			if (
				std::isnan(costs[i][j]) != std::isnan(costs[j][i]) ||
				(!std::isnan(costs[i][j]) && costs[i][j] != costs[j][i])
			) {
				directed = true;
				break;
			}
		}
	}

	// Add edges: For the directed case only those which are present in the matrix.
	// In the undirected case, only one arc is added for each edge.
	if (directed) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i == j || std::isnan(costs[i][j])) {
					continue;
				}
				edge e = G.newEdge(indexToNode[i], indexToNode[j]);
				if (GA) {
					GA->doubleWeight(e) = costs[i][j];
				}
			}
		}
	} else {
		for (int i = 0; i < n; i++) {
			for (int j = i+1; j < n; j++) {
				if (!std::isnan(costs[i][j])) {
					edge e = G.newEdge(indexToNode[i], indexToNode[j]);
					if (GA) {
						GA->doubleWeight(e) = costs[i][j];
					}
				}
			}
		}
	}
	if (GA) {
		GA->directed() = directed;
	}

	return true;
}

bool TsplibXmlParser::read(Graph &G) {
	return read(G, nullptr);
}

bool TsplibXmlParser::read(Graph &G, GraphAttributes &GA) {
	GA.init(G, GraphAttributes::edgeDoubleWeight);
	return read(G, &GA);
}

}
