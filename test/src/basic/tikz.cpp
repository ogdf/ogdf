/** \file
 * \brief Tests for the ogdf::TikzWriter
 *
 * \author Hendrik Brueckler
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

#include <ogdf/basic/graph_generators.h>
#include <ogdf/fileformats/GraphIO.h>

#include <regex>

#include <testing.h>

/**
 * @brief Generates random attributes relevant for TikZ drawing. Only fills enabled attributes.
 *
 * @param GA the attributes to fill with random values
 */
void randomTikzAttributes(GraphAttributes& GA) {
	const Graph& graph = GA.constGraph();
	long attr = GA.attributes();

	GA.directed() = false;
	for (node v : graph.nodes) {
		if (attr & GraphAttributes::nodeLabelPosition) {
			GA.xLabel(v) = randomNumber(1, 10);
			GA.yLabel(v) = randomNumber(1, 10);
		}
		if (attr & GraphAttributes::nodeStyle) {
			GA.strokeColor(v) = randomNumber(0, 1) ? Color::Name::Peru : Color::Name::Whitesmoke;
			GA.strokeType(v) = randomNumber(0, 1) ? StrokeType::Dashdotdot : StrokeType::Solid;
			GA.strokeWidth(v) = randomDouble(0.1, 1.5);
			GA.fillColor(v) = randomNumber(0, 1) ? Color::Name::Lime : Color::Name::Gainsboro;
		}
		if (attr & GraphAttributes::nodeLabel) {
			GA.label(v) = "Node " + to_string(v->index());
		}
		if (attr & GraphAttributes::nodeGraphics) {
			GA.x(v) = randomNumber(1, 100);
			GA.y(v) = randomNumber(1, 100);
			GA.width(v) = randomNumber(5, 10);
			GA.height(v) = randomNumber(5, 10);
			if (GA.width(v) == GA.height(v)) {
				int shape = randomNumber(1, 3);
				switch (shape) {
				case 1:
					GA.shape(v) = Shape::Pentagon;
					break;
				case 2:
					GA.shape(v) = Shape::Hexagon;
					break;
				case 3:
					GA.shape(v) = Shape::Octagon;
					break;
				}
			} else {
				int shape = randomNumber(1, 10);
				switch (shape) {
				case 1:
					GA.shape(v) = Shape::Rect;
					break;
				case 2:
					GA.shape(v) = Shape::RoundedRect;
					break;
				case 3:
					GA.shape(v) = Shape::Ellipse;
					break;
				case 4:
					GA.shape(v) = Shape::InvParallelogram;
					break;
				case 5:
					GA.shape(v) = Shape::Parallelogram;
					break;
				case 6:
					GA.shape(v) = Shape::Triangle;
					break;
				case 7:
					GA.shape(v) = Shape::InvTriangle;
					break;
				case 8:
					GA.shape(v) = Shape::Trapeze;
					break;
				case 9:
					GA.shape(v) = Shape::InvTrapeze;
					break;
				case 10:
					GA.shape(v) = Shape::Rhomb;
					break;
				}
			}
		}
	}
	for (edge e : graph.edges) {
		if (attr & GraphAttributes::edgeLabel) {
			GA.label(e) = "Edge " + to_string(e->index());
		}
		if (attr & GraphAttributes::edgeArrow) {
			GA.arrowType(e) = randomNumber(0, 1) ? EdgeArrow::Both : EdgeArrow::First;
		}
		if (attr & GraphAttributes::edgeStyle) {
			GA.strokeColor(e) = randomNumber(0, 1) ? Color::Name::Papayawhip : Color::Name::Cornsilk;
			GA.strokeType(e) = randomNumber(0, 1) ? StrokeType::Dashdotdot : StrokeType::Dashdot;
			GA.strokeWidth(e) = randomDouble(0.1, 1.5);
		}
	}
}

/**
 * @brief Generates a string that contains the whole drawTikz output
 *
 * @param attr the attributes to use for drawTikz
 * @return std::string drawTikz output
 */
std::string drawTikzToString(GraphAttributes attr) {
	std::ostringstream write;
	AssertThat(GraphIO::drawTikz(attr, write), IsTrue());
	std::string tikzString = write.str();

	// Replace newlines by tabs because regex is awful.
	std::replace(tikzString.begin(), tikzString.end(), '\n', '\t');
	return tikzString;
}

/**
 * @brief Get the number of matches of a regex vs a string.
 *
 * @param regex regular expression
 * @param s string to match against
 * @return int number of matches
 */
int countRegexMatches(const std::regex& regex, const std::string& s) {
	ptrdiff_t numMatches =
			std::distance(std::sregex_iterator(s.begin(), s.end(), regex), std::sregex_iterator());
	return static_cast<int>(numMatches);
}

go_bandit([]() {
	describe("GraphIO", []() {
		describe("drawTikz", []() {
			// regex for anything but newline (which we replace by tabs, so anything)
			const string any = ".*";
			// regex for anything but semicolon (they separate \node's and \path's)
			const string anyButSC = "[^;]*";

			Graph graph;
			// regexes will cause stack overflow for too big graphs, i.e. too long strings
			// this is not a limitation of drawTikz()
			const int numberOfNodes = 10;

			before_each([&]() {
				graph.clear();
				randomSimpleGraph(graph, numberOfNodes, 2 * numberOfNodes);
			});

			it("produces well-formed LaTeX/TikZ output", [&]() {
				std::string doc = drawTikzToString(GraphAttributes(graph));

				std::regex docClass("\\\\documentclass\\{standalone\\}");
				std::regex tikzIncluded("\\\\usepackage\\{tikz\\}");
				std::regex docStructure("\\\\begin\\{document\\}" + any + "\\\\begin\\{tikzpicture\\}"
						+ any + "\\\\end\\{tikzpicture\\}" + any + "\\\\end\\{document\\}");
				std::regex node("\\\\node" + anyButSC + "\\(Node\\d+\\)" + anyButSC + ";");
				std::regex edge("\\\\path" + anyButSC + ";");

				AssertThat(countRegexMatches(docClass, doc), Equals(1));
				AssertThat(countRegexMatches(tikzIncluded, doc), Equals(1));
				AssertThat(countRegexMatches(docStructure, doc), Equals(1));
				AssertThat(countRegexMatches(node, doc), Equals(graph.numberOfNodes()));
				AssertThat(countRegexMatches(edge, doc), Equals(graph.numberOfEdges()));
			});

			it("handles its supported attributes correctly and stores style info in header", [&]() {
				GraphAttributes attr(graph, GraphAttributes::all);
				randomTikzAttributes(attr);
				std::string doc = drawTikzToString(attr);

				std::regex tikzPictureHeader("\\[yscale = -1\\.0" + any + "width/\\.style" + any
						+ "height/\\.style" + any + "size/\\.style" + any + "nodelabel/\\.style"
						+ any + "shiftednodelabel/\\.style" + any + "edgelabel/\\.style" + any
						+ "nodestyle0/\\.style" + any + "edgestyle0/\\.style");
				std::regex node("\\\\node" + anyButSC + "nodestyle\\d+" + anyButSC + ";");
				std::regex edge("\\\\path" + anyButSC + "edgestyle\\d+" + anyButSC + ";");
				AssertThat(countRegexMatches(tikzPictureHeader, doc), Equals(1));
				AssertThat(countRegexMatches(node, doc), Equals(graph.numberOfNodes()));
				AssertThat(countRegexMatches(edge, doc), Equals(graph.numberOfEdges()));
			});

			it("sets correct node positions including label", [&]() {
				GraphAttributes attr(graph, GraphAttributes::all);
				randomTikzAttributes(attr);
				for (node v : graph.nodes) {
					attr.x(v) = v->index();
					attr.y(v) = v->index();
					attr.label(v) = "NodeLabel" + to_string(v->index());
					attr.xLabel(v) = v->index();
					attr.yLabel(v) = v->index();
				}
				std::string doc = drawTikzToString(attr);

				std::regex nodeshifted("\\\\node" + anyButSC + "(\\d+)(pt|mm|cm|in|ex|em|mu)"
						+ anyButSC + "\\1\\2" + anyButSC + "NodeLabel\\1" + anyButSC
						+ "\\(Node\\1\\) at \\(\\1\\2, \\1\\2\\)(" + anyButSC + ");");
				std::regex node("\\\\node" + anyButSC + "NodeLabel0" + anyButSC
						+ "\\(Node0\\) at \\(0(pt|mm|cm|in|ex|em|mu), 0\\1\\)" + anyButSC + ";");
				AssertThat(countRegexMatches(nodeshifted, doc), Equals(graph.numberOfNodes() - 1));
				AssertThat(countRegexMatches(node, doc), Equals(1));
			});

			it("sets correct edge connections and edge label", [&]() {
				GraphAttributes attr(graph, GraphAttributes::all);
				randomTikzAttributes(attr);
				for (edge e : graph.edges) {
					attr.label(e) = "EdgeLabelSource" + to_string(e->source()->index()) + "Target"
							+ to_string(e->target()->index());
				}
				std::string doc = drawTikzToString(attr);

				std::regex edge("\\\\path" + anyButSC + "EdgeLabelSource(\\d+)Target(\\d+)\\}"
						+ anyButSC + "\\(Node\\1\\)" + anyButSC + "\\(Node\\2\\)" + anyButSC + ";");
				AssertThat(countRegexMatches(edge, doc), Equals(graph.numberOfEdges()));
			});

			it("sets correct style and outsources uniform node style into header", [&]() {
				GraphAttributes attr(graph,
						GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics
								| GraphAttributes::nodeStyle | GraphAttributes::edgeStyle);
				for (node v : graph.nodes) {
					attr.shape(v) = Shape::Rect;
					attr.strokeColor(v) = Color(1, 33, 7);
					attr.strokeType(v) = StrokeType::Dashdotdot;
					attr.strokeWidth(v) = static_cast<float>(3.142);
				}
				for (edge e : graph.edges) {
					attr.strokeColor(e) = Color(42, 42, 42);
					attr.strokeType(e) = StrokeType::Dash;
					attr.strokeWidth(e) = static_cast<float>(4.2069);
				}
				std::string doc = drawTikzToString(attr);

				std::regex nodestyle(
						"nodestyle0/.style = \\{rectangle, draw = \\{rgb, 255: red,1; green,33; blue,7\\}, dash dot dot, line width = 3\\.142");
				std::regex everynodestyle("every node/.append style = \\{nodestyle0");
				std::regex edgestyle(
						"edgestyle0/.style = \\{draw = \\{rgb, 255: red,42; green,42; blue,42\\}, dashed, line width = 4\\.2069");
				std::regex node("\\\\node" + anyButSC + "\\(Node\\d+\\)" + anyButSC + ";");
				std::regex nodewithstyle("\\\\node\\[" + anyButSC + "nodestyle" + anyButSC + ";");
				std::regex edge("\\\\path" + anyButSC + "edgestyle" + anyButSC + ";");
				AssertThat(countRegexMatches(nodestyle, doc), Equals(1));
				AssertThat(countRegexMatches(everynodestyle, doc), Equals(1));
				AssertThat(countRegexMatches(edgestyle, doc), Equals(1));
				AssertThat(countRegexMatches(node, doc), Equals(graph.numberOfNodes()));
				AssertThat(countRegexMatches(nodewithstyle, doc), Equals(0));
				AssertThat(countRegexMatches(edge, doc), Equals(graph.numberOfEdges()));
			});

			it("draws clusters", [&]() {
				ClusterGraph clusterGraph(graph);
				randomClusterGraph(clusterGraph, graph, 10);
				ClusterGraphAttributes attr(clusterGraph,
						ClusterGraphAttributes::clusterGraphics | GraphAttributes::nodeGraphics
								| GraphAttributes::nodeStyle);

				for (node v : graph.nodes) {
					attr.shape(v) = Shape::Octagon;
				}

				std::ostringstream write;
				AssertThat(GraphIO::drawTikz(attr, write), IsTrue());
				std::string doc = write.str();

				std::regex cluster("\\\\node" + anyButSC + "\\(Cluster\\d+\\)" + anyButSC + ";");
				AssertThat(countRegexMatches(cluster, doc),
						Equals(clusterGraph.numberOfClusters() - 1));
			});

			it("avoids connecting an edge bend that lies inside the target node", [&]() {
				GraphAttributes attr(graph, GraphAttributes::all);
				randomTikzAttributes(attr);
				for (node v : graph.nodes) {
					attr.x(v) = v->index() * 100;
					attr.y(v) = v->index() * 100;
				}
				for (edge e : graph.edges) {
					DPolyline bends;
					DPoint midpoint = 0.5 * (attr.point(e->source()) + attr.point(e->target()));
					bends.emplaceBack(attr.point(e->source()));
					bends.emplaceBack(midpoint);
					bends.emplaceBack(attr.point(e->target()));
					attr.bends(e) = bends;
				}
				std::string doc = drawTikzToString(attr);

				std::regex edge("\\\\path" + anyButSC + "\\(Node\\d+\\)" + anyButSC + ";");
				AssertThat(countRegexMatches(edge, doc), Equals(0));
			});
		});
	});
});
