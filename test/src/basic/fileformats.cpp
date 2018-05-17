/** \file
 * \brief Tests for fileformat reading and writing using GraphIO,
 *   only graphs without attributes
 *
 * \author Stephan Beyer, Tilo Wiedera
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

#include <algorithm>
#include <string>
#include <regex>
#include <unordered_map>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/EpsilonTest.h>
#include <ogdf/basic/graph_generators.h>
#include <ogdf/fileformats/GraphIO.h>
#include <resources.h>
#include <ogdf/basic/simple_graph_alg.h>

using namespace ogdf;
using namespace bandit;
using namespace std;

void assertSeemsEqual(const Graph &G1, const Graph &G2) {
	AssertThat(G1.numberOfNodes(), Equals(G2.numberOfNodes()));
	AssertThat(G1.numberOfEdges(), Equals(G2.numberOfEdges()));

	Array<int> counter1, counter2;
	degreeDistribution(G1, counter1);
	degreeDistribution(G2, counter2);

	AssertThat(counter1.size(), Equals(counter2.size()));
	AssertThat(counter1.low(), Equals(counter2.low()));

	for(int i = counter1.low(); i < counter1.high(); i++) {
		AssertThat(counter1[i], Equals(counter2[i]));
	}
}
void establishNodeMapping(NodeArray<node> &map1to2, const GraphAttributes &GA1,
                          const GraphAttributes &GA2)
{
	const Graph &G1 = GA1.constGraph();
	const Graph &G2 = GA2.constGraph();
	std::vector<node> mapIndexToNode;
	mapIndexToNode.resize(G1.numberOfNodes());
	for(node v1 : G1.nodes) {
		int x1;
		if(GA1.has(GraphAttributes::nodeGraphics)) {
			x1 = GA1.x(v1) - 1;
		} else {
			AssertThat(GA1.has(GraphAttributes::nodeLabel), IsTrue());
			x1 = atoi(GA1.label(v1).c_str());
		}
		AssertThat(mapIndexToNode[x1], Equals(nullptr));
		mapIndexToNode[x1] = v1;
	}
	for(node v2 : G2.nodes) {
		int x2;
		if(GA1.has(GraphAttributes::nodeGraphics)) {
			x2 = GA2.x(v2) - 1;
		} else {
			AssertThat(GA1.has(GraphAttributes::nodeLabel), IsTrue());
			x2 = atoi(GA2.label(v2).c_str());
		}
		AssertThat(map1to2[mapIndexToNode[x2]], Equals(nullptr));
		map1to2[mapIndexToNode[x2]] = v2;
	}
}

void assertEqualGAs(const GraphAttributes &GA1, const GraphAttributes &GA2) {
	const Graph &G1 = GA1.constGraph();
	const Graph &G2 = GA2.constGraph();
	NodeArray<node> map1to2(G1, nullptr);
	AssertThat(GA1.attributes(), Equals(GA2.attributes()));
	AssertThat(G1.numberOfNodes(), Equals(G2.numberOfNodes()));
	AssertThat(G1.numberOfEdges(), Equals(G2.numberOfEdges()));

	establishNodeMapping(map1to2, GA1, GA2);

	constexpr double delta { 0.5 };

	for(node v : G1.nodes) {
		if(GA1.has(GraphAttributes::nodeGraphics)) {
			AssertThat(GA2.x(map1to2[v]), Equals(GA1.x(v)));
			AssertThat(GA2.y(map1to2[v]), EqualsWithDelta(GA1.y(v), delta));
			if(GA1.has(GraphAttributes::threeD)) {
				AssertThat(GA2.z(map1to2[v]), EqualsWithDelta(GA1.z(v), delta));
			}
			AssertThat(GA2.width(map1to2[v]), EqualsWithDelta(GA1.width(v), delta));
			AssertThat(GA2.height(map1to2[v]), EqualsWithDelta(GA1.height(v), delta));
			AssertThat(GA2.shape(map1to2[v]), Equals(GA1.shape(v)));
		}
		if(GA1.has(GraphAttributes::nodeId)) {
			AssertThat(GA2.idNode(map1to2[v]), Equals(GA1.idNode(v)));
		}
		if(GA1.has(GraphAttributes::nodeLabel)) {
			AssertThat(GA2.label(map1to2[v]), Equals(GA1.label(v)));
		}
		if(GA1.has(GraphAttributes::nodeLabelPosition)) {
			AssertThat(GA2.xLabel(map1to2[v]), EqualsWithDelta(GA1.xLabel(v), delta));
			AssertThat(GA2.yLabel(map1to2[v]), EqualsWithDelta(GA1.yLabel(v), delta));
			if(GA1.has(GraphAttributes::threeD)) {
				AssertThat(GA2.zLabel(map1to2[v]), Equals(GA1.zLabel(v)));
			}
		}
		if(GA1.has(GraphAttributes::nodeStyle)) {
			AssertThat(GA2.fillColor(map1to2[v]), Equals(GA1.fillColor(v)));
			AssertThat(GA2.strokeColor(map1to2[v]), Equals(GA1.strokeColor(v)));
			AssertThat(GA2.strokeType(map1to2[v]), Equals(GA1.strokeType(v)));
			AssertThat(GA2.strokeWidth(map1to2[v]), Equals(GA1.strokeWidth(v)));
			AssertThat(GA2.fillPattern(map1to2[v]), Equals(GA1.fillPattern(v)));
		}
		if(GA1.has(GraphAttributes::nodeTemplate)) {
			AssertThat(GA2.templateNode(map1to2[v]), Equals(GA1.templateNode(v)));
		}
		if(GA1.has(GraphAttributes::nodeType)) {
			AssertThat(int(GA2.type(map1to2[v])), Equals(int(GA1.type(v))));
		}
		if(GA1.has(GraphAttributes::nodeWeight)) {
			AssertThat(GA2.weight(map1to2[v]), EqualsWithDelta(GA1.weight(v), delta));
		}
	}
	for(edge e : G1.edges) {
		edge e2 = G2.searchEdge(map1to2[e->source()], map1to2[e->target()]);
		AssertThat(e2, !Equals(nullptr));

		if(GA1.has(GraphAttributes::edgeArrow)) {
			AssertThat(GA2.arrowType(e2), Equals(GA1.arrowType(e)));
		}
		if(GA1.has(GraphAttributes::edgeGraphics)) {
			AssertThat(GA2.bends(e2), Equals(GA1.bends(e)));
		}
		if(GA1.has(GraphAttributes::edgeLabel)) {
			AssertThat(GA2.label(e2), Equals(GA1.label(e)));
		}
		if(GA1.has(GraphAttributes::edgeType)) {
			AssertThat(GA2.type(e2), Equals(GA1.type(e)));
		}
		if(GA1.has(GraphAttributes::edgeStyle)) {
			AssertThat(GA2.strokeColor(e2), Equals(GA1.strokeColor(e)));
			AssertThat(GA2.strokeType(e2), Equals(GA1.strokeType(e)));
			AssertThat(GA2.strokeWidth(e2), EqualsWithDelta(GA1.strokeWidth(e), delta));
		}
		if(GA1.has(GraphAttributes::edgeDoubleWeight)) {
			AssertThat(GA2.doubleWeight(e2), EqualsWithDelta(GA1.doubleWeight(e), delta));
		} else if(GA1.has(GraphAttributes::edgeIntWeight)) {
			AssertThat(GA2.intWeight(e2), Equals(GA1.intWeight(e)));
		}
	}
}

void describeSTPonlyGraph() {
	describe("unweighted STP", [] {
		for_each_file("fileformats/stp/valid", [&](const ResourceFile* file){
			it("successfully parses " + file->fullPath(), [&] {
				Graph graph;
				stringstream is{file->data()};
				AssertThat(GraphIO::readSTP(graph, is), IsTrue());
			});
		});

		for_each_file("fileformats/stp/invalid", [&](const ResourceFile* file){
			it("detects errors in " + file->fullPath(), [&] {
				Graph graph;
				stringstream is{file->data()};
				AssertThat(GraphIO::readSTP(graph, is), IsFalse());
			});
		});
	});
}

template<typename T>
void describeSTP(const string &typeName) {
	describe("STP for " + typeName, [] {
		for(int i = 4; i < 1024; i *= 2) {
			it("stores and loads an instance of size " + to_string(i), [&] {
				std::ostringstream writeStream;

				EdgeWeightedGraph<T> graph;
				List<node> terminals;
				NodeArray<bool> isTerminal(graph, false);

				randomGraph(graph, i, (i*(i-1))/2);
				for(node v : graph.nodes) {
					if(randomDouble(0, 1) > .5) {
						terminals.pushBack(v);
						isTerminal(v) = true;
					}
				}
				for (edge e : graph.edges) {
					graph.setWeight(e, (T) randomDouble(0, 1000));
				}

				string myComment = "";
				if (randomDouble(0, 1) > .5) {
					myComment += "Name \"MyRandomInstance\"\n";
					myComment += "Creator \"Tilo Wiedera\"\n";
				}
				GraphIO::writeSTP(graph, terminals, writeStream, myComment);

				EdgeWeightedGraph<T> readGraph;
				List<node> readTerminals;
				NodeArray<bool> readIsTerminal;

				std::istringstream readStream(writeStream.str());
				AssertThat(GraphIO::readSTP(readGraph, readTerminals, readIsTerminal, readStream), Equals(true));

				AssertThat(readGraph.numberOfNodes(), Equals(graph.numberOfNodes()));
				AssertThat(readGraph.numberOfEdges(), Equals(graph.numberOfEdges()));
				AssertThat(readTerminals.size(), Equals(terminals.size()));
				for(node v : readGraph.nodes) {
					AssertThat(readIsTerminal[v], Equals(readTerminals.search(v).valid()));
				}
			});
		}

		it("clears the graph", [&](){
			EdgeWeightedGraph<T> writeGraph;
			List<node> terminals;
			std::ostringstream write;
			AssertThat(GraphIO::writeSTP(writeGraph, terminals, write), Equals(true));

			EdgeWeightedGraph<T> readGraph;
			customGraph(readGraph, 2, {{0, 1}});
			NodeArray<bool> isTerminal(readGraph, true);
			terminals.pushBack(readGraph.firstNode());
			std::istringstream read(write.str());
			AssertThat(GraphIO::readSTP(readGraph, terminals, isTerminal, read), Equals(true));
			AssertThat(readGraph.empty(), IsTrue());
			AssertThat(terminals.empty(), IsTrue());
			AssertThat(isTerminal.begin(), Equals(isTerminal.end()));
		});

		for_each_file("fileformats/stp/valid", [&](const ResourceFile* file){
			it("successfully parses " + file->fullPath(), [&] {
				EdgeWeightedGraph<T> graph;
				List<node> terminals;
				NodeArray<bool> isTerminal;
				stringstream is{file->data()};
				AssertThat(GraphIO::readSTP(graph, terminals, isTerminal, is), IsTrue());

				AssertThat(graph.numberOfNodes(), IsGreaterThan(0));
				AssertThat(graph.numberOfEdges(), IsGreaterThan(0));
				AssertThat(terminals.size(), IsGreaterThan(0));

				int terminalCounter = 0;
				for(node v : graph.nodes) {
					terminalCounter += isTerminal[v];
				}

				AssertThat(terminalCounter, Equals(terminals.size()));
			});
		});

		for_each_file("fileformats/stp/invalid", [&](const ResourceFile* file){
			it("detects errors in " + file->fullPath(), [&](){
				EdgeWeightedGraph<T> graph;
				List<node> terminals;
				NodeArray<bool> isTerminal;
				stringstream is{file->data()};
				AssertThat(GraphIO::readSTP(graph, terminals, isTerminal, is), IsFalse());
			});
		});
	});
}

template<typename T>
void describeDMF(const string &typeName) {
	describe("DMF for " + typeName, [] {
		const void* nullPointer = nullptr;

		for_each_file("fileformats/dmf/valid", [&](const ResourceFile* file) {
			it("reads " + file->fullPath(), [&]() {
				Graph graph;
				EdgeArray<T> weights;
				node source;
				node sink;

				stringstream is{file->data()};
				AssertThat(GraphIO::readDMF(graph, weights, source, sink, is), IsTrue());
				AssertThat(graph.numberOfNodes(), IsGreaterThan(1));
				AssertThat(weights.valid(), IsTrue());
				AssertThat(source, Is().Not().EqualTo(nullPointer));
				AssertThat(sink, Is().Not().EqualTo(nullPointer));
#ifdef OGDF_DEBUG
				AssertThat(source->graphOf(), Equals(&graph));
				AssertThat(sink->graphOf(), Equals(&graph));
#endif
				AssertThat(source, Is().Not().EqualTo(sink));

				for(edge e : graph.edges) {
					AssertThat(weights(e) > 0, IsTrue());
				}
			});
		});

		for_each_file("fileformats/dmf/invalid", [&](const ResourceFile* file) {
			it("reads " + file->fullPath(), [&]() {
				Graph graph;
				EdgeArray<T> weights(graph, 0);
				node source;
				node sink;
				stringstream is{file->data()};
				AssertThat(GraphIO::readDMF(graph, weights, source, sink, is), IsFalse());
			});
		});

		it("writes and reads a random graph", [&]() {
			Graph graph;
			EdgeArray<T> weights(graph, 0);
			node source;
			node sink;

			randomGraph(graph, 42, 189);
			source = graph.chooseNode();
			sink = graph.chooseNode([&](node v) { return v != source; });

			T sum = 0;
			for(edge e : graph.edges) {
				T cap = static_cast<T>(randomDoubleNormal(10, 5));
				if(cap < 0) {
					cap *= -1;
				}
				weights(e) = cap;
				sum += cap;
			}

			std::ostringstream writeStream;

			AssertThat(GraphIO::writeDMF(graph, weights, source, sink, writeStream), IsTrue());

			Graph readGraph;
			EdgeArray<T> readWeights(readGraph, 0);
			node readSource = nullptr;
			node readSink = nullptr;

			std::istringstream readStream(writeStream.str());
			AssertThat(GraphIO::readDMF(readGraph, readWeights, readSource, readSink, readStream), IsTrue());

			AssertThat(readGraph.numberOfNodes(), Equals(graph.numberOfNodes()));
			AssertThat(readGraph.numberOfEdges(), Equals(graph.numberOfEdges()));
			AssertThat(readSource, Is().Not().EqualTo(nullPointer));
			AssertThat(readSink, Is().Not().EqualTo(nullPointer));
#ifdef OGDF_DEBUG
				AssertThat(readSource->graphOf(), Equals(&readGraph));
				AssertThat(readSink->graphOf(), Equals(&readGraph));
#endif
			AssertThat(readSource->degree(), Equals(source->degree()));
			AssertThat(readSink->degree(), Equals(sink->degree()));

			T readSum = 0;
			for(edge e : readGraph.edges) {
				readSum += readWeights(e);
			}

			EpsilonTest eps(1.0e-3);
			AssertThat(eps.equal(sum, readSum), IsTrue());
		});

		it("clears the graph", [&]() {
			Graph writeGraph;
			EdgeArray<T> writeWeights(writeGraph, 42);
			completeGraph(writeGraph, 3);
			node source = writeGraph.firstNode();
			node sink = writeGraph.lastNode();

			std::ostringstream write;
			AssertThat(GraphIO::writeDMF(writeGraph, writeWeights, source, sink, write), IsTrue());

			Graph readGraph;
			EdgeArray<T> readWeights(readGraph, 0);
			customGraph(readGraph, 2, {{0, 1}});
			source = nullptr;
			sink = nullptr;

			std::istringstream read(write.str());
			AssertThat(GraphIO::readDMF(readGraph, readWeights, source, sink, read), IsTrue());
			AssertThat(readGraph.numberOfNodes(), Equals(3));
			AssertThat(readGraph.numberOfEdges(), Equals(3));
			AssertThat(readWeights[readGraph.firstEdge()], Equals(42));
			AssertThat(source, !Equals(sink));
			AssertThat(source, !Equals(nullptr));
			AssertThat(sink, !Equals(nullptr));
		});
	});
}

/**
 * Used to describe a format parser and writer.
 *
 * \param name The name of the format.
 * \param reader The parse function to be tested.
 * \param writer The write function to be tested.
 * \param isXml Whether the format is based on XML.
 */
void testFormat(const std::string name, function<bool(Graph &G, istream &is)> reader,
                        function<bool(Graph &G, ostream &os)> writer, bool isXml)
{
	std::string lowerCaseName = name;
	std::transform(lowerCaseName.begin(), lowerCaseName.end(), lowerCaseName.begin(), ::tolower);

	auto errorTest = [&](const ResourceFile* file) {
		it("detects errors in " + file->fullPath(), [&]() {
			Graph graph;
			stringstream ss{file->data()};
			AssertThat(reader(graph, ss), IsFalse());
		});
	};

	auto resourceBasedTest = [&]() {
		for_each_file("fileformats/" + lowerCaseName + "/valid", [&](const ResourceFile* file) {
			it("successfully parses " + file->fullPath(), [&](){
				Graph graph;
				stringstream ss{file->data()};
				AssertThat(reader(graph, ss), IsTrue());
				AssertThat(graph.numberOfNodes(), IsGreaterThan(0));
				AssertThat(graph.numberOfEdges(), IsGreaterThan(0));
			});
		});
		for_each_file("fileformats/" + lowerCaseName + "/valid/skip", [&](const ResourceFile* file) {
			it_skip("successfully parses " + file->fullPath(), [&]() {});
		});

		for_each_file("fileformats/" + lowerCaseName + "/invalid", errorTest);
		for_each_file("fileformats/" + lowerCaseName + "/invalid/skip", [&](const ResourceFile* file) {
			it_skip("detects errors in " + file->fullPath(), [&]() {});
		});

		it("returns false if the files does not exist", [&]() {
			Graph graph;
			ifstream input;
			input.close();
			AssertThat(reader(graph, input), IsFalse());
		});
	};

	if(isXml) {
		for_each_file("fileformats/xml/invalid", errorTest);
	}

	it("detects invalid input streams", [&]() {
		Graph G;
		std::istringstream badStream;
		badStream.setstate(std::istringstream::badbit);
		AssertThat(reader(G, badStream), IsFalse());
	});

	it("detects invalid output streams", [&]() {
		Graph G;
		randomGraph(G, 10, 20);
		std::ostringstream badStream;
		badStream.setstate(std::ostringstream::badbit);
		AssertThat(writer(G, badStream), IsFalse());
	});

	resourceBasedTest();

	it("writes and reads an empty graph", [&]() {
		Graph G, Gtest;
		std::ostringstream write;
		AssertThat(writer(G, write), Equals(true));
		std::istringstream read(write.str());
		AssertThat(reader(Gtest, read), Equals(true));
		assertSeemsEqual(G, Gtest);
	});

	it("clears the graph", [&]() {
		Graph writeGraph;
		std::ostringstream write;
		AssertThat(writer(writeGraph, write), IsTrue());

		Graph readGraph;
		customGraph(readGraph, 2, {{0, 1}});
		std::istringstream read(write.str());
		AssertThat(reader(readGraph, read), IsTrue());
		AssertThat(readGraph.empty(), IsTrue());
	});

	it("writes and reads a graph of isolated nodes", [&]() {
		Graph G, Gtest;
		G.newNode(); G.newNode();
		std::ostringstream write;
		AssertThat(writer(G, write), Equals(true));
		std::istringstream read(write.str());
		AssertThat(reader(Gtest, read), Equals(true));
		assertSeemsEqual(G, Gtest);
	});

	it("writes and reads a Petersen graph", [&]() {
		Graph G, Gtest;
		petersenGraph(G, 5, 2);
		std::ostringstream write;
		AssertThat(writer(G, write), Equals(true));
		std::istringstream read(write.str());
		AssertThat(reader(Gtest, read), Equals(true));
		assertSeemsEqual(G, Gtest);
	});

	it("writes and reads a big complete graph", [&]() {
		Graph G, Gtest;
		completeGraph(G, 243);
		std::ostringstream write;
		AssertThat(writer(G, write), Equals(true));
		std::istringstream read(write.str());
		AssertThat(reader(Gtest, read), Equals(true));
		assertSeemsEqual(G, Gtest);
	});
}

/**
 * Used to describe a format parser and writer.
 *
 * \param name The name of the format.
 * \param reader The parse function to be tested.
 * \param writer The write function to be tested.
 * \param isXml Whether the format is based on XML.
 * \param alreadyDescribed Whether this function is called from inside a describe.
 */
void describeFormat(const std::string name, function<bool(Graph &G, istream &is)> reader, GraphIO::WriterFunc writer,
                    bool isXml, bool alreadyDescribed = false)
{
	if(alreadyDescribed) {
		testFormat(name, reader, writer, isXml);
	} else {
		describe(name, [&]() {
			testFormat(name, reader, writer, isXml);
		});
	}
}

/**
 * Used to describe a format parser and writer that respects GraphAttributes.
 *
 * @copydetails describeFormat
 *
 * @param readerGA The parse function respecting GraphAttributes to be tested.
 * @param writerGA The write function respecting GraphAttributes to be tested.
 * @param attr The chosen GraphAttributes.
 */
void describeGAFormat(const std::string name, GraphIO::AttrReaderFunc readerGA, GraphIO::AttrWriterFunc writerGA,
                      GraphIO::ReaderFunc reader, GraphIO::WriterFunc writer, bool isXml, long attr, bool isGEXF = false)
{
	function<bool(Graph &G, istream &is)> graphOnlyReader = [&](Graph& G,istream& is) {
		GraphAttributes GA(G, 0);
		return readerGA(GA, G, is);
	};

	function<bool(Graph &G, ostream &is)> graphOnlyWriter = [&](Graph& G,ostream& os) {
		GraphAttributes GA(G, 0);
		return writerGA(GA, os);
	};

	describeFormat(name, reader, writer, isXml, true);

	describe(name, [&](){
		describe("with GraphAttributes", [&]() {
			testFormat(name, graphOnlyReader, graphOnlyWriter, isXml);

			it("writes and reads a big graph while maintaining GraphAttributes", [&](){
				Graph graph;
				randomSimpleGraph(graph, 20, 40);
				GraphAttributes GA(graph, attr); // all attributes specified in the call activated
				for(node v : graph.nodes) {
					if(attr & GraphAttributes::nodeLabelPosition) {
						GA.xLabel(v) = v->index();
						GA.yLabel(v) = randomNumber(1, std::numeric_limits<int>::max());
						if(attr & GraphAttributes::threeD) {
							GA.zLabel(v) = randomNumber(1, std::numeric_limits<int>::max());
						}
					}
					if(attr & GraphAttributes::nodeStyle) {
						GA.strokeColor(v) = randomNumber(0,1) ? Color::Name::Peru : Color::Name::Whitesmoke;
						GA.strokeType(v) = randomNumber(0,1) ? StrokeType::Dashdotdot : StrokeType::Solid;
						GA.strokeWidth(v) = randomNumber(1, std::numeric_limits<int>::max());
						GA.fillPattern(v) = randomNumber(0, 1) ? FillPattern::Cross : FillPattern::Dense1;
						GA.fillColor(v) = randomNumber(0, 1) ? Color::Name::Blanchedalmond : Color::Name::Gainsboro;
						GA.fillBgColor(v) = randomNumber(0, 1) ? Color::Name::Mistyrose : Color::Name::Mintcream;
					}
					if(attr & GraphAttributes::threeD) {
						GA.z(v) = randomNumber(1, std::numeric_limits<int>::max());
					}
					if(attr & GraphAttributes::nodeWeight) {
						GA.weight(v) = randomNumber(1, std::numeric_limits<int>::max());
					}
					if(attr & GraphAttributes::nodeTemplate) {
						GA.templateNode(v) = to_string(randomNumber(1, std::numeric_limits<int>::max()));
					}
					if(attr & GraphAttributes::nodeType) {
						GA.type(v) = (randomNumber(0, 1) ? Graph::NodeType::dummy : Graph::NodeType::associationClass);
					}
					if(attr & GraphAttributes::nodeLabel) {
						GA.label(v) = to_string(v->index());
					}
					if(attr & GraphAttributes::nodeId) {
						GA.idNode(v) = v->index();
					}
					if(attr & GraphAttributes::nodeGraphics) {
						GA.x(v) = v->index() + 1;
						GA.y(v) = randomNumber(1, std::numeric_limits<int>::max());
						int size = randomNumber(1, 10);
						GA.width(v) = size;
						GA.height(v) = isGEXF ? size : randomNumber(1, 10);
						GA.shape(v) = (randomNumber(0, 1) ? Shape::Ellipse : Shape::Image);
					}
				}
				for(edge e : graph.edges) {
					if(attr & GraphAttributes::edgeGraphics) {
						DPolyline bends1;
						bends1.emplaceFront(randomNumber(1, std::numeric_limits<int>::max()), randomNumber(1, std::numeric_limits<int>::max()));
						GA.bends(e) = bends1;
					}
					if(attr & GraphAttributes::edgeIntWeight) {
						GA.intWeight(e) = randomNumber(2, std::numeric_limits<int>::max());
					}
					if(attr & GraphAttributes::edgeDoubleWeight) {
						GA.doubleWeight(e) = randomNumber(2, std::numeric_limits<int>::max());
					}
					if(attr & GraphAttributes::edgeLabel) {
						GA.label(e) = to_string(randomNumber(1, std::numeric_limits<int>::max()));
					}
					if(attr & GraphAttributes::edgeType) {
						GA.type(e) = randomNumber(0,1) ? Graph::EdgeType::generalization : Graph::EdgeType::association;
					}
					if(attr & GraphAttributes::edgeArrow) {
						GA.arrowType(e) = randomNumber(0,1) ? EdgeArrow::Both : EdgeArrow::First;
					}
					if(attr & GraphAttributes::edgeStyle) {
						GA.strokeColor(e) = randomNumber(0,1) ? Color::Name::Papayawhip : Color::Name::Cornsilk;
						GA.strokeType(e) = randomNumber(0,1) ? StrokeType::Dashdotdot : StrokeType::Dashdot;
						GA.strokeWidth(e) = randomNumber(1, std::numeric_limits<int>::max());
					}
				}

				std::ostringstream write;
				auto flagsBefore = write.flags();
				AssertThat(writerGA(GA, write), Equals(true));
				AssertThat(write.flags(), Equals(flagsBefore));
				std::istringstream read(write.str());
				Graph G2;
				GraphAttributes GA2(G2, attr);
				flagsBefore = read.flags();
				AssertThat(readerGA(GA2, G2, read), Equals(true));
				AssertThat(read.flags(), Equals(flagsBefore));
				assertEqualGAs(GA, GA2);
			});
		});
	});
}

/**
 * Tests reading a DOT clustergraph, using a simplified version of:
 * https://graphviz.gitlab.io/_pages/Gallery/directed/cluster.html
 */
void describeDOTwithClusters()
{
	describe("DOT with subgraphs as clusters", []() {
		it("reads a cluster graph", []() {
			const ResourceFile* file = ResourceFile::get(
				"fileformats/dot/valid/cluster");
			std::stringstream is{file->data()};

			Graph G;
			ClusterGraph CG(G);

			const bool readStatus = GraphIO::readDOT(CG, G, is);
			AssertThat(readStatus, Equals(true));

			// this graph has two clusters inside the root cluster, each of which
			// has four nodes.
			AssertThat(CG.numberOfClusters(), Equals(3));
			AssertThat(CG.rootCluster()->children.size(), Equals(2));
			for (const auto &cluster : CG.rootCluster()->children) {
				AssertThat(cluster->nodes.size(), Equals(4));
			}
		});
	});
}

function<bool(Graph &G, istream &is)> toLambda(const GraphIO::ReaderFunc &reader) {
	return [&](Graph &G, istream &is) { return reader(G, is); };
}

go_bandit([]() {
describe("GraphIO", []() {
	describeSTP<int>("int");
	describeSTP<double>("double");
	describeSTPonlyGraph();

	describeDMF<int>("int");
	describeDMF<double>("double");

	describeGAFormat("GML", GraphIO::readGML, GraphIO::writeGML, GraphIO::readGML, GraphIO::writeGML, false,
	                 GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics | GraphAttributes::edgeDoubleWeight
	                 | GraphAttributes::edgeLabel | GraphAttributes::nodeLabel | GraphAttributes::edgeType
	                 | GraphAttributes::nodeId | GraphAttributes::edgeArrow | GraphAttributes::edgeStyle
	                 | GraphAttributes::nodeTemplate | GraphAttributes::nodeWeight | GraphAttributes::nodeStyle
	                 | GraphAttributes::threeD);
	describeFormat("Rome", toLambda(GraphIO::readRome), GraphIO::writeRome, false);
	describeFormat("LEDA", toLambda(GraphIO::readLEDA), GraphIO::writeLEDA, false);
	describeFormat("Chaco", toLambda(GraphIO::readChaco), GraphIO::writeChaco, false);
	describeFormat("PMDiss", toLambda(GraphIO::readPMDissGraph), GraphIO::writePMDissGraph, false);
	describeGAFormat("GraphML", GraphIO::readGraphML, GraphIO::writeGraphML, GraphIO::readGraphML, GraphIO::writeGraphML, true,
	                 GraphAttributes::nodeGraphics | GraphAttributes::edgeIntWeight | GraphAttributes::edgeDoubleWeight
	                 | GraphAttributes::edgeLabel | GraphAttributes::nodeLabel | GraphAttributes::nodeType
	                 | GraphAttributes::edgeType | GraphAttributes::edgeArrow | GraphAttributes::nodeTemplate
	                 | GraphAttributes::nodeWeight | GraphAttributes::threeD | GraphAttributes::nodeStyle);
	describeGAFormat("DOT", GraphIO::readDOT, GraphIO::writeDOT, GraphIO::readDOT, GraphIO::writeDOT, false,
	                 GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics | GraphAttributes::edgeIntWeight
	                 | GraphAttributes::edgeLabel | GraphAttributes::nodeLabel | GraphAttributes::nodeType
	                 | GraphAttributes::edgeDoubleWeight | GraphAttributes::edgeArrow | GraphAttributes::nodeTemplate
	                 | GraphAttributes::nodeWeight | GraphAttributes::nodeStyle | GraphAttributes::threeD);
	describeDOTwithClusters();
	describeGAFormat("GEXF", GraphIO::readGEXF, GraphIO::writeGEXF, GraphIO::readGEXF, GraphIO::writeGEXF, true,
	                 GraphAttributes::nodeGraphics | GraphAttributes::edgeIntWeight
	                 | GraphAttributes::edgeDoubleWeight | GraphAttributes::nodeType | GraphAttributes::edgeType
	                 | GraphAttributes::edgeArrow | GraphAttributes::nodeTemplate | GraphAttributes::nodeWeight
	                 | GraphAttributes::threeD | GraphAttributes::nodeStyle, true);
	describeGAFormat("GDF", GraphIO::readGDF, GraphIO::writeGDF, GraphIO::readGDF, GraphIO::writeGDF, false,
	                 GraphAttributes::nodeGraphics | GraphAttributes::edgeIntWeight | GraphAttributes::edgeGraphics
	                 | GraphAttributes::edgeLabel | GraphAttributes::nodeLabel
	                 | GraphAttributes::edgeDoubleWeight | GraphAttributes::nodeTemplate | GraphAttributes::nodeWeight
	                 | GraphAttributes::threeD | GraphAttributes::nodeStyle);
	describeGAFormat("DL", GraphIO::readDL, GraphIO::writeDL, GraphIO::readDL, GraphIO::writeDL, false,
	                 GraphAttributes::edgeIntWeight | GraphAttributes::edgeDoubleWeight | GraphAttributes::nodeLabel);
	describeGAFormat("TLP", GraphIO::readTLP, GraphIO::writeTLP, GraphIO::readTLP, GraphIO::writeTLP, false,
	                 GraphAttributes::nodeGraphics | GraphAttributes::edgeLabel | GraphAttributes::nodeLabel
	                 | GraphAttributes::threeD | GraphAttributes::nodeStyle);
	describeFormat("Graph6", [&](Graph& G, istream &is) { return GraphIO::readGraph6(G, is, false); },
	               GraphIO::writeGraph6, false);

	describe("generic reader", []() {
		std::function<void(const ResourceFile*)> genericTestTrue = [](const ResourceFile* file) {
			it("parses " + file->fullPath(), [&]() {
				Graph graph;
				stringstream ss{file->data()};
				AssertThat(GraphIO::read(graph, ss), IsTrue());
			});
		};

		std::function<void(const ResourceFile*)> genericTestFalse = [](const ResourceFile* file) {
			it("does not recognize " + file->fullPath(), [&]() {
				Graph graph;
				stringstream ss{file->data()};
				AssertThat(GraphIO::read(graph, ss), IsFalse());
			});
		};

		for_each_file("fileformats/gml/valid", genericTestTrue);
		for_each_file("fileformats/gml/invalid", genericTestFalse);

		for_each_file("fileformats/chaco/valid", genericTestTrue);
		for_each_file("fileformats/chaco/invalid", genericTestFalse);

		for_each_file("fileformats/dl/valid", genericTestTrue);
		for_each_file("fileformats/dl/invalid", genericTestFalse);

		for_each_file("fileformats/dot/valid", genericTestTrue);
		for_each_file("fileformats/dot/invalid", genericTestFalse);

		for_each_file("fileformats/gdf/valid", genericTestTrue);

		for_each_file("fileformats/gexf/valid", genericTestTrue);

		for_each_file("fileformats/graphml/valid", genericTestTrue);

		for_each_file("fileformats/leda/valid", genericTestTrue);
		for_each_file("fileformats/leda/invalid", genericTestFalse);

		for_each_file("fileformats/tlp/valid", genericTestTrue);
		for_each_file("fileformats/tlp/invalid", genericTestFalse);

		for_each_file("fileformats/stp/valid", genericTestTrue);

		for_each_file("fileformats/graph6/valid", genericTestTrue);

		for_each_file("fileformats/dmf/invalid", genericTestFalse);
	});
});
});
