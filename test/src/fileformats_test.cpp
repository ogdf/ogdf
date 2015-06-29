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

#include <algorithm>
#include <string>
#include <tuple>
#include <bandit/bandit.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/EpsilonTest.h>
#include <ogdf/basic/graph_generators.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/fileformats/xml/Parser.h>
#include <resources.h>

using namespace ogdf;
using namespace bandit;

typedef bool (*ReaderFunc)(Graph&, istream&);
typedef bool (*WriterFunc)(const Graph&, ostream&);

bool isSameUndirectedGraph(const Graph &G1, const Graph &G2)
{
	if (G1.numberOfNodes() != G2.numberOfNodes()
	 || G1.numberOfEdges() != G2.numberOfEdges()) {
		return false;
	}
	// we assert that node indices of G1 and G2 coincide
	for (node v1 = G1.firstNode(), v2 = G2.firstNode(); v1; v1 = v1->succ(), v2 = v2->succ()) {
		List<int> neighbors1, neighbors2;
		for (adjEntry adj1 = v1->firstAdj(); adj1; adj1 = adj1->succ()) {
			neighbors1.pushBack(adj1->twinNode()->index());
		}
		for (adjEntry adj2 = v2->firstAdj(); adj2; adj2 = adj2->succ()) {
			neighbors2.pushBack(adj2->twinNode()->index());
		}
		neighbors1.quicksort();
		neighbors2.quicksort();
		if (neighbors1 != neighbors2) {
			return false;
		}
	}
	return true;
}

template<typename T>
void describeSTP(const string &typeName) {
	describe(string("STP for " + typeName).c_str(), [](){
		for(int i = 4; i < 1024; i *= 2) {
			it(string("stores and loads an instance of size " + to_string(i)).c_str(), [&](){
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

		for_each_file("fileformats/stp/valid", [&](const string &filename){
			it(string("successfully parses " + filename).c_str(), [&](){
				EdgeWeightedGraph<T> graph;
				List<node> terminals;
				NodeArray<bool> isTerminal;
				std::ifstream input(filename);
				AssertThat(GraphIO::readSTP(graph, terminals, isTerminal, input), IsTrue());

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

		for_each_file("fileformats/stp/invalid", [&](const string &filename){
			it(string("detects errors in " + filename).c_str(), [&](){
				EdgeWeightedGraph<T> graph;
				List<node> terminals;
				NodeArray<bool> isTerminal;
				std::ifstream input(filename);
				AssertThat(GraphIO::readSTP(graph, terminals, isTerminal, input), IsFalse());
			});
		});
	});
}

template<typename T>
void describeDMF(const string &typeName) {
	describe(string("DMF for " + typeName).c_str(), [](){
		const void* nullPointer = nullptr;

		for_each_file("fileformats/dmf/valid", [&](const string &filename){
			it(string("reads " + filename).c_str(), [&](){
				EdgeWeightedGraph<T> graph;
				node source;
				node sink;
				std::ifstream input(filename);

				AssertThat(GraphIO::readDMF(graph, source, sink, input), IsTrue());
				AssertThat(graph.numberOfNodes(), IsGreaterThan(1));
				AssertThat(source, Is().Not().EqualTo(nullPointer));
				AssertThat(sink, Is().Not().EqualTo(nullPointer));
#ifdef OGDF_DEBUG
				AssertThat(source->graphOf(), Equals(&graph));
				AssertThat(sink->graphOf(), Equals(&graph));
#endif
				AssertThat(source, Is().Not().EqualTo(sink));

				for(edge e : graph.edges) {
					AssertThat(graph.weight(e) > 0, IsTrue());
				}
			});
		});

		for_each_file("fileformats/dmf/invalid", [&](const string &filename){
			it(string("reads " + filename).c_str(), [&](){
				EdgeWeightedGraph<T> graph;
				node source;
				node sink;
				std::ifstream input(filename);

				AssertThat(GraphIO::readDMF(graph, source, sink, input), IsFalse());
			});
		});

		it("writes and reads a random graph", [&](){
			EdgeWeightedGraph<T> graph;
			node source;
			node sink;

			randomGraph(graph, 42, 189);
			sink = source = graph.chooseNode();
			while(sink == source) {
				sink = graph.chooseNode();
			}

			T sum = 0;
			for(edge e : graph.edges) {
				T cap = randomDoubleNormal(10, 5);
				if(cap < 0) {
					cap *= -1;
				}
				graph.setWeight(e, cap);
				sum += cap;
			}

			std::ostringstream writeStream;

			AssertThat(GraphIO::writeDMF(graph, source, sink, writeStream), IsTrue());

			EdgeWeightedGraph<T> readGraph;
			node readSource = nullptr;
			node readSink = nullptr;

			std::istringstream readStream(writeStream.str());
			AssertThat(GraphIO::readDMF(readGraph, readSource, readSink, readStream), IsTrue());

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
				readSum += readGraph.weight(e);
			}

			EpsilonTest eps(1.0e-3);
			AssertThat(eps.equal(sum, readSum), IsTrue());
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
void describeFormat(const std::string name, ReaderFunc reader, WriterFunc writer, bool isXml)
{
	std::string lowerCaseName = name;
	std::transform(lowerCaseName.begin(), lowerCaseName.end(), lowerCaseName.begin(), ::tolower);

	auto errorTest = [&](const string &filename){
		it(string("detects errors in " + filename).c_str(), [&](){
			Graph graph;
			std::ifstream input(filename);
			AssertThat(reader(graph, input), IsFalse());
		});
	};

	auto resourceBasedTest = [&](){
		for_each_file("fileformats/" + lowerCaseName + "/valid", [&](const string &filename){
			it(string("successfully parses " + filename).c_str(), [&](){
				Graph graph;
				std::ifstream input(filename);
				AssertThat(reader(graph, input), IsTrue());
				AssertThat(graph.numberOfNodes(), IsGreaterThan(0));
				AssertThat(graph.numberOfEdges(), IsGreaterThan(0));
			});
		});
		for_each_file("fileformats/" + lowerCaseName + "/valid/skip", [&](const string &filename) {
			it_skip(string("successfully parses " + filename).c_str(), [&]() {});
		});

		for_each_file("fileformats/" + lowerCaseName + "/invalid", errorTest);
		for_each_file("fileformats/" + lowerCaseName + "/invalid/skip", [&](const string &filename) {
			it_skip(string("detects errors in " + filename).c_str(), [&]() {});
		});
	};

	describe(name.c_str(), [&](){
		if(isXml) {
			for_each_file("fileformats/xml/invalid", errorTest);
		}

		resourceBasedTest();

		it("writes and reads an empty graph", [&](){
			Graph G, Gtest;
			std::ostringstream write;
			AssertThat(writer(G, write), Equals(true));
			std::istringstream read(write.str());
			AssertThat(reader(Gtest, read), Equals(true));
			AssertThat(isSameUndirectedGraph(G, Gtest), Equals(true));
		});

		it("writes and reads a graph of isolated nodes", [&](){
			Graph G, Gtest;
			G.newNode(); G.newNode();
			std::ostringstream write;
			AssertThat(writer(G, write), Equals(true));
			std::istringstream read(write.str());
			AssertThat(reader(Gtest, read), Equals(true));
			AssertThat(isSameUndirectedGraph(G, Gtest), Equals(true));
		});

		// TODO: Skipping GDF because it fails
		if(!name.compare("GDF")) {
			return;
		}

		it("writes and reads a Petersen graph", [&](){
			Graph G, Gtest;
			petersenGraph(G, 5, 2);
			std::ostringstream write;
			AssertThat(writer(G, write), Equals(true));
			std::istringstream read(write.str());
			AssertThat(reader(Gtest, read), Equals(true));
			AssertThat(isSameUndirectedGraph(G, Gtest), Equals(true));
		});

		// TODO: Skipping OGML because it fails
		if (!name.compare("OGML")) {
			return;
		}

// TODO Skipping big complete graph on windows since the test crashes
#ifdef WIN32
		it_skip(
#else
		it(
#endif
		"writes and reads a big complete graph", [&](){
			Graph G, Gtest;
			completeGraph(G, 243);
			std::ostringstream write;
			AssertThat(writer(G, write), Equals(true));
			std::istringstream read(write.str());
			AssertThat(reader(Gtest, read), Equals(true));
			AssertThat(isSameUndirectedGraph(G, Gtest), Equals(true));
		});
	});
}

go_bandit([](){
describe("GraphIO", [](){
	describeSTP<int>("int");
	describeSTP<double>("double");

	describeDMF<int>("int");
	describeDMF<double>("double");

	describe("XML", [&](){
		for_each_file("fileformats/xml/valid", [&](const string &filename){
			it(string("successfully parses " + filename).c_str(), [&](){
				xml::Parser xmlParser;
				std::ifstream input(filename);
				xmlParser.parse(input);
			});
		});
			for_each_file("fileformats/xml/invalid", [&](const string &filename){
				it(string("detects errors in " + filename).c_str(), [&](){
					xml::Parser xmlParser;
					std::ifstream input(filename);
					AssertThrows(xml::ParseException, xmlParser.parse(input));
				});
			});
	});

	describeFormat("GML", GraphIO::readGML, GraphIO::writeGML, false);
	describeFormat("OGML", GraphIO::readOGML, GraphIO::writeOGML, true);
	describeFormat("Rome", GraphIO::readRome, GraphIO::writeRome, false);
	describeFormat("LEDA", GraphIO::readLEDA, GraphIO::writeLEDA, false);
	describeFormat("Chaco", GraphIO::readChaco, GraphIO::writeChaco, false);
	describeFormat("PMDiss", GraphIO::readPMDissGraph, GraphIO::writePMDissGraph, false);
	describeFormat("GraphML", GraphIO::readGraphML, GraphIO::writeGraphML, true);
	describeFormat("DOT", GraphIO::readDOT, GraphIO::writeDOT, false);
	describeFormat("GEXF", GraphIO::readGEXF, GraphIO::writeGEXF, true);
	describeFormat("GDF", GraphIO::readGDF, GraphIO::writeGDF, false);
	describeFormat("TLP", GraphIO::readTLP, GraphIO::writeTLP, false);
	describeFormat("DL", GraphIO::readDL, GraphIO::writeDL, false);
	describeFormat("Graph6", GraphIO::readGraph6, GraphIO::writeGraph6, false);
});
});
