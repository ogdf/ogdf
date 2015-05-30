/*
 * $Revision$
 *
 * last checkin:
 *   $Author$
 *   $Date$
 ***************************************************************/

/** \file
 * \brief Tests for fileformat reading and writing using GraphIO,
 *   only graphs without attributes
 *
 * \author Stephan Beyer
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

#include "gtest/gtest.h"
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/graph_generators.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/fileformats/GraphIO.h>

// Note: these tests do not do real file testing,
// all file IO is simulated over a stringstream.

using namespace ogdf;

bool
isSameUndirectedGraph(const Graph &G1, const Graph &G2)
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

TEST(FileformatsTest, OgmlSimpleRead)
{
	std::stringstream ss;
	ss <<
	  "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
	  "<ogml>\n"
	  "  <!-- a simple comment with </ogml> and <tag>\n"
	  "       over two lines -->\n"
	  "  <graph>\n"
	  "    <structure>\n"
	  "      <node id=\"v1\" />\n"
	  "      <node id=\"v2\"></node>\n"
	  "      <edge id=\"e1\">\n"
	  "        <source idRef=\"v1\" />\n"
	  "        <target idRef=\"v2\"></target>\n"
	  "      </edge>\n"
	  "    </structure>\n"
	  "\n"
	  "    <layout>\n"
	  "      <styleTemplates>\n"
	  "      </styleTemplates>\n"
	  "      <styles/>\n"
	  "    </layout>\n"
	  "  </graph>\n"
	  "</ogml>\n";

	Graph G;
	EXPECT_TRUE(GraphIO::readOGML(G, ss));
	EXPECT_EQ(2, G.numberOfNodes());
	EXPECT_EQ(1, G.numberOfEdges());
}

TEST(FileformatsTest, OgmlReadFailMissingCompulsiveTags)
{
	std::stringstream ss;
	ss << "<ogml></ogml>";
	Graph G;
	EXPECT_FALSE(GraphIO::readOGML(G, ss));
}

TEST(FileformatsTest, OgmlReadFailNoXml)
{
	std::stringstream ss;
	ss << "no XML at all";
	Graph G;
	EXPECT_FALSE(GraphIO::readOGML(G, ss));
}

TEST(FileformatsTest, OgmlReadFailInvalidXml)
{
	std::stringstream ss;
	ss <<
	  "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
	  "<ogml>\n"
	  "  <graph>\n"
	  "    <structure>\n"
	  "      <node id=\"v1\" />\n"
	  "      <node id=\"v2\"></node>\n"
	  "      <edge id=\"e1\">\n"
	  "        <source idRef=\"v1\" />\n"
	  "        <target idRef=\"v2\"></target>\n"
	  "      </edge>\n"
	  "    </structure>\n";
	Graph G;
	EXPECT_FALSE(GraphIO::readOGML(G, ss));
}

TEST(FileformatsTest, OgmlReadWriteEmptyGraph)
{
	Graph G, Gtest;
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeOGML(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readOGML(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, OgmlReadWriteIsolatedNodes)
{
	Graph G, Gtest;
	G.newNode(); G.newNode();
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeOGML(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readOGML(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, OgmlReadWritePetersenGraph)
{
	Graph G, Gtest;
	petersenGraph(G, 5, 2);
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeOGML(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readOGML(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, GmlReadWriteEmptyGraph)
{
	Graph G, Gtest;
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeGML(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readGML(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, GmlReadWriteIsolatedNodes)
{
	Graph G, Gtest;
	G.newNode(); G.newNode();
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeGML(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readGML(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, GmlReadWritePetersenGraph)
{
	Graph G, Gtest;
	petersenGraph(G, 5, 2);
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeGML(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readGML(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, RomeReadWriteEmptyGraph)
{
	Graph G, Gtest;
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeRome(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readRome(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, RomeReadWriteIsolatedNodes)
{
	Graph G, Gtest;
	G.newNode(); G.newNode();
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeRome(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readRome(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, RomeReadWritePetersenGraph)
{
	Graph G, Gtest;
	petersenGraph(G, 5, 2);
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeRome(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readRome(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, LedaReadWriteEmptyGraph)
{
	Graph G, Gtest;
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeLEDA(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readLEDA(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, LedaReadWriteIsolatedNodes)
{
	Graph G, Gtest;
	G.newNode(); G.newNode();
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeLEDA(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readLEDA(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, LedaReadWritePetersenGraph)
{
	Graph G, Gtest;
	petersenGraph(G, 5, 2);
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeLEDA(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readLEDA(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, ChacoReadFailNodeIndex)
{
	Graph G;
	std::stringstream ss;
	ss << "10 3\n 13 1\n";
	EXPECT_FALSE(GraphIO::readChaco(G, ss));
}

TEST(FileformatsTest, ChacoReadWriteEmptyGraph)
{
	Graph G, Gtest;
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeChaco(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readChaco(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, ChacoReadWriteIsolatedNodes)
{
	Graph G, Gtest;
	G.newNode(); G.newNode();
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeChaco(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readChaco(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, ChacoReadWritePetersenGraph)
{
	Graph G, Gtest;
	petersenGraph(G, 5, 2);
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeChaco(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readChaco(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, PMDissReadWriteEmptyGraph)
{
	Graph G, Gtest;
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writePMDissGraph(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readPMDissGraph(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, PMDissReadWriteIsolatedNodes)
{
	Graph G, Gtest;
	G.newNode(); G.newNode();
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writePMDissGraph(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readPMDissGraph(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, PMDissReadWritePetersenGraph)
{
	Graph G, Gtest;
	petersenGraph(G, 5, 2);
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writePMDissGraph(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readPMDissGraph(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, GraphMLReadFailMissingCompulsiveTags)
{
	std::stringstream ss;
	ss << "<graphml></graphml>";
	Graph G;
	EXPECT_FALSE(GraphIO::readGraphML(G, ss));
}

TEST(FileformatsTest, GraphMLReadFailNoXml)
{
	std::stringstream ss;
	ss << "no XML at all";
	Graph G;
	EXPECT_FALSE(GraphIO::readGraphML(G, ss));
}

TEST(FileformatsTest, GraphMLReadFailInvalidXml)
{
	std::stringstream ss;
	ss <<
	  "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
	  "<graphml>\n"
	  "  <graph>\n"
	  "</graphml>\n";
	Graph G;
	EXPECT_FALSE(GraphIO::readGraphML(G, ss));
}

TEST(FileformatsTest, GraphMLReadWriteEmptyGraph)
{
	Graph G, Gtest;
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeGraphML(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readGraphML(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, GraphMLReadWriteIsolatedNodes)
{
	Graph G, Gtest;
	G.newNode(); G.newNode();
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeGraphML(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readGraphML(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, GraphMLReadWritePetersenGraph)
{
	Graph G, Gtest;
	petersenGraph(G, 5, 2);
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeGraphML(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readGraphML(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, DotReadWriteEmptyGraph)
{
	Graph G, Gtest;
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeDOT(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readDOT(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, DotReadWriteIsolatedNodes)
{
	Graph G, Gtest;
	G.newNode(); G.newNode();
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeDOT(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readDOT(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, DotReadWritePetersenGraph)
{
	Graph G, Gtest;
	petersenGraph(G, 5, 2);
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeDOT(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readDOT(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, GephiReadFailNoXml)
{
	std::stringstream ss;
	ss << "no XML at all";
	Graph G;
	EXPECT_FALSE(GraphIO::readGEXF(G, ss));
}

TEST(FileformatsTest, GephiReadWriteEmptyGraph)
{
	Graph G, Gtest;
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeGEXF(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readGEXF(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, GephiReadWriteIsolatedNodes)
{
	Graph G, Gtest;
	G.newNode(); G.newNode();
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeGEXF(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readGEXF(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, GephiReadWritePetersenGraph)
{
	Graph G, Gtest;
	petersenGraph(G, 5, 2);
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeGEXF(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readGEXF(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, GdfReadFailCorruptedHeader)
{
	std::stringstream ss;
	ss << "nonedef> should,not,work\n";
	Graph G;
	EXPECT_FALSE(GraphIO::readGDF(G, ss));
}

TEST(FileformatsTest, GdfReadWriteEmptyGraph)
{
	Graph G, Gtest;
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeGDF(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readGDF(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, GdfReadWriteIsolatedNodes)
{
	Graph G, Gtest;
	G.newNode(); G.newNode();
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeGDF(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readGDF(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, GdfReadWritePetersenGraph)
{
	Graph G, Gtest;
	petersenGraph(G, 5, 2);
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeGDF(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readGDF(Gtest, read));
	// GDF file specification says, an undirected edge becomes bidirected
	EXPECT_EQ(2 * G.numberOfEdges(), Gtest.numberOfEdges());
	makeParallelFreeUndirected(Gtest);
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, TulipReadWriteEmptyGraph)
{
	Graph G, Gtest;
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeTLP(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readTLP(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, TulipReadWriteIsolatedNodes)
{
	Graph G, Gtest;
	G.newNode(); G.newNode();
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeTLP(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readTLP(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, TulipReadWritePetersenGraph)
{
	Graph G, Gtest;
	petersenGraph(G, 5, 2);
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeTLP(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readTLP(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, DLReadWriteEmptyGraph)
{
	Graph G, Gtest;
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeDL(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readDL(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, DLReadWriteIsolatedNodes)
{
	Graph G, Gtest;
	G.newNode(); G.newNode();
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeDL(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readDL(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}

TEST(FileformatsTest, DLReadWritePetersenGraph)
{
	Graph G, Gtest;
	petersenGraph(G, 5, 2);
	std::ostringstream write;
	ASSERT_TRUE(GraphIO::writeDL(G, write));
	std::istringstream read(write.str());
	ASSERT_TRUE(GraphIO::readDL(Gtest, read));
	EXPECT_TRUE(isSameUndirectedGraph(G, Gtest));
}
