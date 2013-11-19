/*
 * $Revision: 3831 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 10:00:32 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Implements class GraphIO which provides access to all
 *        graph read and write functionality.
 *
 * \author Carsten Gutwenger, Markus Chimani, Karsten Klein
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

#include <ogdf/basic/Logger.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/fileformats/GmlParser.h>
#include <ogdf/fileformats/OgmlParser.h>
#include <ogdf/fileformats/GraphMLParser.h>
#include <ogdf/fileformats/DotParser.h>
#include <ogdf/fileformats/GexfParser.h>
#include <ogdf/fileformats/GdfParser.h>
#include <ogdf/fileformats/TlpParser.h>
#include <ogdf/fileformats/DLParser.h>
#include <sstream>
#include <map>


// we use these data structures from the stdlib
using std::map;
using std::istringstream;


namespace ogdf {


char GraphIO::s_indentChar  = ' ';
int  GraphIO::s_indentWidth = 2;


ostream &GraphIO::indent(ostream &os, int depth)
{
	int n = s_indentWidth * depth;
	for( ; n > 0; --n)
		os.put(s_indentChar);

	return os;
}


//---------------------------------------------------------
// Graph: GML format
//---------------------------------------------------------

bool GraphIO::readGML(Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) return false;
	return readGML(G, is);
}

bool GraphIO::readGML(Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) return false;
	return readGML(G, is);
}

bool GraphIO::readGML(Graph &G, istream &is)
{
	GmlParser parser(is);
	if (parser.error()) return false;
	return parser.read(G);
}


bool GraphIO::writeGML(const Graph &G, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeGML(G, os);
}

bool GraphIO::writeGML(const Graph &G, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeGML(G, os);
}


//---------------------------------------------------------
// Graph: OGML format
//---------------------------------------------------------

bool GraphIO::readOGML(Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) return false;
	return readOGML(G, is);
}

bool GraphIO::readOGML(Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) return false;
	return readOGML(G, is);
}

bool GraphIO::readOGML(Graph &G, istream &is)
{
	OgmlParser parser;
	return parser.read(is, G);
}


bool GraphIO::writeOGML(const Graph &G, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeOGML(G, os);
}

bool GraphIO::writeOGML(const Graph &G, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeOGML(G, os);
}


//---------------------------------------------------------
// Graph: Rome format
//---------------------------------------------------------

bool GraphIO::readRome(Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) return false;
	return readRome(G, is);
}

bool GraphIO::readRome(Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) return false;
	return readRome(G, is);
}


bool GraphIO::readRome(Graph &G, istream &is)
{
	G.clear();  // start with empty graph

	bool readNodes = true;
	map<int,node> indexToNode;

	string buffer;
	istringstream iss;
	while(std::getline(is, buffer))
	{
		if(buffer.size() == 0)
			continue;

		iss.str(buffer);
		iss.clear();

		if(readNodes) {
			if(buffer[0] == '#') {
				readNodes = false;
				continue;
			}

			int index = -1;
			iss >> index;
			if(index < 1 || indexToNode.find(index) != indexToNode.end()) {
				Logger::slout() << "GraphIO::readRome: Illegal node index!\n";
				return false;
			}

			indexToNode[index] = G.newNode();

		} else {

			int index, dummy, srcIndex = -1, tgtIndex = -1;
			iss >> index >> dummy >> srcIndex >> tgtIndex;

			map<int,node>::const_iterator itSrc = indexToNode.find(srcIndex);
			map<int,node>::const_iterator itTgt = indexToNode.find(tgtIndex);

			if(itSrc == indexToNode.end() || itTgt == indexToNode.end()) {
				Logger::slout() << "GraphIO::readRome: Illegal node index in edge specification.\n";
				return false;
			}

			G.newEdge(itSrc->second, itTgt->second);
		}
	}
	return true;
}


bool GraphIO::writeRome(const Graph &G, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeRome(G, os);
}

bool GraphIO::writeRome(const Graph &G, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeRome(G, os);
}


bool GraphIO::writeRome(const Graph &G, ostream &os)
{
	// assign indices 1, 2, 3, ... to nodes
	NodeArray<int> index(G);

	int i = 0;
	node v;
	forall_nodes(v,G) {
		index[v] = ++i;
		// write node v
		os << i << " " << "0\n";
	}

	os << "#\n"; // write node-edge separator

	i = 0;
	edge e;
	forall_edges(e,G) {
		// write edge e
		os << ++i << " 0 " << index[e->source()] << " " << index[e->target()] << "\n";
	}

	return true;
}


//---------------------------------------------------------
// Graph: LEDA format
//---------------------------------------------------------

bool GraphIO::readLEDA(Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) return false;
	return readLEDA(G, is);
}

bool GraphIO::readLEDA(Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) return false;
	return readLEDA(G, is);
}


bool GraphIO::writeLEDA(const Graph &G, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeLEDA(G, os);
}

bool GraphIO::writeLEDA(const Graph &G, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeLEDA(G, os);
}


//---------------------------------------------------------
// Graph: Chaco format
//---------------------------------------------------------

bool GraphIO::readChaco(Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) return false;
	return readChaco(G, is);
}

bool GraphIO::readChaco(Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) return false;
	return readChaco(G, is);
}

bool GraphIO::readChaco(Graph &G, istream &is)
{
	G.clear();

	string buffer;
	istringstream iss;

	int numN = -1, numE = -1;
	if (std::getline(is, buffer)) {
		iss.str(buffer);
		iss >> numN >> numE;
		if(numN < 0 || numE < 0)
			return false;
	} else
		return false;

	if (numN == 0) return true;

	Array<node> indexToNode(1,numN,0);
	for (int i = 1; i <= numN; i++)
		indexToNode[i] = G.newNode();

	int vid = 0;
	while(std::getline(is, buffer))
	{
		if(buffer.empty())
			continue;

		if(vid > numN) {
			Logger::slout() << "GraphIO::readChaco: More lines with adjacency lists than expected.\n";
			return false;
		}

		iss.str(buffer); iss.clear();
		node v = indexToNode[++vid];

		int wid;
		while(iss >> wid) {
			if(wid < 1 || wid > numN) {
				Logger::slout() << "GraphIO::readChaco: Illegal node index in adjacency list.\n";
				return false;
			}
			if(wid >= vid)
				G.newEdge(v, indexToNode[wid]);
		}
	}

	return true;
}


bool GraphIO::writeChaco(const Graph &G, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeChaco(G, os);
}

bool GraphIO::writeChaco(const Graph &G, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeChaco(G, os);
}

bool GraphIO::writeChaco(const Graph &G, ostream &os)
{
	os << G.numberOfNodes() << " " << G.numberOfEdges() << "\n";

	NodeArray<int> index(G);

	node v;
	int count = 0;
	forall_nodes(v,G)
		index[v] = ++count;

	forall_nodes(v,G) {
		adjEntry adj;
		forall_adj(adj,v)
			os << " " << index[adj->twinNode()];
		os << "\n";
	}

	return true;
}


//---------------------------------------------------------
// Graph: Chaco format
//---------------------------------------------------------

bool GraphIO::readYGraph(Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) return false;
	return readYGraph(G, is);
}

bool GraphIO::readYGraph(Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) return false;
	return readYGraph(G, is);
}

bool GraphIO::readYGraph(Graph &G, istream &is)
{
	const char *errorLineTooShort = "GraphIO::readYGraph: line too short!\n";

	G.clear();

	if(!is) {
		Logger::slout() << errorLineTooShort;
		return false;
	}

	int n = is.get();
	if(!is.good() || n == '\n' || n < 0) {
		Logger::slout() << errorLineTooShort;
		return false;
	}
	n &= 0x3F;

	Array<node> indexToNode(n);
	for(int i = n; i-- > 0; )
		indexToNode[i] = G.newNode();

	int s = 0, c;
	for(int i = 1; i < n; ++i)
	{
		for(int j = 0; j < i; ++j) {
			if(!s) {
				c = is.get();
				if(!is.good() || c == '\n') {
					Logger::slout() << errorLineTooShort;
					return false;
				}
				c &= 0x3F;

				s = 5;
			} else --s;
			if(c & (1 << s))
				G.newEdge(indexToNode[i], indexToNode[j]);
		}
	}

	c = is.get();
	if(!is.eof() && c != '\n') {
		Logger::slout(Logger::LL_MINOR) << "GraphIO::readYGraph: Warning: line too long! ignoring...";
	}

	return true;
}


//---------------------------------------------------------
// Graph: Petra Mutzel Diss format
//---------------------------------------------------------

bool GraphIO::readPMDissGraph(Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) return false;
	return readPMDissGraph(G, is);
}

bool GraphIO::readPMDissGraph(Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) return false;
	return readPMDissGraph(G, is);
}

bool GraphIO::readPMDissGraph(Graph &G, istream &is)
{
	const char *errorInFileHeader = "GraphIO::readPMDissGraph: Error in file header.\n";

	G.clear();

	string buffer;
	istringstream iss;

	int numN = -1, numE = -1;

	// first two lines look as follows (example with 20 nodes, 30 edges):
	// *BEGIN unknown_comp.20.30
	// *GRAPH 20 30 UNDIRECTED UNWEIGHTED

	if (std::getline(is, buffer))
	{
		iss.str(buffer); iss.clear();

		string str;
		iss >> str;
		if(str != "*BEGIN") {
			Logger::slout() << "GraphIO::readPMDissGraph: Error in file header, could not find \"*BEGIN\".\n";
			return false;
		}

		if (std::getline(is, buffer)) {
			iss.str(buffer); iss.clear();

			iss >> str >> numN >> numE;

			if(str != "*GRAPH" || numN < 0 || numE < 0) {
				Logger::slout() << errorInFileHeader;
				return false;
			}
		}
		else {
			Logger::slout() << errorInFileHeader;
			return false;
		}
	}
	else {
		Logger::slout() << errorInFileHeader;
		return false;
	}

	if (numN == 0)
		return true;

	Array<node> indexToNode(1,numN,0);
	for (int i = 1; i <= numN; i++)
	{
		indexToNode[i] = G.newNode();
	}

	while(std::getline(is, buffer))
	{
		if(buffer.empty())
			continue;

		if(buffer[0] == '*')
			continue;

		iss.str(buffer); iss.clear();

		int srcIndex = -1, tgtIndex = -1;
		iss >> srcIndex >> tgtIndex;

		if(srcIndex < 1 || srcIndex > numN || tgtIndex < 1 || tgtIndex > numN) {
			Logger::slout() << "GraphIO::readPMDissGraph: Illegal node index in edge specification.\n";
			return false;
		}

		G.newEdge(indexToNode[srcIndex], indexToNode[tgtIndex]);
	}
	return true;
}


bool GraphIO::writePMDissGraph(const Graph &G, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writePMDissGraph(G, os);
}

bool GraphIO::writePMDissGraph(const Graph &G, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writePMDissGraph(G, os);
}

bool GraphIO::writePMDissGraph(const Graph &G, ostream &os)
{
	os << "*BEGIN unknown_name." << G.numberOfNodes() << "." << G.numberOfEdges() << "\n";
	os << "*GRAPH " << G.numberOfNodes() << " " << G.numberOfEdges() << " UNDIRECTED UNWEIGHTED\n";

	NodeArray<int> index(G);
	int nextIndex = 1;
	node v;
	forall_nodes(v,G)
		index[v] = nextIndex++;

	edge e;
	forall_edges(e,G)
		os << index[e->source()] << " " << index[e->target()] << "\n";

	os << "*CHECKSUM -1\n";
	os << "*END unknown_name." << G.numberOfNodes() << "." << G.numberOfEdges() << "\n";

	return true;
}


//---------------------------------------------------------
// ClusterGraph: GML format
//---------------------------------------------------------

bool GraphIO::readGML(ClusterGraph &C, Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) return false;
	return readGML(C, G, is);
}

bool GraphIO::readGML(ClusterGraph &C, Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) return false;
	return readGML(C, G, is);
}

bool GraphIO::readGML(ClusterGraph &C, Graph &G, istream &is)
{
	GmlParser gml(is);
	if (gml.error())
		return false;

	if(!gml.read(G))
		return false;

	return gml.readCluster(G, C);
}

bool GraphIO::writeGML(const ClusterGraph &C, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeGML(C, os);
}

bool GraphIO::writeGML(const ClusterGraph &C, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeGML(C, os);
}


//---------------------------------------------------------
// ClusterGraph: OGML format
//---------------------------------------------------------

bool GraphIO::readOGML(ClusterGraph &C, Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) return false;
	return readOGML(C, G, is);
}

bool GraphIO::readOGML(ClusterGraph &C, Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) return false;
	return readOGML(C, G, is);
}

bool GraphIO::readOGML(ClusterGraph &C, Graph &G, istream &is)
{
	OgmlParser parser;
	return parser.read(is, G, C);
}

bool GraphIO::writeOGML(const ClusterGraph &C, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeOGML(C, os);
}

bool GraphIO::writeOGML(const ClusterGraph &C, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeOGML(C, os);
}


//---------------------------------------------------------
// GraphAttributes: GML format
//---------------------------------------------------------

bool GraphIO::readGML(GraphAttributes &A, Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) return false;
	return readGML(A, G, is);
}

bool GraphIO::readGML(GraphAttributes &A, Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) return false;
	return readGML(A, G, is);
}

bool GraphIO::readGML(GraphAttributes &A, Graph &G, istream &is)
{
	GmlParser parser(is);
	if (parser.error()) return false;
	return parser.read(G, A);
}


bool GraphIO::writeGML(const GraphAttributes &A, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeGML(A, os);
}

bool GraphIO::writeGML(const GraphAttributes &A, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeGML(A, os);
}


//---------------------------------------------------------
// GraphAttributes: OGML format
//---------------------------------------------------------

bool GraphIO::readOGML(GraphAttributes &A, Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) return false;
	return readOGML(A, G, is);
}

bool GraphIO::readOGML(GraphAttributes &A, Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) return false;
	return readOGML(A, G, is);
}

bool GraphIO::readOGML(GraphAttributes &A, Graph &G, istream &is)
{
	OgmlParser parser;
	return parser.read(is, G, A);
}


bool GraphIO::writeOGML(const GraphAttributes &A, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeOGML(A, os);
}

bool GraphIO::writeOGML(const GraphAttributes &A, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeOGML(A, os);
}


//---------------------------------------------------------
// GraphAttributes: Rudy format
//---------------------------------------------------------

bool GraphIO::readRudy(GraphAttributes &A, Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) return false;
	return readRudy(A, G, is);
}

bool GraphIO::readRudy(GraphAttributes &A, Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) return false;
	return readRudy(A, G, is);
}

bool GraphIO::readRudy(GraphAttributes &A, Graph &G, istream &is)
{
	if(!is) return false;

	G.clear();

	int n, m;
	is >> n >> m;

	if(n < 0 || m < 0) {
		Logger::slout() << "GraphIO::readRudy: Illegal number of nodes or edges!\n";
		return false;
	}

	Array<node> mapToNode(0,n-1,0);
	for(int i = 0; i < n; ++i)
		mapToNode[i] = G.newNode();

	bool haveDoubleWeight = (A.attributes() & GraphAttributes::edgeDoubleWeight) != 0;

	for(int i = 0; i < m; i++)
	{
		int src = 0, tgt = 0;
		double weight = 1.0;

		is >> src >> tgt >> weight;
		if(src < 1 || src > n || tgt < 1 || tgt > n) {
			Logger::slout() << "GraphIO::readRudy: Illegal node index!\n";
			return false;
		}

		src--; tgt--;

		edge e = G.newEdge(mapToNode[src],mapToNode[tgt]);
		if (haveDoubleWeight)
			A.doubleWeight(e) = weight;
	}

	return true;
}


bool GraphIO::writeRudy(const GraphAttributes &A, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeRudy(A, os);
}

bool GraphIO::writeRudy(const GraphAttributes &A, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeRudy(A, os);
}

bool GraphIO::writeRudy(const GraphAttributes &A, ostream &os)
{
	const Graph &G = A.constGraph();
	os << G.numberOfNodes() << " " << G.numberOfEdges() << endl;

	// assign indices 1, 2, 3, ... to nodes
	NodeArray<int> index(G);

	int i = 0;
	node v;
	forall_nodes(v,G)
		index[v] = ++i;

	bool haveDoubleWeight = (A.attributes() & GraphAttributes::edgeDoubleWeight) != 0;

	edge e;
	forall_edges(e,G) {
		double w = (haveDoubleWeight) ? A.doubleWeight(e) : 1.0;
		os << index[e->source()] << " " << index[e->target()] << " " << w << "\n";
	}

	return true;
}


//---------------------------------------------------------
// ClusterGraphAttributes: GML format
//---------------------------------------------------------

bool GraphIO::readGML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) return false;
	return readGML(A, C, G, is);
}

bool GraphIO::readGML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) return false;
	return readGML(A, C, G, is);
}

bool GraphIO::readGML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, istream &is)
{
	GmlParser gml(is);
	if (gml.error())
		return false;

	if(!gml.read(G, A))
		return false;

	return gml.readAttributedCluster(G, C, A);
}


bool GraphIO::writeGML(const ClusterGraphAttributes &A, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeGML(A, os);
}

bool GraphIO::writeGML(const ClusterGraphAttributes &A, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeGML(A, os);
}



//---------------------------------------------------------
// ClusterGraphAttributes: OGML format
//---------------------------------------------------------

bool GraphIO::readOGML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) return false;
	return readOGML(A, C, G, is);
}

bool GraphIO::readOGML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) return false;
	return readOGML(A, C, G, is);
}

bool GraphIO::readOGML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, istream &is)
{
	OgmlParser parser;
	return parser.read(is, G, C, A);
}


bool GraphIO::writeOGML(const ClusterGraphAttributes &A, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeOGML(A, os);
}

bool GraphIO::writeOGML(const ClusterGraphAttributes &A, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeOGML(A, os);
}


//---------------------------------------------------------
// Hypergraphs (point-based expansion): PLA format
//---------------------------------------------------------

bool GraphIO::readBENCH(Graph &G, List<node>& hypernodes, List<edge>* shell, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) return false;
	return readBENCH(G, hypernodes, shell, is);
}

bool GraphIO::readBENCH(Graph &G, List<node>& hypernodes, List<edge>* shell, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) return false;
	return readBENCH(G, hypernodes, shell, is);
}


bool GraphIO::readPLA(Graph &G, List<node>& hypernodes, List<edge>* shell, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) return false;
	return readPLA(G, hypernodes, shell, is);
}

bool GraphIO::readPLA(Graph &G, List<node>& hypernodes, List<edge>* shell, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) return false;
	return readPLA(G, hypernodes, shell, is);
}


//---------------------------------------------------------
// Graph Drawing Challenge
//---------------------------------------------------------

bool GraphIO::readChallengeGraph(Graph &G, GridLayout &gl, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) return false;
	return readChallengeGraph(G, gl, is);
}

bool GraphIO::readChallengeGraph(Graph &G, GridLayout &gl, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) return false;
	return readChallengeGraph(G, gl, is);
}

bool GraphIO::readChallengeGraph(Graph &G, GridLayout &gl, istream &is)
{
	G.clear();

	string buffer;
	istringstream iss;

	int n = -1;
	do {
		if(is.eof()) return false;
		std::getline(is, buffer);
		if(!buffer.empty() && buffer[0] != '#') {
			iss.str(buffer); iss.clear();
			iss >> n;
			if(n < 0) return false;
		}
	} while(n < 0);

	Array<node> indexToNode(n);
	for(int i = 0; i < n; ) {
		if(is.eof()) return false;
		std::getline(is, buffer);

		if(!buffer.empty() && buffer[0] != '#') {
			node v = G.newNode();
			iss.str(buffer); iss.clear();
			iss >> gl.x(v) >> gl.y(v);
			indexToNode[i++] = v;
		}
	}

	while(!is.eof()) {
		std::getline(is, buffer);

		if(!buffer.empty() && buffer[0] != '#') {
			iss.str(buffer); iss.clear();
			int srcIndex, tgtIndex;

			if(iss.eof()) return false;
			iss >> srcIndex;
			if(srcIndex < 0 || srcIndex >= n) return false;

			if(iss.eof()) return false;
			iss >> tgtIndex;
			if(tgtIndex < 0 || tgtIndex >= n) return false;

			node src = indexToNode[srcIndex];
			node tgt = indexToNode[tgtIndex];
			edge e = G.newEdge(src,tgt);

			string symbol;
			if(iss.eof()) return false;
			iss >> symbol;
			if(symbol != "[") return false;

			IPolyline &ipl = gl.bends(e);;
			for(;;) {
				if(iss.eof()) return false;
				iss >> symbol;
				if(symbol == "]") break;

				IPoint ip;
				ip.m_x = atoi(symbol.c_str());
				if(iss.eof()) return false;
				iss >> ip.m_y;
				ipl.pushBack(ip);
			}
		}
	}

	return true;
}


bool GraphIO::writeChallengeGraph(const Graph &G, const GridLayout &gl, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeChallengeGraph(G, gl, os);
}

bool GraphIO::writeChallengeGraph(const Graph &G, const GridLayout &gl, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeChallengeGraph(G, gl, os);
}

bool GraphIO::writeChallengeGraph(const Graph &G, const GridLayout &gl, ostream &os)
{
	if(!os.good()) return false;

	os << "# Number of Nodes\n";
	os << G.numberOfNodes() << "\n";

	os << "# Nodes\n";
	NodeArray<int> index(G);
	int i = 0;
	node v;
	forall_nodes(v,G) {
		os << gl.x(v) << " " << gl.y(v) << "\n";
		index[v] = i++;
	}

	os << "# Edges\n";
	edge e;
	forall_edges(e,G) {
		os << index[e->source()] << " " << index[e->target()] << " [";
		const IPolyline &ipl = gl.bends(e);
		for(ListConstIterator<IPoint> it = ipl.begin(); it.valid(); ++it)
			os << " " << (*it).m_x << " " << (*it).m_y;
		os << " ]\n";
	}

	return true;
}


//---------------------------------------------------------
// Edge List Subgraph
//---------------------------------------------------------

bool GraphIO::readEdgeListSubgraph(Graph &G, List<edge> &delEdges, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) return false;
	return readEdgeListSubgraph(G, delEdges, is);
}

bool GraphIO::readEdgeListSubgraph(Graph &G, List<edge> &delEdges, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) return false;
	return readEdgeListSubgraph(G, delEdges, is);
}

bool GraphIO::readEdgeListSubgraph(Graph &G, List<edge> &delEdges, istream &is)
{
	G.clear();
	delEdges.clear();

	string buffer;

	if(is.eof()) return false;
	std::getline(is, buffer);
	istringstream iss(buffer);

	int n = 0, m = 0, m_del = 0;
	iss >> n >> m >> m_del;

	if(n < 0 || m < 0 || m_del < 0)
		return false;

	Array<node> indexToNode(n);
	for(int i = 0; i < n; ++i)
		indexToNode[i] = G.newNode();

	int m_all = m + m_del;
	for(int i = 0; i < m_all; ++i) {
		if(is.eof()) return false;

		std::getline(is, buffer);
		iss.str(buffer);
		iss.clear();

		int src = -1, tgt = -1;
		iss >> src >> tgt;
		if(src < 0 || src >= n || tgt < 0 || tgt >= n)
			return false;

		edge e = G.newEdge(indexToNode[src], indexToNode[tgt]);

		if(i >= m)
			delEdges.pushBack(e);
	}

	return true;
}


bool GraphIO::writeEdgeListSubgraph(const Graph &G, const List<edge> &delEdges, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeEdgeListSubgraph(G, delEdges, os);
}

bool GraphIO::writeEdgeListSubgraph(const Graph &G, const List<edge> &delEdges, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeEdgeListSubgraph(G, delEdges, os);
}

bool GraphIO::writeEdgeListSubgraph(const Graph &G, const List<edge> &delEdges, ostream &os)
{
	if(!os.good()) return false;

	const int m_del = delEdges.size();
	const int n = G.numberOfNodes();
	const int m = G.numberOfEdges() - m_del;

	os << n << " " << m << " " << m_del << "\n";

	EdgeArray<bool> markSub(G,true);
	for(ListConstIterator<edge> it = delEdges.begin(); it.valid(); ++it)
		markSub[*it] = false;

	NodeArray<int> index(G);
	int i = 0;
	node v;
	forall_nodes(v,G)
		index[v] = i++;

	edge e;
	forall_edges(e,G)
		if(markSub[e])
			os << index[e->source()] << " " << index[e->target()] << "\n";

	for(ListConstIterator<edge> it = delEdges.begin(); it.valid(); ++it)
		os << index[(*it)->source()] << " " << index[(*it)->target()] << "\n";


	return true;
}


//---------------------------------------------------------
// SteinLib parser
//---------------------------------------------------------

bool GraphIO::readSTP(
	EdgeWeightedGraph<double> &wG,
	List<node>           &terminals,
	NodeArray<bool>      &isTerminal,
	const char           *filename)
{
	ifstream is(filename);
	if(!is.is_open()) return false;
	return readSTP(wG, terminals, isTerminal, is);
}

bool GraphIO::readSTP(
	EdgeWeightedGraph<int> &wG,
	List<node>           &terminals,
	NodeArray<bool>      &isTerminal,
	const char           *filename)
{
	ifstream is(filename);
	if(!is.is_open()) return false;
	return readSTP(wG, terminals, isTerminal, is);
}

bool GraphIO::readSTP(
	EdgeWeightedGraph<double> &wG,
	List<node>           &terminals,
	NodeArray<bool>      &isTerminal,
	const string         &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) return false;
	return readSTP(wG, terminals, isTerminal, is);
}

bool GraphIO::readSTP(
	EdgeWeightedGraph<int> &wG,
	List<node>           &terminals,
	NodeArray<bool>      &isTerminal,
	const string         &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) return false;
	return readSTP(wG, terminals, isTerminal, is);
}

template<typename T>
static bool read_SteinLib(
	EdgeWeightedGraph<T> &wG,
	List<node>           &terminals,
	NodeArray<bool>      &isTerminal,
	istream              &is)
{
	wG.clear();
	terminals.clear();

	string buffer;
	std::istringstream iss;

	int section = 0;
	int nextSection = 1;

	string key, value;

	int n = 0;
	Array<node> indexToNode;
	//node root; // root terminal (directed case)

	string name, date, creator, remark;

	// 1. line = identifier
	if(!std::getline(is, buffer))
		return false;

	// these three variants are used in the SteinLib
	if (!equalIgnoreCase(buffer, "33D32945 STP File, STP Format Version 1.0")
	 && !equalIgnoreCase(buffer, "33D32945 STP File, STP Format Version 1.00")
	 && !equalIgnoreCase(buffer, "33d32945 STP File, STP Format Version  1.00")) {
		return false;
	}

	while (std::getline(is, buffer))
	{
		if (buffer.empty() || buffer[0] == '#' || buffer[0] == '\n') {
			continue;
		}

		iss.str(buffer); iss.clear();

		switch (section) {
		case 0:
			if (equalIgnoreCase(buffer, "SECTION Comment")
				&& nextSection == 1) {
					section = 1;
			} else if (equalIgnoreCase(buffer, "SECTION Graph")
				&& nextSection == 2) {
					section = 2;
			} else if (equalIgnoreCase(buffer, "SECTION Terminals")
				&& nextSection == 3) {
					section = 3;
			} else if (equalIgnoreCase(buffer, "SECTION Coordinates")
				&& nextSection == 4) {
					section = 4;
			} else if (equalIgnoreCase(buffer, "EOF") && nextSection >= 4) {
				return true;
			}
			break;

		case 1: // comment section
			iss >> key >> value;
			if (key == "Name") {
				name = value;
			} else if (equalIgnoreCase(key, "Date")) {
				date = value;
			} else if (equalIgnoreCase(key, "Creator")) {
				creator = value;
			} else if (equalIgnoreCase(key, "Remark")) {
				remark = value;
			} else if (equalIgnoreCase(key, "END")) {
				nextSection = 2;
				section = 0;
			} else {
				return false;
			}
			break;

		case 2: // graph section
			iss >> key;

			if(equalIgnoreCase(key, "END")) {
				nextSection = 3;
				section = 0;

			} else if(equalIgnoreCase(key, "Nodes")) {
				n = -1; iss >> n;
				if(n < 0)
					return false;

				indexToNode = Array<node>(1, n, 0);
				for (int i = 1; i <= n; i++) {
					indexToNode[i] = wG.newNode();
					isTerminal[indexToNode[i]] = false;
				}

			} else if (equalIgnoreCase(key,"E") || equalIgnoreCase(key,"A")) {
				int src = -1, tgt = -1, w = 1;
				iss >> src >> tgt >> w;

				if(src <= 0 || src > n || tgt <= 0 || tgt > n)
					return false;

				wG.newEdge(indexToNode[src], indexToNode[tgt], w);
			}

			// ignored keys: Edges, Arcs

			break;

		case 3: // terminals section
			iss >> key;
			if(equalIgnoreCase(key, "T")) {
				int v = -1;
				iss >> v;

				if(v <= 0 || v > n)
					return false;

				terminals.pushBack(indexToNode[v]);
				isTerminal[indexToNode[v]] = true;
			} else if (equalIgnoreCase(key, "END")) {
				nextSection = 4;
				section = 0;
			}

			// ignored keys: Root

			break;

		case 4: // coordinates section (omitted)
			iss >> key;
			if (equalIgnoreCase(key, "END")) {
				nextSection = 5;
				section = 0;
			}
			break;
		default:
			return false;
		}
	}
	return false;
}


bool GraphIO::readSTP(
	EdgeWeightedGraph<double> &wG,
	List<node>           &terminals,
	NodeArray<bool>      &isTerminal,
	istream              &is)
{
	return read_SteinLib(wG, terminals, isTerminal, is);
}

bool GraphIO::readSTP(
	EdgeWeightedGraph<int> &wG,
	List<node>           &terminals,
	NodeArray<bool>      &isTerminal,
	istream              &is)
{
	return read_SteinLib(wG, terminals, isTerminal, is);
}


//---------------------------------------------------------
// SVG graphics format
//---------------------------------------------------------

bool GraphIO::drawSVG(const GraphAttributes &A, const char *filename, const SVGSettings &settings)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return drawSVG(A, os, settings);
}

bool GraphIO::drawSVG(const GraphAttributes &A, const string &filename, const SVGSettings &settings)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return drawSVG(A, os, settings);
}

bool GraphIO::drawSVG(const ClusterGraphAttributes &A, const char *filename, const SVGSettings &settings)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return drawSVG(A, os, settings);
}

bool GraphIO::drawSVG(const ClusterGraphAttributes &A, const string &filename, const SVGSettings &settings)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return drawSVG(A, os, settings);
}


//---------------------------------------------------------
// Graph: GraphML format
//---------------------------------------------------------

bool GraphIO::readGraphML(Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) {
		return false;
	}
	return readGraphML(G, is);
}

bool GraphIO::readGraphML(Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) {
		return false;
	}
	return readGraphML(G, is);
}

bool GraphIO::readGraphML(Graph &G, istream &is)
{
	GraphMLParser parser(is);
	return parser.read(G);
}

bool GraphIO::writeGraphML(const Graph &G, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeGraphML(G, os);
}

bool GraphIO::writeGraphML(const Graph &G, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeGraphML(G, os);
}


//---------------------------------------------------------
// ClusterGraph: GraphML format
//---------------------------------------------------------

bool GraphIO::readGraphML(ClusterGraph &C, Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) {
		return false;
	}
	return readGraphML(C, G, is);
}

bool GraphIO::readGraphML(ClusterGraph &C, Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) {
		return false;
	}
	return readGraphML(C, G, is);
}

bool GraphIO::readGraphML(ClusterGraph &C, Graph &G, istream &is)
{
	GraphMLParser parser(is);
	return parser.read(G, C);
}

bool GraphIO::writeGraphML(const ClusterGraph &C, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeGraphML(C, os);
}

bool GraphIO::writeGraphML(const ClusterGraph &C, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeGraphML(C, os);
}


//---------------------------------------------------------
// GraphAttributes: GraphML format
//---------------------------------------------------------

bool GraphIO::readGraphML(GraphAttributes &A, Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) {
		return false;
	}
	return readGraphML(A, G, is);
}

bool GraphIO::readGraphML(GraphAttributes &A, Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) {
		return false;
	}
	return readGraphML(A, G, is);
}

bool GraphIO::readGraphML(GraphAttributes &A, Graph &G, istream &is)
{
	GraphMLParser parser(is);
	return parser.read(G, A);
}

bool GraphIO::writeGraphML(const GraphAttributes &A, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeGraphML(A, os);
}

bool GraphIO::writeGraphML(const GraphAttributes &A, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeGraphML(A, os);
}


//---------------------------------------------------------
// ClusterGraphAttributes: GraphML format
//---------------------------------------------------------

bool GraphIO::readGraphML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) {
		return false;
	}
	return readGraphML(A, C, G, is);
}

bool GraphIO::readGraphML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) {
		return false;
	}
	return readGraphML(A, C, G, is);
}

bool GraphIO::readGraphML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, istream &is)
{
	GraphMLParser parser(is);
	return parser.read(G, C, A);
}

bool GraphIO::writeGraphML(const ClusterGraphAttributes &A, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeGraphML(A, os);
}

bool GraphIO::writeGraphML(const ClusterGraphAttributes &A, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeGraphML(A, os);
}


//---------------------------------------------------------
// Graph: DOT format
//---------------------------------------------------------

bool GraphIO::readDOT(Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) {
		return false;
	}
	return readDOT(G, is);
}

bool GraphIO::readDOT(Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) {
		return false;
	}
	return readDOT(G, is);
}

bool GraphIO::readDOT(Graph &G, istream &is)
{
	dot::Parser parser(is);
	return parser.read(G);
}

bool GraphIO::writeDOT(const Graph &G, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeDOT(G, os);
}

bool GraphIO::writeDOT(const Graph &G, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeDOT(G, os);
}


//---------------------------------------------------------
// ClusterGraph: DOT format
//---------------------------------------------------------

bool GraphIO::readDOT(ClusterGraph &C, Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) {
		return false;
	}
	return readDOT(C, G, is);
}

bool GraphIO::readDOT(ClusterGraph &C, Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) {
		return false;
	}
	return readDOT(C, G, is);
}

bool GraphIO::readDOT(ClusterGraph &C, Graph &G, istream &is)
{
	dot::Parser parser(is);
	return parser.read(G, C);
}

bool GraphIO::writeDOT(const ClusterGraph &C, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeDOT(C, os);
}

bool GraphIO::writeDOT(const ClusterGraph &C, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeDOT(C, os);
}


//---------------------------------------------------------
// GraphAttributes: DOT format
//---------------------------------------------------------

bool GraphIO::readDOT(GraphAttributes &A, Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) {
		return false;
	}
	return readDOT(A, G, is);
}

bool GraphIO::readDOT(GraphAttributes &A, Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) {
		return false;
	}
	return readDOT(A, G, is);
}

bool GraphIO::readDOT(GraphAttributes &A, Graph &G, istream &is)
{
	dot::Parser parser(is);
	return parser.read(G, A);
}

bool GraphIO::writeDOT(const GraphAttributes &A, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeDOT(A, os);
}

bool GraphIO::writeDOT(const GraphAttributes &A, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeDOT(A, os);
}


//---------------------------------------------------------
// ClusterGraphAttributes: DOT format
//---------------------------------------------------------

bool GraphIO::readDOT(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) {
		return false;
	}
	return readDOT(A, C, G, is);
}

bool GraphIO::readDOT(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) {
		return false;
	}
	return readDOT(A, C, G, is);
}

bool GraphIO::readDOT(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, istream &is)
{
	dot::Parser parser(is);
	return parser.read(G, C, A);
}

bool GraphIO::writeDOT(const ClusterGraphAttributes &A, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeDOT(A, os);
}

bool GraphIO::writeDOT(const ClusterGraphAttributes &A, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeDOT(A, os);
}


//---------------------------------------------------------
// Graph: GEXF format
//---------------------------------------------------------

bool GraphIO::readGEXF(Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) {
		return false;
	}
	return readGEXF(G, is);
}

bool GraphIO::readGEXF(Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) {
		return false;
	}
	return readGEXF(G, is);
}

bool GraphIO::readGEXF(Graph &G, istream &is)
{
	gexf::Parser parser(is);
	return parser.read(G);
}

bool GraphIO::writeGEXF(const Graph &G, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeGEXF(G, os);
}

bool GraphIO::writeGEXF(const Graph &G, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeGEXF(G, os);
}


//---------------------------------------------------------
// ClusterGraph: GEXF format
//---------------------------------------------------------

bool GraphIO::readGEXF(ClusterGraph &C, Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) {
		return false;
	}
	return readGEXF(C, G, is);
}

bool GraphIO::readGEXF(ClusterGraph &C, Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) {
		return false;
	}
	return readGEXF(C, G, is);
}

bool GraphIO::readGEXF(ClusterGraph &C, Graph &G, istream &is)
{
	gexf::Parser parser(is);
	return parser.read(G, C);
}

bool GraphIO::writeGEXF(const ClusterGraph &C, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeGEXF(C, os);
}

bool GraphIO::writeGEXF(const ClusterGraph &C, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeGEXF(C, os);
}


//---------------------------------------------------------
// GraphAttributes: GEXF format
//---------------------------------------------------------

bool GraphIO::readGEXF(GraphAttributes &A, Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) {
		return false;
	}
	return readGEXF(A, G, is);
}

bool GraphIO::readGEXF(GraphAttributes &A, Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) {
		return false;
	}
	return readGEXF(A, G, is);
}

bool GraphIO::readGEXF(GraphAttributes &A, Graph &G, istream &is)
{
	gexf::Parser parser(is);
	return parser.read(G, A);
}

bool GraphIO::writeGEXF(const GraphAttributes &A, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeGEXF(A, os);
}

bool GraphIO::writeGEXF(const GraphAttributes &A, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeGEXF(A, os);
}


//---------------------------------------------------------
// ClusterGraphAttributes: GEXF format
//---------------------------------------------------------

bool GraphIO::readGEXF(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) {
		return false;
	}
	return readGEXF(A, C, G, is);
}

bool GraphIO::readGEXF(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) {
		return false;
	}
	return readGEXF(A, C, G, is);
}

bool GraphIO::readGEXF(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, istream &is)
{
	gexf::Parser parser(is);
	return parser.read(G, C, A);
}

bool GraphIO::writeGEXF(const ClusterGraphAttributes &A, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeGEXF(A, os);
}

bool GraphIO::writeGEXF(const ClusterGraphAttributes &A, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeGEXF(A, os);
}


//---------------------------------------------------------
// Graph: GDF format
//---------------------------------------------------------

bool GraphIO::readGDF(Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) {
		return false;
	}
	return readGDF(G, is);
}

bool GraphIO::readGDF(Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) {
		return false;
	}
	return readGDF(G, is);
}

bool GraphIO::readGDF(Graph &G, istream &is)
{
	gdf::Parser parser(is);
	return parser.read(G);
}

bool GraphIO::writeGDF(const Graph &G, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeGDF(G, os);
}

bool GraphIO::writeGDF(const Graph &G, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeGDF(G, os);
}


//---------------------------------------------------------
// GraphAttributes: GDF format
//---------------------------------------------------------

bool GraphIO::readGDF(GraphAttributes &A, Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) {
		return false;
	}
	return readGDF(A, G, is);
}

bool GraphIO::readGDF(GraphAttributes &A, Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) {
		return false;
	}
	return readGDF(A, G, is);
}

bool GraphIO::readGDF(GraphAttributes &A, Graph &G, istream &is)
{
	gdf::Parser parser(is);
	return parser.read(G, A);
}

bool GraphIO::writeGDF(const GraphAttributes &A, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeGDF(A, os);
}

bool GraphIO::writeGDF(const GraphAttributes &A, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeGDF(A, os);
}


//---------------------------------------------------------
// Graph: TLP format
//---------------------------------------------------------

bool GraphIO::readTLP(Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) {
		return false;
	}
	return readTLP(G, is);
}

bool GraphIO::readTLP(Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) {
		return false;
	}
	return readTLP(G, is);
}

bool GraphIO::readTLP(Graph &G, istream &is)
{
	tlp::Parser parser(is);
	return parser.read(G);
}

bool GraphIO::writeTLP(const Graph &G, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeTLP(G, os);
}

bool GraphIO::writeTLP(const Graph &G, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeTLP(G, os);
}


//---------------------------------------------------------
// ClusterGraph: TLP format
//---------------------------------------------------------

bool GraphIO::readTLP(ClusterGraph &C, Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) {
		return false;
	}
	return readTLP(C, G, is);
}

bool GraphIO::readTLP(ClusterGraph &C, Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) {
		return false;
	}
	return readTLP(C, G, is);
}

bool GraphIO::readTLP(ClusterGraph &C, Graph &G, istream &is)
{
	tlp::Parser parser(is);
	return parser.read(G, C);
}

bool GraphIO::writeTLP(const ClusterGraph &C, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeTLP(C, os);
}

bool GraphIO::writeTLP(const ClusterGraph &C, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeTLP(C, os);
}


//---------------------------------------------------------
// GraphAttributes: TLP format
//---------------------------------------------------------

bool GraphIO::readTLP(GraphAttributes &A, Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) {
		return false;
	}
	return readTLP(A, G, is);
}

bool GraphIO::readTLP(GraphAttributes &A, Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) {
		return false;
	}
	return readTLP(A, G, is);
}

bool GraphIO::readTLP(GraphAttributes &A, Graph &G, istream &is)
{
	tlp::Parser parser(is);
	return parser.read(G, A);
}

bool GraphIO::writeTLP(const GraphAttributes &A, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeTLP(A, os);
}

bool GraphIO::writeTLP(const GraphAttributes &A, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeTLP(A, os);
}


//---------------------------------------------------------
// ClusterGraphAttributes: TLP format
//---------------------------------------------------------

bool GraphIO::readTLP(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) {
		return false;
	}
	return readTLP(A, C, G, is);
}

bool GraphIO::readTLP(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) {
		return false;
	}
	return readTLP(A, C, G, is);
}

bool GraphIO::readTLP(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, istream &is)
{
	tlp::Parser parser(is);
	return parser.read(G, C, A);
}

bool GraphIO::writeTLP(const ClusterGraphAttributes &A, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeTLP(A, os);
}

bool GraphIO::writeTLP(const ClusterGraphAttributes &A, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeTLP(A, os);
}


//---------------------------------------------------------
// Graph: UCINET DL format
//---------------------------------------------------------

bool GraphIO::readDL(Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) {
		return false;
	}
	return readDL(G, is);
}

bool GraphIO::readDL(Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) {
		return false;
	}
	return readDL(G, is);
}

bool GraphIO::readDL(Graph &G, istream &is)
{
	DLParser parser(is);
	return parser.read(G);
}

bool GraphIO::writeDL(const Graph &G, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeDL(G, os);
}

bool GraphIO::writeDL(const Graph &G, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeDL(G, os);
}


//---------------------------------------------------------
// GraphAttributes: UCINET DL format
//---------------------------------------------------------

bool GraphIO::readDL(GraphAttributes &A, Graph &G, const char *filename)
{
	ifstream is(filename);
	if(!is.is_open()) {
		return false;
	}
	return readDL(A, G, is);
}

bool GraphIO::readDL(GraphAttributes &A, Graph &G, const string &filename)
{
	ifstream is(OGDF_STRING_OPEN(filename));
	if(!is.is_open()) {
		return false;
	}
	return readDL(A, G, is);
}

bool GraphIO::readDL(GraphAttributes &A, Graph &G, istream &is)
{
	DLParser parser(is);
	return parser.read(G, A);
}

bool GraphIO::writeDL(const GraphAttributes &A, const char *filename)
{
	ofstream os(filename);
	if(!os.is_open()) return false;
	return writeDL(A, os);
}

bool GraphIO::writeDL(const GraphAttributes &A, const string &filename)
{
	ofstream os(OGDF_STRING_OPEN(filename));
	if(!os.is_open()) return false;
	return writeDL(A, os);
}


} // end namespace ogdf

