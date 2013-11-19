/*
 * $Revision: 3816 $
 *
 * last checkin:
 *   $Author: beyer $
 *   $Date: 2013-10-29 17:07:00 +0100 (Di, 29. Okt 2013) $
 ***************************************************************/

/** \file
 * \brief Implements GML write functionality of class GraphIO.
 *
 * \author Carsten Gutwenger
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

#include <ogdf/fileformats/GraphIO.h>


namespace ogdf {


// begin of GML file
static void write_gml_header(ostream &os, bool directed)
{
	os << "Creator \"ogdf::GraphIO::writeGML\"\n";
	os << "graph [\n";
	GraphIO::indent(os,1) << "directed " << ((directed) ? 1 : 0) << "\n";
}


// end of GML file
static void write_gml_footer(ostream &os)
{
	os << "]\n"; // graph
}


// write graph structure
static void write_gml_graph(const Graph &G, ostream &os, NodeArray<int> &index)
{
	int nextId = 0;

	node v;
	forall_nodes(v,G) {
		GraphIO::indent(os,1) << "node [\n";
		GraphIO::indent(os,2) << "id " << (index[v] = nextId++) << "\n";
		GraphIO::indent(os,1) << "]\n"; // node
	}

	edge e;
	forall_edges(e,G) {
		GraphIO::indent(os,1) << "edge [\n";
		GraphIO::indent(os,2) << "source " << index[e->source()] << "\n";
		GraphIO::indent(os,2) << "target " << index[e->target()] << "\n";
		GraphIO::indent(os,1) << "]\n"; // edge
	}
}


const int c_maxLengthPerLine = 200;

static void writeLongString(ostream &os, const string &str)
{
	os << "\"";

	int num = 1;
	for(string::size_type p = 0; p < str.size(); ++p)
	{
		switch(str[p]) {
		case '\\':
			os << "\\\\";
			num += 2;
			break;
		case '\"':
			os << "\\\"";
			num += 2;
			break;

		// ignored white space
		case '\r':
		case '\n':
		case '\t':
			break;

		default:
			os << str[p];
			++num;
		}

		if(num >= c_maxLengthPerLine) {
			os << "\\\n";
			num = 0;
		}
	}

	os << "\"";
}

const char *arrow_str[4] = {
	"none", "last", "first", "both"
};

// write graph structure with attributes
static void write_gml_graph(const GraphAttributes &A, ostream &os, NodeArray<int> &index)
{
	const Graph &G = A.constGraph();
	int nextId = 0;

	os.setf(ios::showpoint);
	os.precision(10);

	node v;
	forall_nodes(v,G) {
		GraphIO::indent(os,1) << "node [\n";
		GraphIO::indent(os,2) << "id " << (index[v] = nextId++) << "\n";

		if (A.attributes() & GraphAttributes::nodeTemplate) {
			GraphIO::indent(os,2) << "template ";
			writeLongString(os, A.templateNode(v));
			os << "\n";
		}
		if (A.attributes() & GraphAttributes::nodeLabel) {
			GraphIO::indent(os,2) << "label ";
			writeLongString(os, A.label(v));
			os << "\n";
		}
		if (A.attributes() & GraphAttributes::nodeWeight) {
			GraphIO::indent(os,2) << "weight "  << A.weight(v) << "\n";
		}
		if (A.attributes() & GraphAttributes::nodeGraphics) {
			GraphIO::indent(os,2) << "graphics [\n";
			GraphIO::indent(os,3) << "x " << A.x     (v) << "\n";
			GraphIO::indent(os,3) << "y " << A.y     (v) << "\n";
			GraphIO::indent(os,3) << "w " << A.width (v) << "\n";
			GraphIO::indent(os,3) << "h " << A.height(v) << "\n";
			if (A.attributes() & GraphAttributes::nodeStyle)
			{
				GraphIO::indent(os,3) << "fill \"" << A.fillColor(v) << "\"\n";
				GraphIO::indent(os,3) << "line \"" << A.strokeColor (v) << "\"\n";
				GraphIO::indent(os,3) << "pattern \"" << A.fillPattern(v) << "\"\n";
				GraphIO::indent(os,3) << "stipple "   << A.strokeType (v) << "\n";
				GraphIO::indent(os,3) << "lineWidth " << A.strokeWidth(v) << "\n";
			}
			GraphIO::indent(os,3) << "type \"";

			switch (A.shape(v)) {
			case shRect:
				os << "rectangle";
				break;
			case shRoundedRect:
				os << "roundedRect";
				break;
			case shEllipse:
				os << "oval";
				break;
			case shTriangle:
				os << "triangle";
				break;
			case shPentagon:
				os << "pentagon";
				break;
			case shHexagon:
				os << "hexagon";
				break;
			case shOctagon:
				os << "octagon";
				break;
			case shRhomb:
				os << "rhomb";
				break;
			case shTrapeze:
				os << "trapeze";
				break;
			case shParallelogram:
				os << "parallelogram";
				break;
			case shInvTriangle:
				os << "invTriangle";
				break;
			case shInvTrapeze:
				os << "invTrapeze";
				break;
			case shInvParallelogram:
				os << "invParallelogram";
				break;
			case shImage:
				os << "image";
				break;
			}

			os << "\"\n";
			GraphIO::indent(os,2) << "]\n"; // graphics
		}

		GraphIO::indent(os,1) << "]\n"; // node
	}

	edge e;
	forall_edges(e,G) {
		GraphIO::indent(os,1) << "edge [\n";
		GraphIO::indent(os,2) << "source " << index[e->source()] << "\n";
		GraphIO::indent(os,2) << "target " << index[e->target()] << "\n";

		if (A.attributes() & GraphAttributes::edgeLabel){
			GraphIO::indent(os,2) << "label ";
			writeLongString(os, A.label(e));
			os << "\n";
		}
		if (A.attributes() & GraphAttributes::edgeType)
			GraphIO::indent(os,2) << "generalization " << A.type(e) << "\n";

		if (A.attributes() & GraphAttributes::edgeSubGraphs)
			GraphIO::indent(os,2) << "subgraph " << A.subGraphBits(e) << "\n";

		if (A.attributes() & GraphAttributes::edgeGraphics) {
			GraphIO::indent(os,2) << "graphics [\n";

			GraphIO::indent(os,3) << "type \"line\"\n";

			if (A.attributes() & GraphAttributes::edgeType) {
				if (A.attributes() & GraphAttributes::edgeArrow) {
					int ae = (int)A.arrowType(e);
					if(0 <= ae && ae < 4)
						GraphIO::indent(os,3) << "arrow \"" << arrow_str[ae] << "\"\n";
				} else {
					GraphIO::indent(os,3) << "arrow ";
					if (A.type(e) == Graph::generalization)
						os << "\"last\"\n";
					else
						os << "\"none\"\n";
				}

			} else { // GraphAttributes::edgeType not used
				GraphIO::indent(os,3) << "arrow ";
				if (A.directed()) {
					os << "\"last\"\n";
				} else {
					os << "\"none\"\n";
				}
			}

			if (A.attributes() & GraphAttributes::edgeStyle) {
				GraphIO::indent(os,3) << "stipple "   << A.strokeType(e) << "\n";
				GraphIO::indent(os,3) << "lineWidth " << A.strokeWidth(e) << "\n";
			}

			if (A.attributes() & GraphAttributes::edgeDoubleWeight) {
				GraphIO::indent(os,3) << "weight " << A.doubleWeight(e) << "\n";
			}

			const DPolyline &dpl = A.bends(e);
			if (!dpl.empty()) {
				GraphIO::indent(os,3) << "Line [\n";

				node v = e->source();
				if( dpl.front().m_x < A.x(v) - A.width (v)/2 ||
					dpl.front().m_x > A.x(v) + A.width (v)/2 ||
					dpl.front().m_y < A.y(v) - A.height(v)/2 ||
					dpl.front().m_y > A.y(v) + A.height(v)/2)
				{
					GraphIO::indent(os,4) << "point [ x " << A.x(e->source()) << " y " << A.y(e->source()) << " ]\n";
				}

				ListConstIterator<DPoint> it;
				for(it = dpl.begin(); it.valid(); ++it)
					GraphIO::indent(os,4) << "point [ x " << (*it).m_x << " y " << (*it).m_y << " ]\n";

				v = e->target();
				if( dpl.back().m_x < A.x(v) - A.width (v)/2 ||
					dpl.back().m_x > A.x(v) + A.width (v)/2 ||
					dpl.back().m_y < A.y(v) - A.height(v)/2 ||
					dpl.back().m_y > A.y(v) + A.height(v)/2)
				{
					GraphIO::indent(os,4) << "point [ x " << A.x(e->target()) << " y " << A.y(e->target()) << " ]\n";
				}

				GraphIO::indent(os,3) << "]\n"; // Line
			}//bends

			//output width and color
			if ((A.attributes() & GraphAttributes::edgeStyle))
				GraphIO::indent(os,3) << "fill \"" << A.strokeColor(e) << "\"\n";

			GraphIO::indent(os,2) << "]\n"; // graphics
		}

		GraphIO::indent(os,1) << "]\n"; // edge
	}
}


// write cluster structure
static void write_gml_cluster(cluster c, int d, ostream &os, const NodeArray<int> &index, int &nextClusterIndex)
{
	if (nextClusterIndex == 0)
		GraphIO::indent(os,d) << "rootcluster [\n";
	else
	{
		GraphIO::indent(os,d) << "cluster [\n";
		GraphIO::indent(os,d+1) << "id " << nextClusterIndex << "\n";
	}

	nextClusterIndex++;

	ListConstIterator<cluster> it;
	for (it = c->cBegin(); it.valid(); it++)
		write_gml_cluster(*it, d+1, os, index, nextClusterIndex);

	ListConstIterator<node> itn;
	for (itn = c->nBegin(); itn.valid(); itn++)
		GraphIO::indent(os,d+1) << "node " << index[*itn] << "\n";

	GraphIO::indent(os,d) << "]\n"; // cluster
}


// write cluster structure with attributes
static void write_gml_cluster(const ClusterGraphAttributes &A, cluster c, int d, ostream &os, const NodeArray<int> &index, int &nextClusterIndex)
{
	if (nextClusterIndex == 0)
		GraphIO::indent(os,d) << "rootcluster [\n";
	else
	{
		GraphIO::indent(os,d) << "cluster [\n";
		GraphIO::indent(os,d+1) << "id " << nextClusterIndex << "\n";

		const string &templStr = A.templateCluster(c);
		if(templStr.length() > 0) {
			// GDE extension: Write cluster template and custom attribute
			GraphIO::indent(os,d+1) << "template ";
			writeLongString(os, templStr);
			os << "\n";

			GraphIO::indent(os,d+1) << "label ";
			writeLongString(os, A.label(c));
			os << "\n";

		} else {
			GraphIO::indent(os,d+1) << "label \"" << A.label(c) << "\"\n";
		}

		GraphIO::indent(os,d+1) << "graphics [\n";

		double shiftPos;
		shiftPos = A.y(c);

		GraphIO::indent(os,d+2) << "x " << A.x(c) << "\n";
		GraphIO::indent(os,d+2) << "y " << shiftPos/*y(c->index())*/ << "\n";

		GraphIO::indent(os,d+2) << "width "     << A.width(c)       << "\n";
		GraphIO::indent(os,d+2) << "height "    << A.height(c)      << "\n";
		GraphIO::indent(os,d+2) << "fill \""    << A.fillColor(c)   << "\"\n";
		GraphIO::indent(os,d+2) << "pattern "   << A.fillPattern(c) << "\n";

		GraphIO::indent(os,d+2) << "color \""   << A.strokeColor(c)       << "\"\n";
		GraphIO::indent(os,d+2) << "lineWidth " << A.strokeWidth(c)   << "\n";

		if (A.strokeType(c) != stSolid)
			GraphIO::indent(os,d+2) << "stipple " << A.strokeType(c) << "\n";

		GraphIO::indent(os,d+2) << "style \"rectangle\"\n";

		GraphIO::indent(os,d+1) << "]\n"; //graphics
	}

	nextClusterIndex++;

	ListConstIterator<cluster> it;
	for (it = c->cBegin(); it.valid(); it++)
		write_gml_cluster(A, *it, d+1, os, index, nextClusterIndex);

	ListConstIterator<node> itn;
	for (itn = c->nBegin(); itn.valid(); itn++)
		GraphIO::indent(os,d+1) << "node " << index[*itn] << "\n";

	GraphIO::indent(os,d) << "]\n"; // cluster
}


// write Graph
bool GraphIO::writeGML(const Graph &G, ostream &os)
{
	write_gml_header(os,true);
	NodeArray<int> index(G);
	write_gml_graph(G, os, index);
	write_gml_footer(os);

	return true;
}


// write ClusterGraph
bool GraphIO::writeGML(const ClusterGraph &C, ostream &os)
{
	const Graph &G = C.constGraph();

	write_gml_header(os,true);
	NodeArray<int> index(G);
	write_gml_graph(G, os, index);
	write_gml_footer(os);

	int nextClusterIndex = 0;
	write_gml_cluster(C.rootCluster(), 1, os, index, nextClusterIndex);

	return true;
}


// write GraphAttributes
bool GraphIO::writeGML(const GraphAttributes &A, ostream &os)
{
	write_gml_header(os, A.directed());
	NodeArray<int> index(A.constGraph());
	write_gml_graph(A, os, index);
	write_gml_footer(os);

	return true;
}


// write ClusterGraphAttributes
bool GraphIO::writeGML(const ClusterGraphAttributes &A, ostream &os)
{
	write_gml_header(os, A.directed());
	NodeArray<int> index(A.constGraph());
	write_gml_graph(A, os, index);
	write_gml_footer(os);

	int nextClusterIndex = 0;
	write_gml_cluster(A, A.constClusterGraph().rootCluster(), 1, os, index, nextClusterIndex);

	return true;
}


}
