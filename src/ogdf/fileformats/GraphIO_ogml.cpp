/*
 * $Revision: 3366 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-04 16:13:53 +0200 (Do, 04. Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Implements OGML write functionality of class GraphIO.
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


// begin of OGML file
static void write_ogml_header(ostream &os)
{
	// XML declaration and OGML tag
	os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"; // Latin-1
	os << "<ogml xmlns=\"http://www.ogdf.net/ogml\">\n";
	GraphIO::indent(os,1) << "<graph>\n";
}


// end of OGML file
static void write_ogml_footer(ostream &os)
{
	GraphIO::indent(os,1) << "</graph>\n";
	os << "</ogml>\n";
}


static void write_ogml_graph_edges(const Graph &G, ostream &os)
{
	edge e;
	forall_edges(e,G) {
		GraphIO::indent(os,3) << "<edge id=\"e" << e->index() << "\">\n";
		GraphIO::indent(os,4) << "<source idRef=\"n" << e->source()->index() << "\" />\n";
		GraphIO::indent(os,4) << "<target idRef=\"n" << e->target()->index() << "\" />\n";
		GraphIO::indent(os,3) << "</edge>\n";
	}
}

// write graph structure
static void write_ogml_graph(const Graph &G, ostream &os)
{
	GraphIO::indent(os,2) << "<structure>\n";

	node v;
	forall_nodes(v,G) {
		GraphIO::indent(os,3) << "<node id=\"n" << v->index() << "\">\n";
		GraphIO::indent(os,3) << "</node>\n";
	}

	write_ogml_graph_edges(G, os);

	GraphIO::indent(os,2) << "</structure>\n";
}


// recursively write clusters and nodes
static void write_ogml_graph(cluster c, int level, ostream &os)
{
	if(level > 0) {
		GraphIO::indent(os,2+level) << "<node id=\"c" << c->index() << "\">\n";
	}

	ListConstIterator<node> itn;
	for (itn = c->nBegin(); itn.valid(); ++itn) {
		GraphIO::indent(os,3+level) << "<node id=\"n" << (*itn)->index() << "\">\n";
		GraphIO::indent(os,3+level) << "</node>\n";
	}

	ListConstIterator<cluster> itc;
	for (itc = c->cBegin(); itc.valid(); ++itc) {
		write_ogml_graph(*itc, level+1, os);
	}

	if(level > 0) {
		GraphIO::indent(os,2+level) << "</node>\n";
	}
}


// write cluster graph structure with clusters
static void write_ogml_graph(const ClusterGraph &C, ostream &os)
{
	GraphIO::indent(os,2) << "<structure>\n";

	write_ogml_graph(C.rootCluster(), 0, os);
	write_ogml_graph_edges(C.constGraph(), os);

	GraphIO::indent(os,2) << "</structure>\n";
}


// static helper method for exchanging X(HT)ML-tag specific chars
static string formatLabel(const string& labelText)
{
	size_t length = labelText.length();
	string formattedString;

	for (size_t i = 0; i < length; ++i) {
		char c = labelText[i];
		if (c == '<') {
			formattedString += "&lt;";
		} else {
			if (c == '>') {
				formattedString += "&gt;";
				if ((i+1 < length) && (labelText[i+1] != '\n'))
					formattedString += '\n';
			} else {
				formattedString += c;
			}
		}
	}
	return formattedString;
}


static void write_ogml_graph_edges(const GraphAttributes &A, ostream &os)
{
	const Graph &G = A.constGraph();

	edge e;
	forall_edges(e,G) {
		GraphIO::indent(os,3) << "<edge id=\"e" << e->index() << "\">\n";
		if (A.attributes() & GraphAttributes::edgeLabel) {
			GraphIO::indent(os,4) << "<label id=\"le" << e->index() << "\">\n";
			GraphIO::indent(os,5) << "<content>" << formatLabel(A.label(e)) << "</content>\n";
			GraphIO::indent(os,4) << "</label>\n";
		}
		GraphIO::indent(os,4) << "<source idRef=\"n" << e->source()->index() << "\" />\n";
		GraphIO::indent(os,4) << "<target idRef=\"n" << e->target()->index() << "\" />\n";
		GraphIO::indent(os,3) << "</edge>\n";
	}
}


// write graph structure with attributes
static void write_ogml_graph(const GraphAttributes &A, ostream &os)
{
	const Graph &G = A.constGraph();

	GraphIO::indent(os,2) << "<structure>\n";

	node v;
	forall_nodes(v,G) {
		GraphIO::indent(os,3) << "<node id=\"n" << v->index() << "\">\n";
		if (A.attributes() & GraphAttributes::nodeLabel) {
			GraphIO::indent(os,4) << "<label id=\"ln" << v->index() << "\">\n";
			GraphIO::indent(os,5) << "<content>" << formatLabel(A.label(v)) << "</content>\n";
			GraphIO::indent(os,4) << "</label>\n";
		}
		GraphIO::indent(os,3) << "</node>\n";
	}

	write_ogml_graph_edges(A, os);

	GraphIO::indent(os,2) << "</structure>\n";
}


// recursively write clusters and nodes
static void write_ogml_graph(const ClusterGraphAttributes &A, cluster c, int level, ostream &os)
{
	if(level > 0) {
		GraphIO::indent(os,2+level) << "<node id=\"c" << c->index() << "\">\n";
		if (A.attributes() & GraphAttributes::nodeLabel) {
			GraphIO::indent(os,4) << "<label id=\"lc" << c->index() << "\">\n";
			GraphIO::indent(os,5) << "<content>" << formatLabel(A.label(c)) << "</content>\n";
			GraphIO::indent(os,4) << "</label>\n";
		}
	}

	ListConstIterator<node> itn;
	for (itn = c->nBegin(); itn.valid(); ++itn) {
		node v = *itn;
		GraphIO::indent(os,3+level) << "<node id=\"n" << v->index() << "\">\n";
		if (A.attributes() & GraphAttributes::nodeLabel) {
			GraphIO::indent(os,4) << "<label id=\"ln" << v->index() << "\">\n";
			GraphIO::indent(os,5) << "<content>" << formatLabel(A.label(v)) << "</content>\n";
			GraphIO::indent(os,4) << "</label>\n";
		}
		GraphIO::indent(os,3+level) << "</node>\n";
	}

	ListConstIterator<cluster> itc;
	for (itc = c->cBegin(); itc.valid(); ++itc) {
		write_ogml_graph(*itc, level+1, os);
	}

	if(level > 0) {
		GraphIO::indent(os,2+level) << "</node>\n";
	}
}


// write cluster structure with attributes
static void write_ogml_graph(const ClusterGraphAttributes &A, ostream &os)
{
	GraphIO::indent(os,2) << "<structure>\n";

	write_ogml_graph(A, A.constClusterGraph().rootCluster(), 0, os);
	write_ogml_graph_edges(A, os);

	GraphIO::indent(os,2) << "</structure>\n";
}


static const char *edgeStyleToOGML(StrokeType edgeStyle)
{
	switch (edgeStyle)  {
	case stNone:
		return "none";
	case stSolid:
		return "solid";
	case stDash:
		return "dash";
	case stDot:
		return "dot";
	case stDashdot:
		return "dashDot";
	case stDashdotdot:
		return "dashDotDot";
	default:
		return "solid";
	}
}


static const char *fillPatternToOGML(FillPattern brushPattern)
{
	switch (brushPattern) {
	case fpNone:
		return "noFill";
	case fpSolid:
		return "solid";
	case fpDense1:
		return "dense1";
	case fpDense2:
		return "dense2";
	case fpDense3:
		return "dense3";
	case fpDense4:
		return "dense4";
	case fpDense5:
		return "dense5";
	case fpDense6:
		return "dense6";
	case fpDense7:
		return "dense7";
	case fpHorizontal:
		return "hor";
	case fpVertical:
		return "ver";
	case fpCross:
		return "cross";
	case fpBackwardDiagonal:
		return "bDiag";
	case fpForwardDiagonal:
		return "fDiag";
	case fpDiagonalCross:
		return "diagCross";
	default:
		return "solid";
	}//switch
}



static void write_ogml_layout_nodes_edges(const GraphAttributes &A, ostream &os)
{
	const Graph &G = A.constGraph();

	if (A.attributes() & (GraphAttributes::nodeGraphics | GraphAttributes::nodeStyle))
	{
		node v;
		forall_nodes(v,G) {
			GraphIO::indent(os,4) << "<nodeStyle idRef=\"n" << v->index() << "\">\n";

			if(A.attributes() & GraphAttributes::nodeGraphics) {
				GraphIO::indent(os,5) << "<location x=\"" << A.x(v)-0.5*A.width(v) << "\" y=\""<< A.y(v)-0.5*A.height(v) << "\" />\n";
				GraphIO::indent(os,5) << "<shape type=\"";
				switch (A.shape(v)) {
				case shRect:
					os << "rect";
					break;
				case shRoundedRect:
					os << "roundedRect";
					break;
				case shEllipse:
					os << "ellipse";
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
				os << "\" width=\"" << A.width(v) << "\" height=\"" << A.height(v) << "\" />\n";
			}

			if(A.attributes() & (GraphAttributes::nodeStyle)) {
				// fill-tag
				GraphIO::indent(os,5) << "<fill";

				// color-attribute of fill-tag
				os << " color=\"" << A.fillColor(v) << "\"";

				// pattern- and patternColor-attribute of fill-tag (closing)
				os << " pattern=\"" << fillPatternToOGML(A.fillPattern(v)) << "\" patternColor=\"" << A.fillBgColor(v) << "\" />\n";
				// line-tag
				GraphIO::indent(os,5) << "<line type=\"" << edgeStyleToOGML(A.strokeType(v)) <<  "\" width=\"" << A.strokeWidth(v) << "\""
					<< " color=\"" << A.strokeColor(v) << "\"";

				// closing fill-tag
				os << " />\n";
			}

			GraphIO::indent(os,4) << "</nodeStyle>\n";
		}
	}

	int pointId = 0;

	if (A.attributes() & (GraphAttributes::edgeGraphics | GraphAttributes::edgeStyle))
	{
		edge e;
		forall_edges(e,G) {
			GraphIO::indent(os,4) << "<edgeStyle idRef=\"e" << e->index() << "\">\n";

			if(A.attributes() & GraphAttributes::edgeStyle) {
				GraphIO::indent(os,5) << "<line ";
				if (A.attributes() & GraphAttributes::edgeStyle) {
					os << "type=\"" << edgeStyleToOGML(A.strokeType(e)) << "\" width=\"" << A.strokeWidth(e) << "\" ";
					os << "color=\"" << A.strokeColor(e) << "\" />\n";
				} else 	{
					os << " />\n";
				}
			}

			// TODO review the handling of edge arrows
			if(A.attributes() & GraphAttributes::edgeArrow)
			{
				switch(A.arrowType(e)) {
				case eaNone:
					GraphIO::indent(os,5) << "<sourceStyle type=\"none\" color=\"#000000\" size=\"1\" />\n";
					GraphIO::indent(os,5) << "<targetStyle type=\"none\" color=\"#000000\" size=\"1\" />\n";
					break;
				case eaLast:
					GraphIO::indent(os,5) << "<sourceStyle type=\"none\" color=\"#000000\" size=\"1\" />\n";
					GraphIO::indent(os,5) << "<targetStyle type=\"arrow\" color=\"#000000\" size=\"1\" />\n";
					break;
				case eaFirst:
					GraphIO::indent(os,5) << "<sourceStyle type=\"arrow\" color=\"#000000\" size=\"1\" />\n";
					GraphIO::indent(os,5) << "<targetStyle type=\"none\" color=\"#000000\" size=\"1\" />\n";
					break;
				case eaBoth:
					GraphIO::indent(os,5) << "<sourceStyle type=\"arrow\" color=\"#000000\" size=\"1\" />\n";
					GraphIO::indent(os,5) << "<targetStyle type=\"arrow\" color=\"#000000\" size=\"1\" />\n";
					break;
				case eaUndefined:
					// do nothing
					break;
				default:
					// do nothing
					break;
				}
			}

			// handling of points
			// TODO: Revise for new OGML specification
			const DPolyline &dpl = A.bends(e);
			if (!dpl.empty()) {
				// handle source
				node v = e->source();
				if(dpl.front().m_x < A.x(v) - A.width(v)/2 ||
					dpl.front().m_x > A.x(v) + A.width(v)/2 ||
					dpl.front().m_y < A.y(v) - A.height(v)/2 ||
					dpl.front().m_y > A.y(v) + A.height(v)/2)	{
						GraphIO::indent(os,5) << "<point id=\"p" << pointId++ << "\" x=\"" << A.x(e->source()) << "\" y=\"" << A.y(e->source()) << "\" />\n";
				}
				// handle points
				ListConstIterator<DPoint> it;
				for(it = dpl.begin(); it.valid(); ++it) {
					GraphIO::indent(os,5) << "<point id=\"p" << pointId++ << "\" x=\"" << (*it).m_x << "\" y=\"" << (*it).m_y << "\" />\n";
				}
				// handle target
				v = e->target();
				if(dpl.back().m_x < A.x(v) - A.width(v)/2 ||
					dpl.back().m_x > A.x(v) + A.width(v)/2 ||
					dpl.back().m_y < A.y(v) - A.height(v)/2 ||
					dpl.back().m_y > A.y(v) + A.height(v)/2) {
						GraphIO::indent(os,5) << "<point id=\"p" << pointId++ << "\" x=\"" << A.x(e->target()) << "\" y=\"" << A.y(e->target()) << "\" />\n";
				}
			}

			GraphIO::indent(os,4) << "</edgeStyle>\n";
		}
	}
}


// write graph layout with attributes
static void write_ogml_layout(const GraphAttributes &A, ostream &os)
{
	GraphIO::indent(os,2) << "<layout>\n";
	GraphIO::indent(os,3) << "<styles>\n";

	write_ogml_layout_nodes_edges(A,os);

	GraphIO::indent(os,3) << "</styles>\n";
	GraphIO::indent(os,2) << "</layout>\n";
}


// write cluster layout with attributes
static void write_ogml_layout(const ClusterGraphAttributes &A, ostream &os)
{
	const ClusterGraph &C = A.constClusterGraph();

	GraphIO::indent(os,2) << "<layout>\n";
	GraphIO::indent(os,3) << "<styles>\n";

	cluster c;
	forall_clusters(c, C) {
		if(c != C.rootCluster()) {
			GraphIO::indent(os,4) << "<nodeStyle idRef=\"c" << c->index() << "\">\n";

			GraphIO::indent(os,5) << "<location x=\"" << A.x(c) << "\" y=\"" << A.y(c) << "\" />\n";
			GraphIO::indent(os,5) << "<shape type=\"rect\" width=\"" << A.width(c) << "\" height=\"" << A.height(c) << "\" />\n";
			GraphIO::indent(os,5) << "<fill color=\"" << A.fillColor(c) << "\""
				<< " pattern=\"" << fillPatternToOGML(A.fillPattern(c))
				<< "\" patternColor=\""
				<< A.fillBgColor(c) << "\" />\n";
			GraphIO::indent(os,5) << "<line type=\"" << edgeStyleToOGML(A.strokeType(c)) << "\" width=\"" << A.strokeWidth(c) << "\" color=\"" << A.strokeColor(c) << "\" />\n";

			GraphIO::indent(os,4) << "</nodeStyle>\n";
		}
	}

	write_ogml_layout_nodes_edges(A,os);

	GraphIO::indent(os,3) << "</styles>\n";
	GraphIO::indent(os,2) << "</layout>\n";
}


// write Graph
bool GraphIO::writeOGML(const Graph &G, ostream &os)
{
	write_ogml_header(os);
	write_ogml_graph(G, os);
	write_ogml_footer(os);

	return true;
}


// write ClusterGraph
bool GraphIO::writeOGML(const ClusterGraph &C, ostream &os)
{
	write_ogml_header(os);
	write_ogml_graph(C, os);
	write_ogml_footer(os);

	return true;
}


// write GraphAttributes
bool GraphIO::writeOGML(const GraphAttributes &A, ostream &os)
{
	write_ogml_header(os);
	write_ogml_graph(A, os);
	write_ogml_layout(A, os);
	write_ogml_footer(os);

	return true;
}


// write ClusterGraphAttributes
bool GraphIO::writeOGML(const ClusterGraphAttributes &A, ostream &os)
{
	write_ogml_header(os);
	write_ogml_graph(A, os);
	write_ogml_layout(A, os);
	write_ogml_footer(os);

	return true;
}


}
