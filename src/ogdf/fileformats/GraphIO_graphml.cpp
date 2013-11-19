/*
 * $Revision: 3837 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 15:19:30 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Implements GraphML write functionality of class GraphIO.
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

#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/fileformats/GraphML.h>


namespace ogdf {


static inline void writeGraphMLHeader(std::ostream &out)
{
	const std::string xmlns = "http://graphml.graphdrawing.org/xmlns";
	out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
	out << "<graphml xmlns=\"" << xmlns << "\"\n"
	       "         xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
	       "         xsi:schemaLocation=\"" << xmlns << "\n"
	       "                             " << xmlns << "/1.0/graphml.xsd\">\n";
}


static inline void writeGraphMLFooter(std::ostream &out)
{
	out << "</graphml>\n";
}


static inline void defineGraphMLAttribute(
	std::ostream &out,
	const std::string &kind, const std::string &name, const std::string &type)
{
	GraphIO::indent(out, 1) << "<key "
	                        << "for=\"" << kind << "\" "
	                        << "attr.name=\"" << name << "\" "
	                        << "attr.type=\"" << type << "\" "
	                        << "id=\"" << name << "\" />\n";
}


static inline void defineGraphMLAttributes(std::ostream &out, long attributes)
{
	using namespace graphml;

	// Gephi-compatible attributes
	if(attributes & GraphAttributes::nodeLabel) {
		defineGraphMLAttribute(out, "node", toString(a_nodeLabel), "string");
	}

	if(attributes & GraphAttributes::nodeGraphics) {
		defineGraphMLAttribute(out, "node", toString(a_x), "float");
		defineGraphMLAttribute(out, "node", toString(a_y), "float");
		defineGraphMLAttribute(out, "node", toString(a_size), "float");
	}

	if(attributes & GraphAttributes::nodeStyle) {
		defineGraphMLAttribute(out, "node", toString(a_r), "int");
		defineGraphMLAttribute(out, "node", toString(a_g), "int");
		defineGraphMLAttribute(out, "node", toString(a_b), "int");
	}

	if(attributes & GraphAttributes::edgeLabel) {
		defineGraphMLAttribute(out, "edge", toString(a_edgeLabel), "string");
	}

	if(attributes & GraphAttributes::edgeDoubleWeight) {
		defineGraphMLAttribute(out, "edge", toString(a_edgeWeight), "double");
	} else if(attributes & GraphAttributes::edgeIntWeight) {
		defineGraphMLAttribute(out, "edge", toString(a_edgeWeight), "int");
	}

	// OGDF-specific attributes.
	if (attributes & GraphAttributes::nodeGraphics) {
		defineGraphMLAttribute(out, "node", toString(a_width), "double");
		defineGraphMLAttribute(out, "node", toString(a_height), "double");
		defineGraphMLAttribute(out, "node", toString(a_shape), "string");
	}

	if(attributes & GraphAttributes::nodeStyle) {
		defineGraphMLAttribute(out, "node", toString(a_nodeStroke), "string");
	}

	if(attributes & GraphAttributes::nodeWeight) {
		defineGraphMLAttribute(out, "node", toString(a_nodeWeight), "int");
	}

	if(attributes & GraphAttributes::nodeType) {
		defineGraphMLAttribute(out, "node", toString(a_nodeType), "string");
	}

	if(attributes & GraphAttributes::nodeTemplate) {
		defineGraphMLAttribute(out, "node", toString(a_template), "string");
	}

	if(attributes & GraphAttributes::threeD) {
		defineGraphMLAttribute(out, "node", toString(a_z), "float");
	}

	if(attributes & GraphAttributes::edgeGraphics) {
		// Currently bending points are printed as list. More XML-ish way has
		// to be adopted in future (it will probably involve writing custom
		// XML schema...).
		defineGraphMLAttribute(out, "edge", toString(a_edgeBends), "string");
	}

	if(attributes & GraphAttributes::edgeType) {
		defineGraphMLAttribute(out, "edge", toString(a_edgeType), "string");
	}

	if(attributes & GraphAttributes::edgeArrow) {
		defineGraphMLAttribute(out, "edge", toString(a_edgeArrow), "string");
	}

	if(attributes & GraphAttributes::edgeStyle) {
		defineGraphMLAttribute(out, "edge", toString(a_edgeStroke), "string");
	}

	if(attributes & GraphAttributes::edgeSubGraphs) {
		defineGraphMLAttribute(out, "edge", toString(a_edgeSubGraph), "int");
	}
}


template <typename T>
static inline void writeGraphMLAttribute(
	std::ostream &out, int depth,
	const std::string &name, const T &value)
{
	GraphIO::indent(out, depth) << "<data key=\"" << name << "\">"
	                            << value
	                            << "</data>\n";
}


static inline void writeGraphMLNode(
	std::ostream &out, int depth,
	const Graph &G, const node &v)
{
	GraphIO::indent(out, depth) << "<node id=\"" << v->index() << "\" />\n";
}


static inline void writeGraphMLEdge(
	std::ostream &out, int depth,
	const Graph &G, const edge &e)
{
	const node &s = e->source();
	const node &t = e->target();
	GraphIO::indent(out, depth) << "<edge "
	                            << "id=\"" << e->index() << "\" "
	                            << "source=\"" << s->index() << "\" "
	                            << "target=\"" << t->index() << "\" "
	                            << "/>\n";
}


static inline void writeGraphMLNode(
	std::ostream &out, int depth,
	const GraphAttributes &GA, const node &v)
{
	using namespace graphml;

	// Use attribute id if avaliable, node index if not.
	GraphIO::indent(out, depth++) << "<node id=\"";
	if(GA.attributes() & GraphAttributes::nodeId) {
		out << GA.idNode(v);
	} else {
		out << v->index();
	}
	out << "\">\n";

	if(GA.attributes() & GraphAttributes::nodeLabel && GA.label(v) != "") {
		writeGraphMLAttribute(out, depth, toString(a_nodeLabel), GA.label(v));
	}

	if(GA.attributes() & GraphAttributes::nodeGraphics) {
		writeGraphMLAttribute(out, depth, toString(a_x), GA.x(v));
		writeGraphMLAttribute(out, depth, toString(a_y), GA.y(v));
		writeGraphMLAttribute(out, depth, toString(a_width), GA.width(v));
		writeGraphMLAttribute(out, depth, toString(a_height), GA.height(v));

		const double size = std::max(GA.width(v), GA.height(v));
		writeGraphMLAttribute(out, depth, toString(a_size), size);

		writeGraphMLAttribute(
			out, depth,
			toString(a_shape), toString(GA.shape(v)));
	}

	if(GA.attributes() & GraphAttributes::threeD) {
		writeGraphMLAttribute(out, depth, toString(a_z), GA.z(v));
	}

	if(GA.attributes() & GraphAttributes::nodeStyle) {
		const Color &col = GA.fillColor(v);
		writeGraphMLAttribute(
			out, depth,
			toString(a_r), static_cast<int>(col.red()));
		writeGraphMLAttribute(
			out, depth,
			toString(a_g), static_cast<int>(col.green()));
		writeGraphMLAttribute(
			out, depth,
			toString(a_b), static_cast<int>(col.blue()));

		writeGraphMLAttribute(
			out, depth,
			toString(a_nodeStroke), GA.strokeColor(v));
	}

	if(GA.attributes() & GraphAttributes::nodeType) {
		writeGraphMLAttribute(
			out, depth,
			toString(a_nodeType), toString(GA.type(v)));
	}

	if(GA.attributes() & GraphAttributes::nodeTemplate &&
	  GA.templateNode(v).length() > 0)
	{
		writeGraphMLAttribute(
			out, depth,
			toString(a_template), GA.templateNode(v));
	}

	if(GA.attributes() & GraphAttributes::nodeWeight) {
		writeGraphMLAttribute(out, depth, toString(a_nodeWeight), GA.weight(v));
	}

	GraphIO::indent(out, --depth) << "</node>\n";
}


static inline void writeGraphMLEdge(
	std::ostream &out, int depth,
	const GraphAttributes &GA, const edge &e)
{
	using namespace graphml;

	const node &s = e->source();
	const node &t = e->target();
	GraphIO::indent(out, depth++) << "<edge "
	                              << "id=\"" << e->index() << "\" "
	                              << "source=\"" << s->index() << "\" "
	                              << "target=\"" << t->index() << "\""
	                              << ">\n";

	if(GA.attributes() & GraphAttributes::edgeLabel && GA.label(e) != "") {
		writeGraphMLAttribute(out, depth, toString(a_edgeLabel), GA.label(e));
	}

	if(GA.attributes() & GraphAttributes::edgeDoubleWeight) {
		writeGraphMLAttribute(
			out, depth,
			toString(a_edgeWeight), GA.doubleWeight(e));
	} else if(GA.attributes() & GraphAttributes::edgeIntWeight) {
		writeGraphMLAttribute(
			out, depth,
			toString(a_edgeWeight), GA.intWeight(e));
	}

	if(GA.attributes() & GraphAttributes::edgeGraphics) {
		std::stringstream sstream; // For code consistency.

		forall_listiterators(DPoint, it, GA.bends(e)) {
			const DPoint &p = *it;
			sstream << p.m_x << " " << p.m_y << " ";
		}

		// Call to GA.bends(e).length() causes a "Segmentation fault".
		const std::string str = sstream.str();
		if (str.length() > 0) {
			writeGraphMLAttribute(
				out, depth,
				toString(a_edgeBends), sstream.str());
		}
	}

	if(GA.attributes() & GraphAttributes::edgeType) {
		writeGraphMLAttribute(
			out, depth,
			toString(a_edgeType), toString(GA.type(e)));
	}

	if(GA.attributes() & GraphAttributes::edgeArrow) {
		const EdgeArrow &arrow = GA.arrowType(e);
		if (arrow != eaUndefined) {
			writeGraphMLAttribute(
				out, depth,
				toString(a_edgeArrow), toString(arrow));
		}
	}

	if(GA.attributes() & GraphAttributes::edgeStyle) {
		writeGraphMLAttribute(
			out, depth,
			toString(a_edgeStroke), GA.strokeColor(e));
	}

	if(GA.attributes() & GraphAttributes::edgeSubGraphs) {
		const __uint32 mask = GA.subGraphBits(e);

		// Iterate over all subgraphs and print avaliable.
		for(size_t sg = 0; sg < sizeof(mask) * 8; sg++) {
			if((1 << sg) & mask) {
				writeGraphMLAttribute(out, depth, toString(a_edgeSubGraph), sg);
			}
		}
	}

	GraphIO::indent(out, --depth) << "</edge>\n";
}


static void writeGraphMLCluster(
	std::ostream &out, int depth,
	const ClusterGraph &C, const cluster &c, int &clusterId)
{
	const bool isRoot = C.rootCluster() == c;

	if(!isRoot) {
		GraphIO::indent(out, depth++) << "<node "
		                              << "id=\"cluster" << clusterId << "\""
		                              << ">\n";
		GraphIO::indent(out, depth++) << "<graph "
		                              << "id=\"cluster" << clusterId << ":\" "
		                              << "edgedefault=\"directed\""
		                              << ">\n";
	}
	clusterId++;

	for(ListConstIterator<cluster> cit = c->cBegin(); cit.valid(); cit++) {
		writeGraphMLCluster(out, depth, C, *cit, clusterId);
	}

	for(ListConstIterator<node> nit = c->nBegin(); nit.valid(); nit++) {
		writeGraphMLNode(out, depth, C, *nit);
	}

	if(!isRoot) {
		GraphIO::indent(out, --depth) << "</graph>\n";
		GraphIO::indent(out, --depth) << "</node>\n";
	}
}


static void writeGraphMLCluster(
	std::ostream &out, int depth,
	const ClusterGraphAttributes &CA, const cluster &c, int &clusterId)
{
	const bool isRoot = CA.constClusterGraph().rootCluster() == c;
	const std::string edgeDefault = CA.directed() ? "directed" : "undirected";

	if(!isRoot) {
		GraphIO::indent(out, depth++) << "<node "
		                              << "id=\"cluster" << clusterId << "\""
		                              << ">\n";
		GraphIO::indent(out, depth++) << "<graph "
		                              << "id=\"cluster" << clusterId << ":\" "
		                              << "edgedefault=\"" << edgeDefault << "\""
		                              << ">\n";
	}
	clusterId++;

	for(ListConstIterator<cluster> cit = c->cBegin(); cit.valid(); cit++) {
		writeGraphMLCluster(out, depth, CA, *cit, clusterId);
	}

	for(ListConstIterator<node> nit = c->nBegin(); nit.valid(); nit++) {
		writeGraphMLNode(out, depth, CA, *nit);
	}

	// There should be no attributes for root cluster.
	if(isRoot) {
		return;
	}

	GraphIO::indent(out, --depth) << "</graph>\n";

	using namespace graphml;

	// Writing cluster attributes (defined as cluster-node attributes).
	if(CA.label(c).length() > 0) {
		writeGraphMLAttribute(
			out, depth,
			toString(a_nodeLabel), CA.label(c));
	}
	writeGraphMLAttribute(out, depth, toString(a_x), CA.x(c));
	writeGraphMLAttribute(out, depth, toString(a_y), CA.y(c));

	const Color &col = CA.fillColor(c);
	writeGraphMLAttribute(
		out, depth,
		toString(a_r), static_cast<int>(col.red()));
	writeGraphMLAttribute(
		out, depth,
		toString(a_g), static_cast<int>(col.green()));
	writeGraphMLAttribute(
		out, depth,
		toString(a_b), static_cast<int>(col.blue()));
	writeGraphMLAttribute(
		out, depth,
		toString(a_clusterStroke), CA.strokeColor(c));

	if(CA.templateCluster(c).length() > 0) {
		writeGraphMLAttribute(
			out, depth,
			toString(a_template), CA.templateCluster(c));
	}

	// TODO: not important cluster attributes like fill patterns, background
	// color, stroke width, etc.

	GraphIO::indent(out, --depth) << "</node>\n";
}


bool GraphIO::writeGraphML(const Graph &G, std::ostream &out)
{
	writeGraphMLHeader(out);
	GraphIO::indent(out, 1) << "<graph id=\"G\" edgedefault=\"directed\">\n";

	node v;
	forall_nodes(v, G) {
		writeGraphMLNode(out, 2, G, v);
	}

	edge e;
	forall_edges(e, G) {
		writeGraphMLEdge(out, 2, G, e);
	}

	GraphIO::indent(out, 1) << "</graph>\n";
	writeGraphMLFooter(out);

	return true;
}


bool GraphIO::writeGraphML(const ClusterGraph &C, std::ostream &out)
{
	const Graph &G = C.constGraph();

	writeGraphMLHeader(out);
	GraphIO::indent(out, 1) << "<graph id=\"G\" edgedefault=\"directed\">\n";

	int clusterId = 0;
	writeGraphMLCluster(out, 2, G, C.rootCluster(), clusterId);

	edge e;
	forall_edges(e, G) {
		writeGraphMLEdge(out, 2, G, e);
	}

	GraphIO::indent(out, 1) << "</graph>\n";
	writeGraphMLFooter(out);

	return true;
}


bool GraphIO::writeGraphML(const GraphAttributes &GA, std::ostream &out)
{
	const Graph &G = GA.constGraph();
	const std::string edgeDefault = GA.directed() ? "directed" : "undirected";

	writeGraphMLHeader(out);
	defineGraphMLAttributes(out, GA.attributes());
	GraphIO::indent(out, 1) << "<graph "
	                        << "id=\"G\" "
	                        << "edgedefault=\"" << edgeDefault << "\""
	                        << ">\n";

	node v;
	forall_nodes(v, G) {
		writeGraphMLNode(out, 2, GA, v);
	}

	edge e;
	forall_edges(e, G) {
		writeGraphMLEdge(out, 2, GA, e);
	}

	GraphIO::indent(out, 1) << "</graph>\n";
	writeGraphMLFooter(out);

	return true;
}


bool GraphIO::writeGraphML(const ClusterGraphAttributes &CA, std::ostream &out)
{
	const Graph &G = CA.constGraph();
	const ClusterGraph &C = CA.constClusterGraph();

	writeGraphMLHeader(out);
	defineGraphMLAttributes(out, CA.attributes());
	defineGraphMLAttribute(out, "node", toString(graphml::a_clusterStroke), "string");

	GraphIO::indent(out, 1) << "<graph id=\"G\" edgedefault=\"directed\">\n";

	int clusterId = 0;
	writeGraphMLCluster(out, 2, CA, C.rootCluster(), clusterId);

	edge e;
	forall_edges(e, G) {
		writeGraphMLEdge(out, 2, CA, e);
	}

	GraphIO::indent(out, 1) << "</graph>\n";
	writeGraphMLFooter(out);

	return true;
}


}
