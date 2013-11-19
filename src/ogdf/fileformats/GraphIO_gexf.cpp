/*
 * $Revision: 3837 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 15:19:30 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Implements GEXF write functionality of class GraphIO.
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
#include <ogdf/fileformats/GEXF.h>
#include <ogdf/fileformats/GraphML.h>


namespace ogdf {

namespace gexf {


static inline void writeHeader(std::ostream &out, bool viz = false)
{
	out << "<?xml "
	    << "version=\"1.0\" "
	    << "encoding=\"UTF-8\"?"
	    << ">\n";
	out << "<gexf "
	    << "xmlns=\"http://www.gexf.net/1.2draft\" "
	    << (viz ? "xmlns:viz=\"http://www.gexf.net/1.2draft/viz\"" : "")
	    << "version=\"1.2\""
	    << ">\n";

	// TODO: creator, description and date information.
}


static inline void writeFooter(std::ostream &out)
{
	out << "</gexf>\n";
}


template <typename T>
static inline void writeAttValue(
	std::ostream &out, int depth,
	const graphml::Attribute &attr, const T &value)
{
	const std::string &title = graphml::toString(attr);
	GraphIO::indent(out, depth) << "<attvalue "
	                            << "for=\"" << title << "\" "
	                            << "value=\"" << value << "\" "
	                            << "/>\n";
}


static inline void defineAttribute(
	std::ostream &out, int depth,
	const std::string &name, const std::string &type)
{
	GraphIO::indent(out, depth) << "<attribute "
	                            << "id=\"" << name << "\" "
	                            << "title=\"" << name << "\" "
	                            << "type=\"" << type << "\" "
	                            << "/>\n";
}


static inline void defineAttributes(
	std::ostream &out, int depth,
	const GraphAttributes &GA)
{
	const long attrs = GA.attributes();

	// Declare node attributes.
	GraphIO::indent(out, depth) << "<attributes class=\"node\">\n";
	if(attrs & GraphAttributes::nodeType) {
		defineAttribute(
			out, depth + 1,
			graphml::toString(graphml::a_nodeType), "string");
	}
	if(attrs & GraphAttributes::nodeTemplate) {
		defineAttribute(
			out, depth + 1,
			graphml::toString(graphml::a_template), "string");
	}
	if(attrs & GraphAttributes::nodeWeight) {
		defineAttribute(
			out, depth + 1,
			graphml::toString(graphml::a_nodeWeight), "float");
	}
	GraphIO::indent(out, depth) << "</attributes>\n";

	// Declare edge attributes.
	GraphIO::indent(out, depth) << "<attributes class=\"edge\">\n";
	if(attrs & GraphAttributes::edgeType) {
		defineAttribute(
			out, depth + 1,
			graphml::toString(graphml::a_edgeType), "string");
	}
	if(attrs & GraphAttributes::edgeArrow) {
		defineAttribute(
			out, depth + 1,
			graphml::toString(graphml::a_edgeArrow), "string");
	}
	GraphIO::indent(out, depth) << "</attributes>\n";
}


static inline void writeAttributes(
	std::ostream &out, int depth,
	const GraphAttributes &GA, node v)
{
	const long attrs = GA.attributes();

	if(attrs & GraphAttributes::nodeGraphics) {
		const double z = (attrs & GraphAttributes::threeD) ? GA.z(v) : 0.0;
		GraphIO::indent(out, depth) << "<viz:position "
		                            << "x=\"" << GA.x(v) << "\" "
		                            << "y=\"" << GA.y(v) << "\" "
		                            << "z=\"" << z << "\" "
		                            << "/>\n";

		// TODO: size is a scale here, so we have to know average size first.
		// const double size = std::max(GA.width(v), GA.height(v));
		// GraphIO::indent(out, depth) << "<viz:size "
		//                             << "value=\"" << size << "\" "
		//                             << "/>\n";

		const Shape shape = GA.shape(v);
		GraphIO::indent(out, depth) << "<viz:shape "
		                            << "value=\"" << toString(shape) << "\" "
		                            << "/>\n";
	}

	if(attrs & GraphAttributes::nodeStyle) {
		const Color &color = GA.fillColor(v);

		const int red = color.red();
		const int green = color.green();
		const int blue = color.blue();
		const int alpha = color.alpha();

		GraphIO::indent(out, depth) << "<viz:color "
		                            << "red=\"" << red << "\" "
		                            << "green=\"" << green << "\" "
		                            << "blue=\"" << blue << "\" "
		                            << "alpha=\"" << alpha << "\" "
		                            << "/>\n";

	}

	/*
	 * Node type, template and weight are not supported by VIZ module. So, they
	 * need to be written using <attvalues> tag (for estetic reasons, we write
	 * them only if either of them is present). For convenience reasons, we use
	 * the same names and values as in GraphML format.
	 */
	if(!(attrs & (GraphAttributes::nodeType |
	              GraphAttributes::nodeTemplate |
	              GraphAttributes::nodeWeight))) {
		return;
	}

	GraphIO::indent(out, depth) << "<attvalues>\n";

	if(attrs & GraphAttributes::nodeType) {
		writeAttValue(
			out, depth + 1,
			graphml::a_nodeType, graphml::toString(GA.type(v)));
	}

	if(attrs & GraphAttributes::nodeTemplate) {
		writeAttValue(out, depth + 1, graphml::a_template, GA.templateNode(v));
	}

	if(attrs & GraphAttributes::nodeWeight) {
		writeAttValue(out, depth + 1, graphml::a_nodeWeight, GA.weight(v));
	}

	GraphIO::indent(out, depth) << "</attvalues>\n";
}


static inline void writeAttributes(
	std::ostream &out, int depth,
	const GraphAttributes &GA, edge e)
{
	const long attrs = GA.attributes();

	if(attrs & GraphAttributes::edgeStyle) {
		const Color &color = GA.strokeColor(e);

		const int red = color.red();
		const int green = color.green();
		const int blue = color.blue();
		const int alpha = color.alpha();

		GraphIO::indent(out, depth) << "<viz:color "
		                            << "red=\"" << red << "\" "
		                            << "green=\"" << green << "\" "
		                            << "blue=\"" << blue << "\" "
		                            << "alpha=\"" << alpha << "\" "
		                            << "/>\n";
	}

	if(attrs & GraphAttributes::edgeDoubleWeight) {
		const double weight = GA.doubleWeight(e);
		GraphIO::indent(out, depth) << "<viz:thickness "
		                            << "value=\"" << weight << "\" "
		                            << "/>\n";
	} else if(attrs & GraphAttributes::edgeIntWeight) {
		const int weight = GA.intWeight(e);
		GraphIO::indent(out, depth) << "<viz:thickness "
		                            << "value=\"" << weight << "\" "
		                            << "/>\n";
	}

	/*
	 * Edge type and arrow are not supported by VIZ module. Therefore, they
	 * need to be written using <attvalues> tag (for estetic reasons, we write
	 * them only if either of them is present). For convenience reasons, we use
	 * the same names and values as in GraphML format.
	 */
	if(!(attrs & (GraphAttributes::edgeType | GraphAttributes::edgeArrow))) {
		return;
	}

	GraphIO::indent(out, depth) << "<attvalues>\n";

	if(attrs & GraphAttributes::edgeType) {
		writeAttValue(
			out, depth + 1,
			graphml::a_edgeType, graphml::toString(GA.type(e)));
	}
	if(attrs & GraphAttributes::edgeArrow) {
		writeAttValue(
			out, depth + 1,
			graphml::a_edgeArrow, graphml::toString(GA.arrowType(e)));
	}

	GraphIO::indent(out, depth) << "</attvalues>\n";
}


static inline void writeNode(
	std::ostream &out, int depth,
	const GraphAttributes *GA, node v)
{
	if(GA) {
		GraphIO::indent(out, depth) << "<node id=\"" << v->index() << "\"";
		if(GA->attributes() & GraphAttributes::nodeLabel) {
			out << " label=\"" << GA->label(v) << "\"";
		}
		out << ">\n";

		writeAttributes(out, depth + 1, *GA, v);

		GraphIO::indent(out, depth) << "</node>\n";
	} else {
		GraphIO::indent(out, depth) << "<node "
		                            << "id=\"" << v->index() << "\" "
		                            << "/>\n";
	}
}


static inline void writeEdge(
	std::ostream &out, int depth,
	const GraphAttributes *GA, edge e)
{
	if(GA) {
		GraphIO::indent(out, depth) << "<edge id=\"" << e->index() << "\"";
		if(GA->attributes() & GraphAttributes::edgeLabel) {
			out << " label=\"" << GA->label(e) << "\"";
		}
		out << ">\n";

		writeAttributes(out, depth + 1, *GA, e);

		GraphIO::indent(out, depth) << "</edge>\n";
	} else {
		GraphIO::indent(out, depth) << "<edge "
		                            << "id=\"" << e->index() << "\" "
		                            << "source=\"" << e->source() << "\" "
		                            << "target=\"" << e->target() << "\" "
		                            << "/>\n";
	}
}


static inline void writeEdges(
	std::ostream &out,
	const Graph &G, const GraphAttributes *GA)
{
	GraphIO::indent(out, 2) << "<edges>\n";

	edge e;
	forall_edges(e, G) {
		writeEdge(out, 3, GA, e);
	}

	GraphIO::indent(out, 2) << "</edges>\n";
}


static void writeCluster(
	std::ostream &out, int depth,
	const ClusterGraph &C, const ClusterGraphAttributes *CA, cluster c)
{
	if(C.rootCluster() != c) {
		GraphIO::indent(out, depth) << "<node "
		                            << "id=\"cluster" << c->index() << "\""
		                            << ">\n";
	} else {
		const std::string dir =
			(CA && !CA->directed()) ? "undirected" : "directed";
		GraphIO::indent(out, depth) << "<graph "
		                            << "mode=\"static\""
		                            << "defaultedgetype=\"" << dir << "\""
		                            << ">\n";

		if(CA) {
			defineAttributes(out, depth + 1, *CA);
		}
	}

	GraphIO::indent(out, depth + 1) << "<nodes>\n";

	for(ListConstIterator<cluster> cit = c->cBegin(); cit.valid(); cit++) {
		writeCluster(out, depth + 2, C, CA, *cit);
	}

	for(ListConstIterator<node> nit = c->nBegin(); nit.valid(); nit++) {
		writeNode(out, depth + 2, CA, *nit);
	}

	GraphIO::indent(out, depth + 1) << "</nodes>\n";

	if(C.rootCluster() != c) {
		GraphIO::indent(out, depth) << "</node>\n";
	} else {
		writeEdges(out, C.constGraph(), CA);
		GraphIO::indent(out, depth) << "</graph>\n";
	}
}


static void writeGraph(
	std::ostream &out,
	const Graph &G, const GraphAttributes *GA)
{
	const std::string dir =
		(GA && !GA->directed()) ? "undirected" : "directed";
	GraphIO::indent(out, 1) << "<graph "
	                        << "mode=\"static\" "
	                        << "defaultedgetype=\"" << dir << "\""
	                        << ">\n";

	if(GA) {
		defineAttributes(out, 2, *GA);
	}

	GraphIO::indent(out, 2) << "<nodes>\n";
	node v;
	forall_nodes(v, G) {
		writeNode(out, 3, GA, v);
	}
	GraphIO::indent(out, 2) << "</nodes>\n";

	gexf::writeEdges(out, G, GA);

	GraphIO::indent(out, 1) << "</graph>\n";
}


} // end namespace gexf


bool GraphIO::writeGEXF(const Graph &G, std::ostream &out)
{
	gexf::writeHeader(out);
	gexf::writeGraph(out, G, NULL);
	gexf::writeFooter(out);

	return true;
}


bool GraphIO::writeGEXF(const ClusterGraph &C, std::ostream &out)
{
	gexf::writeHeader(out);
	gexf::writeCluster(out, 1, C, NULL, C.rootCluster());
	gexf::writeFooter(out);

	return true;
}


bool GraphIO::writeGEXF(const GraphAttributes &GA, std::ostream &out)
{
	const Graph &G = GA.constGraph();

	gexf::writeHeader(out, true);
	gexf::writeGraph(out, G, &GA);
	gexf::writeFooter(out);

	return true;
}


bool GraphIO::writeGEXF(const ClusterGraphAttributes &CA, std::ostream &out)
{
	const ClusterGraph &C = CA.constClusterGraph();

	gexf::writeHeader(out, true);
	gexf::writeCluster(out, 1, C, &CA, C.rootCluster());
	gexf::writeFooter(out);

	return true;
}


} // end namespace ogdf

