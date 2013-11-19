/*
 * $Revision: 3837 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 15:19:30 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of GEXF format parsing utilities.
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

#include <ogdf/fileformats/GEXF.h>
#include <ogdf/fileformats/GexfParser.h>
#include <ogdf/fileformats/GraphML.h>


namespace ogdf {

namespace gexf {


Parser::Parser(std::istream &is) : m_xml(is), m_nodeId(NULL), m_clusterId(NULL)
{
}


static inline bool readAttrDefs(
	HashArray<std::string, std::string> &attrMap,
	const XmlTagObject &attrsTag)
{
	List<XmlTagObject *> attrTags;
	attrsTag.findSonXmlTagObjectByName("attributes", attrTags);
	forall_listiterators(XmlTagObject *, jt, attrTags) {
		const XmlTagObject &attrTag = **jt;

		XmlAttributeObject *idAttr, *titleAttr;
		attrTag.findXmlAttributeObjectByName("id", idAttr);
		attrTag.findXmlAttributeObjectByName("title", titleAttr);

		if(!idAttr || !titleAttr) {
			std::cerr << "ERROR: \"id\" or \"title\" not found for attribute "
			          << "(line " << attrTag.getLine() << ").\n";
			return false;
		}

		attrMap[idAttr->getValue()] = titleAttr->getValue();
	}

	return true;
}

bool Parser::init()
{
	m_nodeId.clear();
	m_clusterId.clear();
	m_nodeAttr.clear();
	m_edgeAttr.clear();

	m_xml.createParseTree();

	const XmlTagObject gexfTag = m_xml.getRootTag();
	if(gexfTag.getName() != "gexf") {
		std::cerr << "ERROR: Root tag must be \"gexf\".\n";
		return false;
	}

	gexfTag.findSonXmlTagObjectByName("graph", m_graphTag);
	if(!m_graphTag) {
		std::cerr << "ERROR: Expected \"graph\" tag.\n";
		return false;
	}

	m_graphTag->findSonXmlTagObjectByName("nodes", m_nodesTag);
	if(!m_nodesTag) {
		std::cerr << "ERROR: No \"nodes\" tag found in graph.\n";
		return false;
	}

	m_graphTag->findSonXmlTagObjectByName("edges", m_edgesTag);
	if(!m_edgesTag) {
		std::cerr << "ERROR: No \"edges\" tag found in graph.\n";
		return false;
	}

	// Read attributes definitions. Could be lazily read later only
	// if GraphAttributes is given.
	List<XmlTagObject *> attrsTags;
	m_graphTag->findSonXmlTagObjectByName("attributes", attrsTags);
	forall_listiterators(XmlTagObject *, it, attrsTags) {
		const XmlTagObject &attrsTag = **it;

		XmlAttributeObject *classAttr;
		attrsTag.findXmlAttributeObjectByName("class", classAttr);
		if(!classAttr) {
			std::cerr << "ERROR: attributes tag is missing a class "
			          << "(line " << attrsTag.getLine() << ").\n";
			return false;
		}

		HashArray<std::string, std::string> *attrMap;
		if(classAttr->getValue() == "node") {
			attrMap = &m_nodeAttr;
		} else if(classAttr->getValue() == "edge") {
			attrMap = &m_edgeAttr;
		} else {
			std::cerr << "ERROR: incorrect attributes tag class "
			          << "(line " << attrsTag.getLine() << ").\n";
			return false;
		}

		if(!readAttrDefs(*attrMap, attrsTag)) {
			return false;
		}
	}

	return true;
}


bool Parser::readNodes(Graph &G, GraphAttributes *GA)
{
	List<XmlTagObject *> nodeTags;
	m_nodesTag->findSonXmlTagObjectByName("node", nodeTags);

	forall_listiterators(XmlTagObject *, it, nodeTags) {
		const XmlTagObject &nodeTag = **it;

		XmlAttributeObject *idAttr;
		nodeTag.findXmlAttributeObjectByName("id", idAttr);
		if(!idAttr) {
			std::cerr << "ERROR: node is missing an attribute "
			          << "(line " << nodeTag.getLine() << ").\n";
			return false;
		}

		const node v = G.newNode();
		m_nodeId[idAttr->getValue()] = v;

		if(GA) {
			readAttributes(*GA, v, nodeTag);
		}
	}

	return true;
}


bool Parser::readCluster(
	Graph &G, ClusterGraph &C, ClusterGraphAttributes *CA, cluster rootCluster,
	const XmlTagObject &rootTag)
{
	List<XmlTagObject *> nodeTags;
	rootTag.findSonXmlTagObjectByName("node", nodeTags);

	forall_listiterators(XmlTagObject *, it, nodeTags) {
		const XmlTagObject &nodeTag = **it;

		XmlAttributeObject *idAttr;
		nodeTag.findXmlAttributeObjectByName("id", idAttr);
		if(!idAttr) {
			std::cerr << "ERROR: node is missing an attribute "
			          << "(line " << nodeTag.getLine() << ").\n";
		}

		// Node is a cluster iff it contains other nodes.
		XmlTagObject *nodesTag;
		nodeTag.findSonXmlTagObjectByName("nodes", nodesTag);
		if(nodesTag) {
			// Node tag found, therefore it is a cluster.
			const cluster c = C.newCluster(rootCluster);
			m_clusterId[idAttr->getValue()] = c;

			if(!readCluster(G, C, CA, c, *nodesTag)) {
				return false;
			}
		} else {
			// Node tag not found, therefore it is "normal" node.
			const node v = G.newNode();
			C.reassignNode(v, rootCluster);
			m_nodeId[idAttr->getValue()] = v;

			if(CA) {
				readAttributes(*CA, v, nodeTag);
			}
		}
	}

	return true;
}


/*
 * Just a helper method to avoid ugly code in Parser#readEdges method. It just
 * populates \a nodes list with either a given \a v node (if not NULL) or all
 * nodes in certain cluster found by performing a lookup with given \a id in
 * \a clusterId association.
 */
static inline bool edgeNodes(
	node v,
	const std::string &id, const HashArray<std::string, cluster> clusterId,
	List<node> &nodes)
{
	if(v) {
		nodes.clear();
		nodes.pushBack(v);
	} else {
		const cluster c = clusterId[id];
		if(!c) {
			return false;
		}

		c->getClusterNodes(nodes);
	}

	return true;
}


bool Parser::readEdges(Graph &G, ClusterGraph *C, GraphAttributes *GA)
{
	List<XmlTagObject *> edgeTags;
	m_edgesTag->findSonXmlTagObjectByName("edge", edgeTags);


	List<node> sourceNodes, targetNodes;

	forall_listiterators(XmlTagObject *, it, edgeTags) {
		const XmlTagObject &edgeTag = **it;

		XmlAttributeObject *sourceAttr, *targetAttr;

		edgeTag.findXmlAttributeObjectByName("source", sourceAttr);
		if(!sourceAttr) {
			std::cerr << "ERROR: edge is missing a source attribute "
			          << "(line " << edgeTag.getLine() << ").\n";
			return false;
		}

		edgeTag.findXmlAttributeObjectByName("target", targetAttr);
		if(!targetAttr) {
			std::cerr << "ERROR: edge is missing a target attribute "
			          << "(line " << edgeTag.getLine() << ").\n";
			return false;
		}

		const std::string &sourceId = sourceAttr->getValue();
		const std::string &targetId = targetAttr->getValue();

		const node source = m_nodeId[sourceId];
		const node target = m_nodeId[targetId];

		if(source && target) {
			const edge e = G.newEdge(source, target);
			if(GA) {
				readAttributes(*GA, e, edgeTag);
			}
		} else if(C && edgeNodes(source, sourceId, m_clusterId, sourceNodes)
		            && edgeNodes(target, targetId, m_clusterId, targetNodes))
		{
			// So, we perform cartesian product on two sets with Graph#newEdge.
			forall_listiterators(node, sit, sourceNodes) {
				forall_listiterators(node, tit, targetNodes) {
					 const edge e = G.newEdge(*sit, *tit);
					 if(GA) {
					 	readAttributes(*GA, e, edgeTag);
					 }
				}
			}
		} else {
			std::cerr << "ERROR: source or target node doesn't exist "
			          << "(line " << edgeTag.getLine() << ").\n";
			return false;
		}
	}

	return true;
}


static inline bool readColor(Color &color, const XmlTagObject &tag)
{
	XmlAttributeObject *redAttr, *greenAttr, *blueAttr, *alphaAttr;
	tag.findXmlAttributeObjectByName("red", redAttr);
	tag.findXmlAttributeObjectByName("green", greenAttr);
	tag.findXmlAttributeObjectByName("blue", blueAttr);
	tag.findXmlAttributeObjectByName("alpha", alphaAttr);

	if(!redAttr || !greenAttr || !blueAttr) {
		std::cerr << "Missing compound attrribute on color tag "
		          << "(line " << tag.getLine() << ").\n";
		return false;
	}

	int compound;

	std::istringstream is;

	is.clear();
	is.str(redAttr->getValue());
	is >> compound;
	color.red(static_cast<__uint8>(compound));

	is.clear();
	is.str(greenAttr->getValue());
	is >> compound;
	color.green(static_cast<__uint8>(compound));

	is.clear();
	is.str(blueAttr->getValue());
	is >> compound;
	color.blue(static_cast<__uint8>(compound));

	if(alphaAttr) {
		is.clear();
		is.str(alphaAttr->getValue());
		is >> compound;
		color.alpha(static_cast<__uint8>(compound));
	}

	return true;
}


static inline bool readVizAttribute(
	GraphAttributes &GA, node v,
	const XmlTagObject &tag)
{
	const long attrs = GA.attributes();

	if(tag.getName() == "viz:position") {
		if(attrs & GraphAttributes::nodeGraphics) {
			XmlAttributeObject *xAttr, *yAttr, *zAttr;
			tag.findXmlAttributeObjectByName("x", xAttr);
			tag.findXmlAttributeObjectByName("y", yAttr);
			tag.findXmlAttributeObjectByName("z", zAttr);

			if(!xAttr || !yAttr) {
				std::cerr << "ERROR: Missing \"x\" or \"y\" on position tag "
				          << "(line " << tag.getLine() << ").\n";
				return false;
			}

			std::istringstream is;

			is.clear();
			is.str(xAttr->getValue());
			is >> GA.x(v);

			is.clear();
			is.str(yAttr->getValue());
			is >> GA.y(v);

			// z attribute is optional and avaliable only in \a threeD mode.
			if(zAttr && (attrs & GraphAttributes::threeD)) {
				is.clear();
				is.str(zAttr->getValue());
				is >> GA.z(v);
			}
		}
	} else if(tag.getName() == "viz:size") {
		XmlAttributeObject *valueAttr;
		tag.findXmlAttributeObjectByName("value", valueAttr);
		if(!valueAttr) {
			std::cerr << "ERROR: \"size\" attribute is missing a value "
			          << "(line " << tag.getLine() << ").\n";
			return false;
		}

		/*
		 * Because size is just a scale here, I assume that all nodes have
		 * some default width and height value. Then, I just rescale them
		 * using given size. Things can go wrong if viz:size is declared twice
		 * for the same node but this is not our problem (just fix this file!).
		 */
		double size;
		std::istringstream is(valueAttr->getValue());
		is >> size;

		GA.width(v) *= size;
		GA.height(v) *= size;
	} else if(tag.getName() == "viz:shape") {
		if(attrs & GraphAttributes::nodeGraphics) {
			XmlAttributeObject *valueAttr;
			tag.findXmlAttributeObjectByName("value", valueAttr);
			if(!valueAttr) {
				std::cerr << "ERROR: \"shape\" attribute is missing a value "
				          << "(line " << tag.getLine() << ").\n";
				return false;
			}

			GA.shape(v) = toShape(valueAttr->getValue());
		}
	} else if(tag.getName() == "viz:color") {
		if(attrs & GraphAttributes::nodeStyle) {
			return readColor(GA.fillColor(v), tag);
		}
	} else {
		std::cerr << "ERROR: Incorrect tag \"" << tag.getName() << "\" "
		          << "(line " << tag.getLine() << ").\n";
		return false;
	}

	return true;
}


static inline bool readVizAttribute(
	GraphAttributes &GA, edge e,
	const XmlTagObject &tag)
{
	const long attrs = GA.attributes();

	if(tag.getName() == "viz:color") {
		if(attrs & GraphAttributes::edgeStyle) {
			return readColor(GA.strokeColor(e), tag);
		}
	} else if(tag.getName() == "viz:thickness") {
		XmlAttributeObject *thickAttr;
		tag.findXmlAttributeObjectByName("value", thickAttr);

		if(!thickAttr) {
			cerr << "ERROR: Missing \"value\" on thickness tag "
			     << "(line " << tag.getLine() << ").\n";
			return false;
		}

		std::istringstream is(thickAttr->getValue());
		if(attrs & GraphAttributes::edgeDoubleWeight) {
			is >> GA.doubleWeight(e);
		} else if(attrs & GraphAttributes::edgeIntWeight) {
			is >> GA.intWeight(e);
		}
	} else if(tag.getName() == "viz:shape") {
		// Values: solid, dotted, dashed, double. Not supported in OGDF.
	} else {
		std::cerr << "ERROR: Incorrect tag \"" << tag.getName() << "\" "
		          << "(line " << tag.getLine() << ").\n";
		return false;
	}

	return true;
}


static inline void readAttValue(
	GraphAttributes &GA, node v,
	const std::string &name, const std::string &value)
{
	const long attrs = GA.attributes();

	// For not "viz" attributes, we use GraphML ones.
	switch(graphml::toAttribute(name)) {
	case graphml::a_nodeType:
		if(attrs & GraphAttributes::nodeType) {
			GA.type(v) = graphml::toNodeType(value);
		}
		break;
	case graphml::a_template:
		if(attrs & GraphAttributes::nodeTemplate) {
			GA.templateNode(v) = value;
		}
		break;
	case graphml::a_nodeWeight:
		if(attrs & GraphAttributes::nodeWeight) {
			std::istringstream ss(value);
			ss >> GA.weight(v);
		}
		break;
	default:
		// Not supported attribute, just ignore.
		break;
	}
}


static inline void readAttValue(
	GraphAttributes &GA, edge e,
	const std::string &name, const std::string &value)
{
	const long attrs = GA.attributes();

	// For not "viz" attributes, we use GraphML ones.
	switch(graphml::toAttribute(name)) {
	case graphml::a_edgeType:
		if(attrs & GraphAttributes::edgeType) {
			GA.type(e) = graphml::toEdgeType(value);
		}
		break;
	case graphml::a_edgeArrow:
		if(attrs & GraphAttributes::edgeArrow) {
			GA.arrowType(e) = graphml::toArrow(value);
		}
		break;
	default:
		// Not supported attribute, just ignore.
		break;
	}
}


template <typename T>
static inline bool readAttValues(
	GraphAttributes &GA, T element,
	const XmlTagObject &tag,
	const HashArray<std::string, std::string> &attrMap)
{
	List<XmlTagObject *> attVals;
	tag.findSonXmlTagObjectByName("attvalue", attVals);

	forall_listiterators(XmlTagObject *, it, attVals) {
		const XmlTagObject &attVal = **it;

		XmlAttributeObject *forAttr, *valueAttr;
		attVal.findXmlAttributeObjectByName("for", forAttr);
		attVal.findXmlAttributeObjectByName("value", valueAttr);

		if(!forAttr || !valueAttr) {
		std::cerr << "ERROR: \"for\" or \"value\" not found for attvalue "
		          << "(line " << attVal.getLine() << ").\n";
			return false;
		}

		const std::string &attrName = attrMap[forAttr->getValue()];
		readAttValue(GA, element, attrName, valueAttr->getValue());
	}

	return true;
}


bool Parser::readAttributes(
	GraphAttributes &GA, node v,
	const XmlTagObject &nodeTag)
{
	for(XmlTagObject *tag = nodeTag.m_pFirstSon; tag; tag = tag->m_pBrother) {
		if(tag->getName() == "nodes") {
			continue;
		} else if(tag->getName() == "attvalues") {
			return readAttValues(GA, v, *tag, m_nodeAttr);
		} else if(!readVizAttribute(GA, v, *tag)) {
			return false;
		}
	}

	return true;
}


bool Parser::readAttributes(
	GraphAttributes &GA, edge e,
	const XmlTagObject &edgeTag)
{
	for(XmlTagObject *tag = edgeTag.m_pFirstSon; tag; tag = tag->m_pBrother) {
		if(tag->getName() == "attvalues") {
			return readAttValues(GA, e, *tag, m_edgeAttr);
		} else if(!readVizAttribute(GA, e, *tag)) {
			return false;
		}
	}

	return true;
}


bool Parser::read(Graph &G)
{
	if(!init()) {
		return false;
	}
	OGDF_ASSERT(m_graphTag);

	G.clear();

	return readNodes(G, NULL) && readEdges(G, NULL, NULL);
}


bool Parser::read(Graph &G, GraphAttributes &GA)
{
	if(!init()) {
		return false;
	}
	OGDF_ASSERT(m_graphTag);

	G.clear();

	return readNodes(G, &GA) && readEdges(G, NULL, &GA);
}


bool Parser::read(Graph &G, ClusterGraph &C)
{
	if(!init()) {
		return false;
	}
	OGDF_ASSERT(m_graphTag);

	G.clear();

	return readCluster(G, C, NULL, C.rootCluster(), *m_nodesTag) &&
	       readEdges(G, &C, NULL);
}


bool Parser::read(Graph &G, ClusterGraph &C, ClusterGraphAttributes &CA)
{
	if(!init()) {
		return false;
	}
	OGDF_ASSERT(m_graphTag);

	G.clear();

	return readCluster(G, C, &CA, C.rootCluster(), *m_nodesTag) &&
	       readEdges(G, &C, &CA);
}


} // end namespace gexf

} // end namespace ogdf
