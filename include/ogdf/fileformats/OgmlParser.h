/*
 * $Revision: 3831 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 10:00:32 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of auxiliary classes OgmlAttributeValue,
 *        OgmlAttribute and OgmlTag.
 *
 * \author Christian Wolf and Bernd Zey
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


#ifdef _MSC_VER
#pragma once
#endif

//KK: Commented out the constraint stuff using //o
//CG: compound graph stuff has been removed with commit 2465

#ifndef OGDF_OGML_PARSER_H
#define OGDF_OGML_PARSER_H

#include <ogdf/fileformats/Ogml.h>
#include <ogdf/fileformats/XmlParser.h>
#include <ogdf/basic/Hashing.h>
#include <ogdf/cluster/ClusterGraph.h>
#include <ogdf/cluster/ClusterGraphAttributes.h>

// constraints
//o#include <ogdf/Constraints.h>


namespace ogdf {

//
// ---------- O g m l P a r s e r ------------------------
//

/**Objects of this class represent a validating parser for files in Ogml.
*/
class OgmlParser
{
private:

	// struct definitions for mapping of templates
	struct OgmlNodeTemplate;
	struct OgmlEdgeTemplate;
	//struct OgmlLabelTemplate;

	struct OgmlSegment;

	class OgmlAttributeValue;
	class OgmlAttribute;
	class OgmlTag;

	friend ostream& operator<<(ostream& os, const OgmlParser::OgmlAttribute& oa);
	friend ostream& operator<<(ostream& os, const OgmlParser::OgmlTag& ot);

	static Hashing<int, OgmlTag>            *s_tags;       //!< Hashtable for saving all ogml tags.
	static Hashing<int, OgmlAttribute>      *s_attributes; //!< Hashtable for saving all ogml attributes.
	static Hashing<int, OgmlAttributeValue> *s_attValues;  //!< Hashtable for saving all values of ogml attributes.


	//! Builds hashtables for tags and attributes.
	static void buildHashTables();

	mutable Ogml::GraphType m_graphType; //!< Saves a graph type. Is set by checkGraphType.

	Hashing<string, const XmlTagObject*> m_ids; //!< Saves all ids of an ogml-file.

	/**
	 * Checks if all tags (XmlTagObject), their attributes (XmlAttributeObject) and
	 * their values are valid (are tags expected, do they own the rigth attributes...)
	 * and sets a valid flag to these. Furthermore it checks if ids of tags are
	 * unique and if id references are valid.
	 * See OgmlTag.h for semantics of the encodings.
	 * Returns the validity state of the current processed tag.
	 */
	int validate(const XmlTagObject *xmlTag, int ogmlTag);

	/**
	 * Wrapper method for validate method above.
	 * Returns true when validation is successfull, false otherwise.
	 */
	//bool validate(const char* fileName);

	//! Prints some useful information about un-/successful validation.
	void printValidityInfo(const OgmlTag &ot,
		const XmlTagObject &xto,
		int valStatus,
		int line);

	/**
	 * Finds the OGML-tag in the parse tree with the specified id,
	 * stores the tag in xmlTag
	 * recTag is the tag for recursive calls
	 * returns false if something goes wrong
	 */
	//bool getXmlTagObjectById(XmlTagObject *recTag, string id, XmlTagObject *&xmlTag);

	/**
	 * Checks the graph type and stores it in the member variable m_graphType
	 * xmlTag has to be the root or the graph or the structure Ogml-tag
	 * returns false if something goes wrong
	 */
	bool checkGraphType(const XmlTagObject *xmlTag) const;

	//! Returns true iff subgraph is an hierarchical graph.
	bool isGraphHierarchical(const XmlTagObject *xmlTag) const;

	//! Returns true iff node contains other nodes.
	bool isNodeHierarchical(const XmlTagObject *xmlTag) const;

	Ogml::GraphType getGraphType() { return m_graphType; };


	// id hash tables
	// required variables for building
	// hash table with id from file and node
	Hashing<string, node> m_nodes;
	Hashing<string, edge> m_edges;
	Hashing<string, cluster> m_clusters;
	// hash table for bend-points
	Hashing<string, DPoint> m_points;

	// hash table for checking uniqueness of ids
	// (key:) int = id in the created graph
	// (info:) string = id in the ogml file
	Hashing<int, string> m_nodeIds;
	Hashing<int, string> m_edgeIds;
	Hashing<int, string> m_clusterIds;

	// build methods

	//! Builds a graph; ignores nodes which have hierarchical structure.
	bool buildGraph(Graph &G);

	//! Builds a cluster graph.
	bool buildCluster(
		const XmlTagObject *rootTag,
		Graph &G,
		ClusterGraph &CG);

	//! Recursive part of buildCluster.
	bool buildClusterRecursive(
		const XmlTagObject *xmlTag,
		cluster parent,
		Graph &G,
		ClusterGraph &CG);

	//! Build a cluster graph with style/layout attributes.
	bool addAttributes(
		Graph &G,
		GraphAttributes &GA,
		ClusterGraphAttributes *pCGA,
		const XmlTagObject *root);

	//! Recursive method for setting labels of clusters and nodes.
	bool setLabelsRecursive(
		Graph &G,
		GraphAttributes &GA,
		ClusterGraphAttributes *pCGA,
		XmlTagObject *root);

	// helping pointer for constraints-loading
	// this pointer is set in the building methods
	// so we don't have to traverse the tree in buildConstraints
	XmlTagObject* m_constraintsTag;

	// hashing lists for templates
	//  string = id
	Hashing<string, OgmlNodeTemplate*> m_ogmlNodeTemplates;
	Hashing<string, OgmlEdgeTemplate*> m_ogmlEdgeTemplates;
	//Hashing<string, OgmlLabelTemplate> m_ogmlLabelTemplates;

	// auxiliary methods for mapping graph attributes

	//! Returns fill pattern of string \a s.
	FillPattern getFillPattern(string s);

	//! Returns the shape as an integer value.
	Shape getShape(string s);

	//! Maps the OGML attribute values to corresponding GDE values.
	string getNodeTemplateFromOgmlValue(string s);

	//! Returns the line type as an integer value.
	StrokeType getStrokeType(string s);

	// arrow style, actually a "boolean" function
	// because it returns only 0 or 1 according to GDE
	// sot <=> source or target
	int getArrowStyleAsInt(string s);

	// the matching method to getArrowStyleAsInt
	EdgeArrow getArrowStyle(int i);

	// function that operates on a string
	// the input string contains "&lt;" instead of "<"
	//  and "&gt;" instead of ">"
	//  to disable interpreting the string as xml-tags (by XmlParser)
	// so this function substitutes  "<" for "&lt;"
	string getLabelCaptionFromString(string str);

	//! Returns the integer value of the id at the end of the string (if it exists).
	bool getIdFromString(string str, int &id);

	//! Validiation method.
	void validate(istream &is);


	//! Unified read method for graphs.
	bool doRead(
		istream &is,
		Graph &G,
		ClusterGraph *pCG,
		GraphAttributes *pGA,
		ClusterGraphAttributes *pCGA);

public:

	//! Constructs an OGML parser.
	OgmlParser() { }

	~OgmlParser() { }


	//! Reads a graph \a G from file \a fileName in OGML format.
	/**
	 * @param is is the input stream to be parsed as OGML file.
	 * @param G is the graph to be build from the OGML file.
	 * @return true if succesfull, false otherwise.
	 */
	bool read(istream &is, Graph &G) {
		return doRead(is, G, 0, 0, 0);
	}

	//! Reads a cluster graph \a CG from file \a fileName in OGML format.
	/**
	 * @param is is the input stream to be parsed as OGML file.
	 * @param G is the graph to be build from the OGML file; must be the graph associated with \a CG.
	 * @param CG is the cluster graph to be build from the OGML file.
	 * @return true if succesfull, false otherwise.
	 */
	bool read(istream &is, Graph &G, ClusterGraph &CG) {
		return doRead(is, G, &CG, 0, 0);
	}

	//! Reads a cluster graph \a CG with attributes \a CGA from file \a fileName in OGML format.
	/**
	 * @param is is the input stream to be parsed as OGML file.
	 * @param G is the graph to be build from the OGML file.
	 * @param GA are the graph attributes (associated with \a G) in which layout and style information are stored.
	 * @return true if succesfull, false otherwise.
	 */
	bool read(
		istream &is,
		Graph &G,
		GraphAttributes &GA)
	{
		return doRead(is, G, 0, &GA, 0);
	}

	//! Reads a cluster graph \a CG with attributes \a CGA from file \a fileName in OGML format.
	/**
	 * @param is is the input stream to be parsed as OGML file.
	 * @param G is the graph to be build from the OGML file; must be the graph associated with \a CG.
	 * @param CG is the cluster graph to be build from the OGML file.
	 * @param CGA are the cluster graph attributes (associated with \a CG) in which layout and style information are stored.
	 * @return true if succesfull, false otherwise.
	 */
	bool read(
		istream &is,
		Graph &G,
		ClusterGraph &CG,
		ClusterGraphAttributes &CGA)
	{
		return doRead(is, G, &CG, &CGA, &CGA);
	}
};//end class OGMLParser

}//end namespace ogdf

#endif

