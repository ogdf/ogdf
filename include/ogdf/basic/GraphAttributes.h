/** \file
 * \brief Declaration of class GraphAttributes which extends a Graph
 *        by additional attributes.
 *
 * \author Carsten Gutwenger
 *         Karsten Klein
 *         Joachim Kupke
 *         Sebastian Leipert
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

#pragma once

#include <ogdf/basic/NodeArray.h>
#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/geometry.h>
#include <ogdf/basic/LayoutStandards.h>
#include <ogdf/basic/exceptions.h>


namespace ogdf {

//---------------------------------------------------------
// GraphAttributes
// graph topology + graphical attributes
//---------------------------------------------------------
//! Stores additional attributes of a graph (like layout information).
/**
 * @ingroup graph-drawing
 *
 * It is frequently necessary to associate additional attributes with a graph.
 * The class GraphAttributes provides various such attributes and is the
 * central place were such attributes are stored.
 *
 * Attributes are simply stored in node or edge arrays; for memory consumption
 * reasons, only a subset of these arrays is in fact initialized for the graph;
 * non-initialized arrays require only a few bytes of extra memory.
 *
 * Which arrays are initialized is specified by a bit vector; each bit in this
 * bit vector corresponds to one or more attributes. E.g., \a #nodeGraphics
 * corresponds to the attributes \a #m_x, \a #m_y, \a #m_width, and \a #m_height;
 * whereas \a #edgeDoubleWeight only corresponds to the attribute \a #m_doubleWeight.
 *
 * Attributes can be initialized by the constructor GraphAttributes(const Graph &,long)
 * or the function initAttributes(); attributes can also be deinitialized by
 * calling destroyAttributes().
 */

class OGDF_EXPORT GraphAttributes {

protected:
	const Graph *m_pGraph; //!< associated graph

	bool m_directed; //!< whether or not the graph is directed

	// graphical representation of nodes
	NodeArray<double>       m_x;				//!< x-coordinate of a node
	NodeArray<double>       m_y;				//!< y-coordinate of a node
	NodeArray<double>       m_z;				//!< z-coordinate of a node
	NodeArray<double>       m_nodeLabelPosX;		//!< x-coordinate of a node label
	NodeArray<double>       m_nodeLabelPosY;		//!< y-coordinate of a node label
	NodeArray<double>       m_nodeLabelPosZ;		//!< z-coordinate of a node label
	NodeArray<double>       m_width;			//!< width of a node's bounding box
	NodeArray<double>       m_height;			//!< height of a nodes's bounding box
	NodeArray<Shape>        m_nodeShape;		//!< shape of a node
	NodeArray<string>       m_nodeLabel;		//!< label of a node
	NodeArray<Stroke>       m_nodeStroke;		//!< stroke of a node
	NodeArray<Fill>         m_nodeFill;			//!< fill of a node
	NodeArray<string>       m_nodeTemplate;		//!< name of template of a node

	// other node attributes
	NodeArray<int>             m_nodeId;		//!< user ID of a node
	NodeArray<int>             m_nodeIntWeight;	//!< (integer) weight of a node
	NodeArray<Graph::NodeType> m_vType;			//!< type (vertex, dummy, generalizationMerger)

	// graphical representation of edges
	EdgeArray<DPolyline>       m_bends;			//!< list of bend points of an edge
	EdgeArray<string>          m_edgeLabel;		//!< label of an edge
	EdgeArray<EdgeArrow>       m_edgeArrow;		//!< arrow type of an edge
	EdgeArray<Stroke>          m_edgeStroke;	//!< stroke of an edge

	// other edge attributes
	EdgeArray<int>             m_intWeight;		//!< (integer) weight of an edge
	EdgeArray<double>          m_doubleWeight;	//!< (real number) weight of an edge
	EdgeArray<Graph::EdgeType> m_eType;			//!< type of an edge (association or generalization)
	EdgeArray<uint32_t>        m_subGraph;		//!< is element of subgraphs given by bitvector

	long m_attributes;	//!< bit vector of currently used attributes

public:
	//! Bits for specifying attributes.
	enum {
		nodeGraphics     = 0x00001, //!< node attributes m_x, m_y, m_width, m_height, m_nodeShape
		edgeGraphics     = 0x00002, //!< edge attribute  m_bends
		edgeIntWeight    = 0x00004, //!< edge attribute  m_intWeight
		edgeDoubleWeight = 0x00008, //!< edge attribute  m_doubleWeight
		edgeLabel        = 0x00010, //!< edge attribute  m_edgeLabel
		nodeLabel        = 0x00020, //!< node attribute  m_nodeLabel
		edgeType         = 0x00040, //!< edge attribute  m_eType
		nodeType         = 0x00080, //!< node attribute  m_vType
		nodeId           = 0x00100, //!< node attribute  m_nodeId
		edgeArrow        = 0x00200, //!< edge attribute  m_edgeArrow
		edgeStyle        = 0x00400, //!< edge attribute  m_edgeStroke
		nodeStyle        = 0x00800, //!< node attributes m_nodeStroke, m_nodeFill
		nodeTemplate     = 0x01000, //!< node attribute  m_nodeTemplate
		edgeSubGraphs    = 0x02000, //!< edge attribute  m_subGraph
		nodeWeight       = 0x04000, //!< node attribute  m_nodeIntWeight
		threeD           = 0x08000, //!< node attribute  m_z, note that all methods work on 2D coordinates only
		nodeLabelPosition = 0x10000, //!< node attribute  m_nodeLabelPosX,Y,Z
	};

	/**
	 * @name Construction and management of attributes
	 */
	//@{

	//! Constructs graph attributes for no associated graph (default constructor).
	/**
	 * The associated graph can be set later with the init() function.
	 */
	GraphAttributes();

	//! Constructs graph attributes associated with the graph \a G.
	/**
	 * @param G is the associated graph.
	 * @param initAttributes specifies the set of attributes that can be accessed.
	 */
	explicit GraphAttributes(const Graph &G, long initAttributes = nodeGraphics | edgeGraphics);

	virtual ~GraphAttributes() {
	}

	//! Returns currently accessible attributes.
	long attributes() const {
		return m_attributes;
	}

	//! Returns true iff all attributes in \a attr are available.
	inline bool has(long attr) const {
		return (m_attributes & attr) == attr;
	}

	//! Initializes the graph attributes for graph \a G.
	/**
	 * @param G is the new associated graph.
	 * @param initAttr specifies the set of attributes that can be accessed.
	 *
	 * \warning All attributes that were allocated before are destroyed by this function!
	 *  If you wish to extend the set of allocated attributes, use initAttributes().
	 */
	virtual void init(const Graph &G, long initAttr);

	//! Initializes attributes in \a attr for usage.
	void initAttributes(long attr);

	//! Destroys attributes in attr.
	void destroyAttributes(long attr);

	//! Returns a reference to the associated graph.
	const Graph& constGraph() const {
		return *m_pGraph;
	}

	//@}
	/**
	 * @name General attributes
	 */
	//@{

	//! Returns if the graph is directed.
	bool directed() const {
		return m_directed;
	}

	//! Sets if the graph is directed to \a directed.
	void setDirected(bool directed) {
		m_directed = directed;
	}

	//@}
	/**
	 * @name Node attributes
	 */
	//@{


	//! Returns the x-coordinate of node \a v.
	/**
	 * \pre \a nodeGraphics is enabled
	 */
	double x(node v) const {
#ifdef OGDF_DEBUG
		if(!has(nodeGraphics)){
			throw PreconditionViolatedException();
		}
#endif
		return m_x[v];
	}


	//! Returns the x-coordinate of node \a v.
	/**
	 * \pre \a nodeGraphics is enabled
	 */
	double &x(node v) {
#ifdef OGDF_DEBUG
		if(!has(nodeGraphics)){
			throw PreconditionViolatedException();
		}
#endif
		return m_x[v];
	}


	//! Returns the y-coordinate of node \a v.
	/**
	 * \pre \a nodeGraphics is enabled
	 */
	double y(node v) const {
#ifdef OGDF_DEBUG
		if(!has(nodeGraphics)){
			throw PreconditionViolatedException();
		}
#endif
		return m_y[v];
	}


	//! Returns the y-coordinate of node \a v.
	/**
	 * \pre \a nodeGraphics is enabled
	 */
	double &y(node v) {
#ifdef OGDF_DEBUG
		if(!has(nodeGraphics)){
			throw PreconditionViolatedException();
		}
#endif
		return m_y[v];
	}


	//! Returns the z-coordinate of node \a v.
	/**
	 * \pre \a threeD is enabled
	 */
	double z(node v) const {
#ifdef OGDF_DEBUG
		if(!has(threeD)){
			throw PreconditionViolatedException();
		}
#endif
		return m_z[v];
	}


	//! Returns the z-coordinate of node \a v.
	/**
	 * \pre \a threeD is enabled
	 */
	double &z(node v) {
#ifdef OGDF_DEBUG
		if(!has(threeD)){
			throw PreconditionViolatedException();
		}
#endif
		return m_z[v];
	}


	//! Returns the label x-coordinate of node \a v.
	/**
	 * \pre \a nodeLabelPosition is enabled
	 */
	double xLabel(node v) const {
#ifdef OGDF_DEBUG
		if(!has(nodeLabelPosition)){
			throw PreconditionViolatedException();
		}
#endif
		return m_nodeLabelPosX[v];
	}

	//! Returns the label x-coordinate of node \a v.
	/**
	 * \pre \a nodeLabelPosition is enabled
	 */
	double &xLabel(node v) {
#ifdef OGDF_DEBUG
		if(!has(nodeLabelPosition)){
			throw PreconditionViolatedException();
		}
#endif
		return m_nodeLabelPosX[v];
	}


	//! Returns the label y-coordinate of node \a v.
	/**
	 * \pre \a nodeLabelPosition is enabled
	 */
	double yLabel(node v) const {
#ifdef OGDF_DEBUG
		if(!has(nodeLabelPosition)){
			throw PreconditionViolatedException();
		}
#endif
		return m_nodeLabelPosY[v];
	}

	//! Returns the label y-coordinate of node \a v.
	/**
	 * \pre \a nodeLabelPosition is enabled
	 */
	double &yLabel(node v) {
#ifdef OGDF_DEBUG
		if(!has(nodeLabelPosition)){
			throw PreconditionViolatedException();
		}
#endif
		return m_nodeLabelPosY[v];
	}


	//! Returns the label z-coordinate of node \a v.
	/**
	 * \pre \a nodeLabelPosition and \a threeD are enabled
	 */
	double zLabel(node v) const {
#ifdef OGDF_DEBUG
		if(!has(nodeLabelPosition | threeD)){
			throw PreconditionViolatedException();
		}
#endif
		return m_nodeLabelPosZ[v];
	}

	//! Returns the label z-coordinate of node \a v.
	/**
	 * \pre \a nodeLabelPosition and \a threeD are enabled
	 */
	double &zLabel(node v) {
#ifdef OGDF_DEBUG
		if(!has(nodeLabelPosition | threeD)){
			throw PreconditionViolatedException();
		}
#endif
		return m_nodeLabelPosZ[v];
	}


	//! Returns the width of the bounding box of node \a v.
	/**
	 * \pre \a nodeGraphics is enabled
	 */
	double width(node v) const {
#ifdef OGDF_DEBUG
		if(!has(nodeGraphics)){
			throw PreconditionViolatedException();
		}
#endif
		return m_width[v];
	}

	//! Returns the width of the bounding box of node \a v.
	/**
	 * \pre \a nodeGraphics is enabled
	 */
	double &width(node v) {
#ifdef OGDF_DEBUG
		if(!has(nodeGraphics)){
			throw PreconditionViolatedException();
		}
#endif
		return m_width[v];
	}


	//! Returns a reference to the node array \a m_width.
	/**
	 * \pre \a nodeGraphics is enabled
	 */
	const NodeArray<double> &width() const {
#ifdef OGDF_DEBUG
		if(!has(nodeGraphics)){
			throw PreconditionViolatedException();
		}
#endif
		return m_width;
	}

	//! Returns a reference to the node array \a m_width.
	/**
	 * \pre \a nodeGraphics is enabled
	 */
	NodeArray<double> &width() {
#ifdef OGDF_DEBUG
		if(!has(nodeGraphics)){
			throw PreconditionViolatedException();
		}
#endif
		return m_width;
	}


	//! Returns the height of the bounding box of node \a v.
	/**
	 * \pre \a nodeGraphics is enabled
	 */
	double height(node v) const {
#ifdef OGDF_DEBUG
		if(!has(nodeGraphics)){
			throw PreconditionViolatedException();
		}
#endif
		return m_height[v];
	}

	//! Returns the height of the bounding box of node \a v.
	/**
	 * \pre \a nodeGraphics is enabled
	 */
	double &height(node v) {
#ifdef OGDF_DEBUG
		if(!has(nodeGraphics)){
			throw PreconditionViolatedException();
		}
#endif
		return m_height[v];
	}


	//! Returns a reference to the node array \a m_height.
	/**
	 * \pre \a nodeGraphics is enabled
	 */
	const NodeArray<double> &height() const {
#ifdef OGDF_DEBUG
		if(!has(nodeGraphics)){
			throw PreconditionViolatedException();
		}
#endif
		return m_height;
	}

	//! Returns a reference to the node array \a m_height.
	/**
	 * \pre \a nodeGraphics is enabled
	 */
	NodeArray<double> &height() {
#ifdef OGDF_DEBUG
		if(!has(nodeGraphics)){
			throw PreconditionViolatedException();
		}
#endif
		return m_height;
	}


	//! Returns the shape type of node \a v.
	/**
	 * \pre \a nodeGraphics is enabled
	 */
	Shape shape(node v) const {
#ifdef OGDF_DEBUG
		if(!has(nodeGraphics)){
			throw PreconditionViolatedException();
		}
#endif
		return m_nodeShape[v];
	}

	//! Returns the shape type of node \a v.
	/**
	 * \pre \a nodeGraphics is enabled
	 */
	Shape &shape(node v) {
#ifdef OGDF_DEBUG
		if(!has(nodeGraphics)){
			throw PreconditionViolatedException();
		}
#endif
		return m_nodeShape[v];
	}


	//! Returns the stroke type of node \a v.
	/**
	 * \pre \a nodeStyle is enabled
	 */
	StrokeType strokeType(node v) const {
#ifdef OGDF_DEBUG
		if(!has(nodeStyle)){
			throw PreconditionViolatedException();
		}
#endif
		return m_nodeStroke[v].m_type;
	}

	//! Sets the stroke type of node \a v to \a st.
	/**
	 * \pre \a nodeStyle is enabled
	 */
	void setStrokeType(node v, StrokeType st) {
#ifdef OGDF_DEBUG
		if(!has(nodeStyle)){
			throw PreconditionViolatedException();
		}
#endif
		m_nodeStroke[v].m_type = st;
	}


	//! Returns the stroke color of node \a v.
	/**
	 * \pre \a nodeStyle is enabled
	 */
	const Color &strokeColor(node v) const {
#ifdef OGDF_DEBUG
		if(!has(nodeStyle)){
			throw PreconditionViolatedException();
		}
#endif
		return m_nodeStroke[v].m_color;
	}

	//! Returns the stroke color of node \a v.
	/**
	 * \pre \a nodeStyle is enabled
	 */
	Color &strokeColor(node v) {
#ifdef OGDF_DEBUG
		if(!has(nodeStyle)){
			throw PreconditionViolatedException();
		}
#endif
		return m_nodeStroke[v].m_color;
	}


	//! Returns the stroke width of node \a v.
	/**
	 * \pre \a nodeStyle is enabled
	 */
	float strokeWidth(node v) const {
#ifdef OGDF_DEBUG
		if(!has(nodeStyle)){
			throw PreconditionViolatedException();
		}
#endif
		return m_nodeStroke[v].m_width;
	}

	//! Returns the stroke width of node \a v.
	/**
	 * \pre \a nodeStyle is enabled
	 */
	float &strokeWidth(node v) {
#ifdef OGDF_DEBUG
		if(!has(nodeStyle)){
			throw PreconditionViolatedException();
		}
#endif
		return m_nodeStroke[v].m_width;
	}


	//! Returns the fill pattern of node \a v.
	/**
	 * \pre \a nodeStyle is enabled
	 */
	FillPattern fillPattern(node v) const {
#ifdef OGDF_DEBUG
		if(!has(nodeStyle)){
			throw PreconditionViolatedException();
		}
#endif
		return m_nodeFill[v].m_pattern;
	}

	//! Sets the fill pattern of node \a v to \a fp.
	/**
	 * \pre \a nodeStyle is enabled
	 */
	void setFillPattern(node v, FillPattern fp) {
#ifdef OGDF_DEBUG
		if(!has(nodeStyle)){
			throw PreconditionViolatedException();
		}
#endif
		m_nodeFill[v].m_pattern = fp;
	}


	//! Returns the fill color of node \a v.
	/**
	 * \pre \a nodeStyle is enabled
	 */
	const Color &fillColor(node v) const {
#ifdef OGDF_DEBUG
		if(!has(nodeStyle)){
			throw PreconditionViolatedException();
		}
#endif
		return m_nodeFill[v].m_color;
	}

	//! Returns the fill color of node \a v.
	/**
	 * \pre \a nodeStyle is enabled
	 */
	Color &fillColor(node v) {
#ifdef OGDF_DEBUG
		if(!has(nodeStyle)){
			throw PreconditionViolatedException();
		}
#endif
		return m_nodeFill[v].m_color;
	}


	//! Returns the background color of fill patterns for node \a v.
	/**
	 * \pre \a nodeStyle is enabled
	 */
	const Color &fillBgColor(node v) const {
#ifdef OGDF_DEBUG
		if(!has(nodeStyle)){
			throw PreconditionViolatedException();
		}
#endif
		return m_nodeFill[v].m_bgColor;
	}

	//! Returns the background color of fill patterns for node \a v.
	/**
	 * \pre \a nodeStyle is enabled
	 */
	Color &fillBgColor(node v) {
#ifdef OGDF_DEBUG
		if(!has(nodeStyle)){
			throw PreconditionViolatedException();
		}
#endif
		return m_nodeFill[v].m_bgColor;
	}


	//! Returns the label of node \a v.
	/**
	 * \pre \a nodeLabel is enabled
	 */
	const string &label(node v) const {
#ifdef OGDF_DEBUG
		if(!has(nodeLabel)){
			throw PreconditionViolatedException();
		}
#endif
		return m_nodeLabel[v];
	}

	//! Returns the label of node \a v.
	/**
	 * \pre \a nodeLabel is enabled
	 */
	string &label(node v) {
#ifdef OGDF_DEBUG
		if(!has(nodeLabel)){
			throw PreconditionViolatedException();
		}
#endif
		return m_nodeLabel[v];
	}


	//! Returns the template name of node \a v.
	/**
	 * \pre \a nodeTemplate is enabled
	 */
	const string &templateNode(node v) const {
#ifdef OGDF_DEBUG
		if(!has(nodeTemplate)){
			throw PreconditionViolatedException();
		}
#endif
		return m_nodeTemplate[v];
	}

	//! Returns the template name of node \a v.
	/**
	 * \pre \a nodeTemplate is enabled
	 */
	string &templateNode(node v) {
#ifdef OGDF_DEBUG
		if(!has(nodeTemplate)){
			throw PreconditionViolatedException();
		}
#endif
		return m_nodeTemplate[v];
	}


	//! Returns the weight of node \a v.
	/**
	 * \pre \a nodeWeight is enabled
	 */
	int weight(node v) const {
#ifdef OGDF_DEBUG
		if(!has(nodeWeight)){
			throw PreconditionViolatedException();
		}
#endif
		return m_nodeIntWeight[v];
	}

	//! Returns the weight of node \a v.
	/**
	 * \pre \a nodeWeight is enabled
	 */
	int &weight(node v) {
#ifdef OGDF_DEBUG
		if(!has(nodeWeight)){
			throw PreconditionViolatedException();
		}
#endif
		return m_nodeIntWeight[v];
	}


	//! Returns the type of node \a v.
	/**
	 * \pre \a nodeType is enabled
	 */
	Graph::NodeType type(node v) const {
#ifdef OGDF_DEBUG
		if(!has(nodeType)){
			throw PreconditionViolatedException();
		}
#endif
		return m_vType.valid() ? m_vType[v] : Graph::vertex;
	}

	//! Returns the type of node \a v.
	/**
	 * \pre \a nodeType is enabled
	 */
	Graph::NodeType &type(node v) {
#ifdef OGDF_DEBUG
		if(!has(nodeType)){
			throw PreconditionViolatedException();
		}
#endif
		return m_vType[v];
	}


	//! Returns the user ID of node \a v.
	/**
	 * \pre \a nodeId is enabled
	 */
	int idNode(node v) const {
#ifdef OGDF_DEBUG
		if(!has(nodeId)){
			throw PreconditionViolatedException();
		}
#endif
		return m_nodeId[v];
	}

	//! Returns the user ID of node \a v.
	/**
	 * \pre \a nodeId is enabled
	 */
	int &idNode(node v) {
#ifdef OGDF_DEBUG
		if(!has(nodeId)){
			throw PreconditionViolatedException();
		}
#endif
		return m_nodeId[v];
	}

	//@}
	/**
	 * @name Edge attributes
	 */
	//@{


	//! Returns the list of bend points of edge \a e.
	/**
	 * \pre \a edgeGraphics is enabled
	 */
	const DPolyline &bends(edge e) const {
#ifdef OGDF_DEBUG
		if(!has(edgeGraphics)){
			throw PreconditionViolatedException();
		}
#endif
		return m_bends[e];
	}

	//! Returns the list of bend points of edge \a e.
	/**
	 * \pre \a edgeGraphics is enabled
	 */
	DPolyline &bends(edge e) {
#ifdef OGDF_DEBUG
		if(!has(edgeGraphics)){
			throw PreconditionViolatedException();
		}
#endif
		return m_bends[e];
	}


	//! Returns the arrow type of edge \a e.
	/**
	 * \pre \a edgeArrow is enabled
	 */
	EdgeArrow arrowType(edge e) const {
#ifdef OGDF_DEBUG
		if(!has(edgeArrow)){
			throw PreconditionViolatedException();
		}
#endif
		return m_edgeArrow[e];
	}

	//! Returns the arrow type of edge \a e.
	/**
	 * \pre \a edgeArrow is enabled
	 */
	EdgeArrow &arrowType(edge e) {
#ifdef OGDF_DEBUG
		if(!has(edgeArrow)){
			throw PreconditionViolatedException();
		}
#endif
		return m_edgeArrow[e];
	}


	//! Returns the stroke type of edge \a e.
	/**
	 * \pre \a edgeStyle is enabled
	 */
	StrokeType strokeType(edge e) const {
#ifdef OGDF_DEBUG
		if(!has(edgeStyle)){
			throw PreconditionViolatedException();
		}
#endif
		return m_edgeStroke[e].m_type;
	}

	//! Sets the stroke type of edge \a e to \a st.
	/**
	 * \pre \a edgeStyle is enabled
	 */
	void setStrokeType(edge e, StrokeType st) {
#ifdef OGDF_DEBUG
		if(!has(edgeStyle)){
			throw PreconditionViolatedException();
		}
#endif
		m_edgeStroke[e].m_type = st;
	}


	//! Returns the stroke color of edge \a e.
	/**
	 * \pre \a edgeStyle is enabled
	 */
	const Color &strokeColor(edge e) const {
#ifdef OGDF_DEBUG
		if(!has(edgeStyle)){
			throw PreconditionViolatedException();
		}
#endif
		return m_edgeStroke[e].m_color;
	}

	//! Returns the stroke color of edge \a e.
	/**
	 * \pre \a edgeStyle is enabled
	 */
	Color &strokeColor(edge e) {
#ifdef OGDF_DEBUG
		if(!has(edgeStyle)){
			throw PreconditionViolatedException();
		}
#endif
		return m_edgeStroke[e].m_color;
	}


	//! Returns the stroke width of edge \a e.
	/**
	 * \pre \a edgeStyle is enabled
	 */
	float strokeWidth(edge e) const {
#ifdef OGDF_DEBUG
		if(!has(edgeStyle)){
			throw PreconditionViolatedException();
		}
#endif
		return m_edgeStroke[e].m_width;
	}

	//! Returns the stroke width of edge \a e.
	/**
	 * \pre \a edgeStyle is enabled
	 */
	float &strokeWidth(edge e) {
#ifdef OGDF_DEBUG
		if(!has(edgeStyle)){
			throw PreconditionViolatedException();
		}
#endif
		return m_edgeStroke[e].m_width;
	}


	//! Returns the label of edge \a e.
	/**
	 * \pre \a edgeLabel is enabled
	 */
	const string &label(edge e) const {
#ifdef OGDF_DEBUG
		if(!has(edgeLabel)){
			throw PreconditionViolatedException();
		}
#endif
		return m_edgeLabel[e];
	}

	//! Returns the label of edge \a e.
	/**
	 * \pre \a edgeLabel is enabled
	 */
	string &label(edge e) {
#ifdef OGDF_DEBUG
		if(!has(edgeLabel)){
			throw PreconditionViolatedException();
		}
#endif
		return m_edgeLabel[e];
	}


	//! Returns the (integer) weight of edge \a e.
	/**
	 * \pre \a edgeIntWeight is enabled
	 */
	int intWeight(edge e) const {
#ifdef OGDF_DEBUG
		if(!has(edgeIntWeight)){
			throw PreconditionViolatedException();
		}
#endif
		return m_intWeight[e];
	}

	//! Returns the (integer) weight of edge \a e.
	/**
	 * \pre \a edgeIntWeight is enabled
	 */
	int &intWeight(edge e) {
#ifdef OGDF_DEBUG
		if(!has(edgeIntWeight)){
			throw PreconditionViolatedException();
		}
#endif
		return m_intWeight[e];
	}


	//! Returns the (real number) weight of edge \a e.
	/**
	 * \pre \a edgeDoubleWeight is enabled
	 */
	double doubleWeight(edge e) const {
#ifdef OGDF_DEBUG
		if(!has(edgeDoubleWeight)){
			throw PreconditionViolatedException();
		}
#endif
		return m_doubleWeight[e];
	}

	//! Returns the (real number) weight of edge \a e.
	/**
	 * \pre \a edgeDoubleWeight is enabled
	 */
	double &doubleWeight(edge e) {
#ifdef OGDF_DEBUG
		if(!has(edgeDoubleWeight)){
			throw PreconditionViolatedException();
		}
#endif
		return m_doubleWeight[e];
	}


	//! Returns the type of edge \a e.
	/**
	 * \pre \a edgeType is enabled
	 */
	Graph::EdgeType type(edge e) const {
#ifdef OGDF_DEBUG
		if(!has(edgeType)){
			throw PreconditionViolatedException();
		}
#endif
		return m_eType.valid() ? m_eType[e] : Graph::association;
	}

	//! Returns the type of edge \a e.
	/**
	 * \pre \a edgeType is enabled
	 */
	Graph::EdgeType &type(edge e) {
#ifdef OGDF_DEBUG
		if(!has(edgeType)){
			throw PreconditionViolatedException();
		}
#endif
		return m_eType[e];
	}


	//! Returns the edgesubgraph value of an edge \a e.
	/**
	 * \pre \a edgeSubGraphs is enabled
	 */
	uint32_t subGraphBits(edge e) const {
#ifdef OGDF_DEBUG
		if(!has(edgeSubGraphs)){
			throw PreconditionViolatedException();
		}
#endif
		return m_subGraph[e];
	}

	//! Returns the edgesubgraph value of an edge \a e.
	/**
	 * \pre \a edgeSubGraphs is enabled
	 */
	uint32_t &subGraphBits(edge e) {
#ifdef OGDF_DEBUG
		if(!has(edgeSubGraphs)){
			throw PreconditionViolatedException();
		}
#endif
		return m_subGraph[e];
	}


	//! Checks whether edge \a e belongs to basic graph \a n.
	/**
	 * \pre \a edgeSubGraphs is enabled
	 */
	bool inSubGraph(edge e, int n) const {
#ifdef OGDF_DEBUG
		if(!has(edgeSubGraphs)){
			throw PreconditionViolatedException();
		}
#endif
		OGDF_ASSERT(n>=0);
		OGDF_ASSERT(n<32);
		return (m_subGraph[e] & (1 << n)) != 0;
	}


	//! Adds edge \a e to basic graph \a n.
	/**
	 * \pre \a edgeSubGraphs is enabled
	 */
	void addSubGraph(edge e, int n) {
#ifdef OGDF_DEBUG
		if(!has(edgeSubGraphs)){
			throw PreconditionViolatedException();
		}
#endif
		OGDF_ASSERT(n>=0);
		OGDF_ASSERT(n<32);
		m_subGraph[e] |= (1 << n);
	}


	//! Removes edge \a e from basic graph \a n.
	/**
	 * \pre \a edgeSubGraphs is enabled
	 */
	void removeSubGraph(edge e, int n) {
#ifdef OGDF_DEBUG
		if(!has(edgeSubGraphs)){
			throw PreconditionViolatedException();
		}
#endif
		OGDF_ASSERT(n>=0);
		OGDF_ASSERT(n<32);
		m_subGraph[e] &= ~(1 << n);
	}

	//@}
	/**
	* @name Layout transformations
	*/
	//@{

	//! Scales the layout by (\a sx,\a sy).
	/**
	 * If \a scaleNodes is true, node sizes are scaled as well.
	 *
	 * \param sx         is the scaling factor for x-coordinates.
	 * \param sy         is the scaling factor for y-coordinates.
	 * \param scaleNodes determines if nodes size are scaled as well (true) or not.
	 *
	 * \pre \a nodeGraphics and \a edgeGraphics are enabled
	 */
	virtual void scale(double sx, double sy, bool scaleNodes = true);

	//! Scales the layout by \a s.
	/**
	* If \a scaleNodes is true, node sizes are scaled as well.
	*
	* \param s          is the scaling factor for both x- and y-coordinates.
	* \param scaleNodes determines if nodes size are scaled as well (true) or not.
	*
	* \pre \a nodeGraphics and \a edgeGraphics are enabled
	*/
	virtual void scale(double s, bool scaleNodes = true) { scale(s, s, scaleNodes); }

	//! Translates the layout by (\a dx,\a dy).
	/**
	 * \param dx is the translation in x-direction.
	 * \param dy is the translation in y-direction.
	 */
	virtual void translate(double dx, double dy);

	//! Translates the layout such that the lower left corner is at (0,0).
	virtual void translateToNonNeg();

	//! Flips the layout vertically within its bounding box.
	/**
	 * If preserving the bounding box is not required, the layout can also be flipped vertically
	 * by calling <tt>scale(1.0, -1.0, false)</tt>.
	 */
	virtual void flipVertical() { flipVertical(boundingBox()); }

	//! Flips the (whole) layout vertically such that the part in \a box remains in this area.
	/**
	 * The whole layout is flipped and then moved such that the part that was in \a box before
	 * flipping is moved to this area.
	 */
	virtual void flipVertical(const DRect &box);

	//! Flips the layout horizontally within its bounding box.
	/**
	* If preserving the bounding box is not required, the layout can also be flipped horizontally
	* by calling <tt>scale(-1.0, 1.0, false)</tt>.
	*/
	virtual void flipHorizontal() { flipHorizontal(boundingBox()); }

	//! Flips the (whole) layout horizontally such that the part in \a box remains in this area.
	/**
	* The whole layout is flipped and then moved such that the part that was in \a box before
	* flipping is moved to this area.
	*/
	virtual void flipHorizontal(const DRect &box);

	//! Scales the layout by (\a sx,\a sy) and then translates it by (\a dx,\a dy).
	/**
	* If \a scaleNodes is true, node sizes are scaled as well. A point (\a x,\a y) is
	* moved to (\a sx &sdot; x + \a dx, \a sy &sdot; y + \a dy).
	*
	* \param sx         is the scaling factor for x-coordinates.
	* \param sy         is the scaling factor for y-coordinates.
	* \param dx         is the translation in x-direction.
	* \param dy         is the translation in y-direction.
	* \param scaleNodes determines if nodes size are scaled as well (true) or not.
	*/
	virtual void scaleAndTranslate(double sx, double sy, double dx, double dy, bool scaleNodes = true);

	//! Scales the layout by \a s and then translates it by (\a dx,\a dy).
	/**
	* If \a scaleNodes is true, node sizes are scaled as well. A point (\a x,\a y) is
	* moved to (\a s &sdot; x + \a dx, \a s &sdot; y + \a dy).
	*
	* \param s          is the scaling factor for both x- and y-coordinates.
	* \param dx         is the translation in x-direction.
	* \param dy         is the translation in y-direction.
	* \param scaleNodes determines if nodes size are scaled as well (true) or not.
	*/
	virtual void scaleAndTranslate(double s, double dx, double dy, bool scaleNodes = true) {
		scaleAndTranslate(s, s, dx, dy, scaleNodes);
	}

	//! Rotates the layout by 90 degree (in clockwise direction) around the origin.
	virtual void rotateRight90();

	//! Rotates the layout by 90 degree (in counter-clockwise direction) around the origin.
	virtual void rotateLeft90();

	//@}
	/**
	 * @name Utility functions
	 */
	//@{

	//! Returns the bounding box of the graph.
	/**
	 * \pre \a nodeGraphics and \a edgeGraphics is enabled
	 */
	virtual DRect boundingBox() const;

	//! Sets the width of all nodes to \a w.
	/**
	 * \pre \a nodeGraphics is enabled
	 */
	void setAllWidth(double w);

	//! Sets the height of all nodes to \a h.
	/**
	 * \pre \a nodeGraphics is enabled
	 */
	void setAllHeight(double h);

	//! Removes all edge bends.
	/**
	 * \pre \a edgeGraphics is enabled
	 */
	void clearAllBends();

	//! Removes unnecessary bend points in orthogonal segements.
	/**
	 * Processes all edges and removes unnecessary bend points in the bend point list
	 * of the edge, i.e., bend points such that the preceding and succeeding bend point
	 * form a horizontal or vertical segement containing this bend point. This function
	 * is useful to remove redundant bend points in an orthogonal layout.
	 *
	 * \pre \a edgeGraphics is enabled
	 *
	 */
	void removeUnnecessaryBendsHV();

	//! Adds additional bend points to all edges for connecting their endpoints.
	/**
	 * According to \a mode switch add either the node center points to
	 * the bends or the anchor point on the node boundary
	 *   - \a mode = 0: only add node center
	 *   - \a mode = 1: compute intersection with the line segment to the center
	 *     and the boundary of the rectangular node
	 *   - \a mode = 2: compute intersection with the first/last line segment
	 *     and the boundary of the rectangular node
	 *
	 * \pre \a nodeGraphics and \a edgeGraphics is enabled
	 */
	void addNodeCenter2Bends(int mode = 1);

	//! Returns true iff \a v represents an association class.
	/**
	 * We hide the internal representation of semantic node types from
	 * the user to be able to change this later (semantic node type member array).
	 * We are not allowed to set association classes manually, only by calling
	 * createAssociationClass().
	 *
	 * \pre \a nodeGraphics are enabled
	 */
	bool isAssociationClass(node v) const {
		return (type(v) == Graph::associationClass);
	}

	//! Returns a list of all inheritance hierarchies in the graph.
	/**
	 * Inheritance hierarchies are identified by edges with type Graph::generalization.
	 *
	 * @param list is a list of all hierarchies; each hierarchie is itself a list
	 *        of all nodes in this hierarchy.
	 *
	 * \return Returns the number of generalization hierarchies.
	 */
	int hierarchyList(List<List<node>*> &list) const;

	//! Returns a list of all inheritance hierarchies in the graph.
	/**
	 * Inheritance hierarchies are identified by edges with type Graph::generalization.
	 *
	 * @param list is a list of all hierarchies; each hierarchie is itself a list
	 *        of all edges in this hierarchy.
	 *
	 * \return Returns the number of generalization hierarchies.
	 */
	int hierarchyList(List<List<edge>*> &list) const;

	//@}
};

} // end namespace ogdf
