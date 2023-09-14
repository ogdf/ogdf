/** \file
 * \brief Generator for visualizing graphs using LaTeX/TikZ.
 *
 * \author Hendrik Brueckler
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.md in the OGDF root directory for details.
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
 * License along with this program; if not, see
 * http://www.gnu.org/copyleft/gpl.html
 */

#pragma once

#include <ogdf/fileformats/GraphIO.h>

namespace ogdf {

/**
 * \brief LaTeX+TikZ Writer
 *
 * Generates and outputs standalone LaTeX documents/graphics using Tikz.
 */
class TikzWriter {
public:
	/**
	 * @brief Units of length available in LaTeX
	 */
	enum class LengthUnit {
		PT, //!< 0.3515 mm
		MM, //!< a millimeter
		CM, //!< a centimeter
		IN, //!< an inch
		EX, //!< roughly the height of an 'x' (lowercase) in the current font (it depends on the font used)
		EM, //!< roughly the width of an 'M' (uppercase) in the current font (it depends on the font used)
		MU //!< math unit equal to 1/18 em, where em is taken from the math symbols family
	};

	/**
	 * Construct a new writer for drawing graphs using graph attributes
	 *
	 * @param attr the attributes of the graph to draw
	 * @param unit the unit to use for LaTeX lengths
	 */
	TikzWriter(const GraphAttributes& attr, LengthUnit unit = LengthUnit::MM)
		: m_attr(attr), m_clsAttr(nullptr), m_unit(unit), m_nodeStyles(), m_edgeStyles() { }

	/**
	 * Construct a new writer for drawing cluster graphs using graph attributes
	 *
	 * @param attr the attributes of the cluster graph to draw
	 * @param unit the unit to use for LaTeX lengths
	 */
	TikzWriter(const ClusterGraphAttributes& attr, LengthUnit unit = LengthUnit::MM)
		: m_attr(attr), m_clsAttr(&attr), m_unit(unit), m_nodeStyles(), m_edgeStyles() { }

	/**
	 * Output the member graph to an output stream in LaTeX+TikZ format
	 *
	 * @param os output stream
	 * @return true iff successfull
	 */
	bool draw(std::ostream& os);

private:
	/**
	 * Wrap LaTeX/TikZ header/footer around a tikzpicture and write to output stream.
	 *
	 * @param os output stream
	 * @param tikzPic the content of the tikzPicture
	 * @param uniformStyle whether not to set style per node
	 * @param uniformWidth whether not to set width per node
	 * @param uniformHeight whether not to set height per node
	 */
	void wrapHeaderFooter(std::ostream& os, std::string tikzPic, bool uniformStyle,
			bool uniformWidth, bool uniformHeight) const;

	/**
	 * Draws a rectangle for each cluster in the ogdf::ClusterGraph.
	 *
	 * @param os output stream
	 */
	void drawAllClusters(std::ostream& os);

	/**
	 * Draws each node of the graph.
	 *
	 * @param os output stream
	 * @param uniformStyle whether not to set style per node
	 * @param uniformWidth whether not to set width per node
	 * @param uniformHeight whether not to set height per node
	 */
	void drawAllNodes(std::ostream& os, bool uniformStyle, bool uniformWidth, bool uniformHeight);

	/**
	 * Draws a sequence of lines for each edge in the graph.
	 *
	 * @param os output stream
	 */
	void drawAllEdges(std::ostream& os);

	/**
	 * Draws a cluster as a rectangle, using its size and style properties.
	 *
	 * @param os output stream
	 * @param c cluster to draw
	 */
	void drawCluster(std::ostream& os, cluster c);

	/**
	 * Draws a node using its shape, size and style properties.
	 *
	 * @param os output stream
	 * @param v node to draw
	 * @param uniformStyle whether not to set style per node
	 * @param uniformWidth whether not to set width per node
	 * @param uniformHeight whether not to set height per node
	 */
	void drawNode(std::ostream& os, node v, bool uniformStyle, bool uniformWidth, bool uniformHeight);

	/**
	 * Draws a sequence of lines for each edge in the graph.
	 *
	 * @param os output stream
	 * @param e edge to draw
	 */
	void drawEdge(std::ostream& os, edge e);

	/**
	 * Get the style of a cluster in TikZ syntax.
	 *
	 * @param c cluster
	 * @return std::string style of \p c in TikZ syntax
	 */
	std::string getClusterStyle(cluster c) const;

	/**
	 * Get the shape of a node in TikZ syntax.
	 *
	 * @param v node
	 * @return std::string shape of \p v in TikZ syntax
	 */
	std::string getNodeShape(node v) const;

	/**
	 * Get the style of a node in TikZ syntax.
	 *
	 * @param v node
	 * @return std::string style of \p v in TikZ syntax
	 */
	std::string getNodeStyle(node v) const;

	/**
	 * Get the label of a node in TikZ syntax.
	 *
	 * @param v node
	 * @return std::string label of \p v in TikZ syntax
	 */
	std::string getNodeLabel(node v) const;

	/**
	 * Get the total width the node text may occupy.
	 *
	 * @param v node
	 * @return double width that may be occupied by label text for node \p v
	*/
	double getTextWidth(node v) const;

	/**
	 * Get the style of an edge in TikZ syntax.
	 *
	 * @param e edge
	 * @return std::string style of \p e in TikZ syntax
	 */
	std::string getEdgeStyle(edge e) const;

	/**
	 * Get the arrows of an edge in TikZ syntax.
	 *
	 * @param e edge
	 * @return std::string arrows of \p e in TikZ syntax
	 */
	std::string getEdgeArrows(edge e) const;

	/**
	 * Get the label of an edge in TikZ syntax, positioned as a node on the edge path.
	 *
	 * @param e edge
	 * @param previousPoint previous point on edge path
	 * @param labelPoint where to put the edge label
	 * @return std::string label of \p e in TikZ syntax
	 */
	std::string getEdgeLabel(edge e, const DPoint& previousPoint, const DPoint& labelPoint) const;

	/**
	 * Check whether a point (e.g. edge bend point) lies within a node (using
	 * node shapes with same size and aspect as in TikZ).
	 *
	 * @param p point to check
	 * @param v node to check
	 * @return true iff \p lies within the border of \p v
	 */
	bool isCoveredBy(const DPoint& p, node v) const;

	/**
	 * @brief Calculates the arrow size to be used for TikZ arrows
	 *
	 * @return double
	 */
	double calcArrowSize() const;

	/**
	 * @brief Mainly avoid scientific notation (not handled by LaTeX) and add length unit mm.
	 *
	 * @param f floating point number
	 * @return std::string \p f in mm as formatted string
	 */
	std::string texLength(double f) const;

	/**
	 * @brief Convert an ogdf::StrokeType, strokeWidth and ogdf::Color to a a line style string in TikZ syntax
	 *
	 * @param strokeType stroke type
	 * @param strokeWidth stroke width
	 * @param strokeColor stroke color
	 * @return std::string line style in TikZ syntax
	 */
	std::string getLineStyle(StrokeType strokeType, double strokeWidth, Color strokeColor) const;

	/**
	 * Convert an ogdf::Color to a string in TikZ syntax
	 *
	 * @param c color
	 * @return std::string \p c in TikZ syntax
	 */
	static std::string getColorString(Color c);

	//! attributes of the graph to be visualized
	const GraphAttributes& m_attr;

	//! attributes of the cluster graph (\c nullptr if no cluster graph given)
	const ClusterGraphAttributes* m_clsAttr;

	//! The LaTeX unit to use for all ocurring lengths
	LengthUnit m_unit;

	//! to avoid as much redundancy as possible, any occurring node style will be predefined and reused
	std::vector<std::string> m_nodeStyles;

	//! to avoid as much redundancy as possible, any occurring edge style will be predefined and reused
	std::vector<std::string> m_edgeStyles;
};

}
