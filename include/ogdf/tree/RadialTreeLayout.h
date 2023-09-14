/** \file
 * \brief Declaration of linear time layout algorithm for free
 *        trees (class RadialTreeLayout).
 *
 * Based on chapter 3.1.1 Radial Drawings of Graph Drawing by
 * Di Battista, Eades, Tamassia, Tollis.
 *
 * \author Carsten Gutwenger, Mirko H. Wagner
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

#include <ogdf/basic/LayoutModule.h>
#include <ogdf/basic/SList.h>

namespace ogdf {

//! The radial tree layout algorithm.
/**
 * <H3>Optional parameters</H3>
 * Radial tree layout provides the following optional parameters.
 *
 * <table>
 *   <tr>
 *     <th><i>Option</i><th><i>Type</i><th><i>Default</i><th><i>Description</i>
 *   </tr><tr>
 *     <td><i>levelDistance</i><td>double<td>30.0
 *     <td>The minimal vertical distance required between levels.
 *   </tr><tr>
 *     <td><i>rootSelection</i><td> #RootSelectionType <td> RootSelectionType::Center
 *     <td>Specifies how to select the root of the tree.
 *   </tr>
 * </table>
*/
class OGDF_EXPORT RadialTreeLayout : public LayoutModule {
public:
	//! Selection strategies for root of the tree.
	enum class RootSelectionType {
		Source, //!< Select a source in the graph.
		Sink, //!< Select a sink in the graph.
		Center //!< Select the center of the tree.
	};

private:
	double m_levelDistance; //!< The minimal distance between levels.

	RootSelectionType m_selectRoot; //!< Specifies how to determine the root.

	node m_root; //!< The root of the tree.

	int m_numLevels; //!< The number of levels (root is on level 0).
	NodeArray<int> m_level; //!< The level of a node.
	NodeArray<node> m_parent; //!< The parent of a node (nullptr if root).
	Array<SListPure<node>> m_nodes; //!< The nodes at a level.
	NodeArray<SListPure<node>> m_children; //!< The children of a node.

	NodeArray<double> m_relWidth; //!< The relative width of the subtree.
	//!< A nodes relative width is the greater of the sum of its childrens widths
	//!< and its diameter divided by its childrens level.

	NodeArray<double> m_absWidth; //!< the absolute width of the subtree.

	NodeArray<double> m_angle; //!< The angle of node center (for placement).
	NodeArray<double> m_wedge; //!< The wedge reserved for subtree.

	NodeArray<double> m_diameter; //!< The diameter of a circle bounding a node.
	Array<double> m_maxDiameter; //!< The maximal diameter on a level.

	Array<double> m_radius; //!< The width of a level.

public:
	//! Creates an instance of radial tree layout and sets options to default values.
	RadialTreeLayout();

	//! Copy constructor.
	RadialTreeLayout(const RadialTreeLayout& tl);

	//! Destructor.
	~RadialTreeLayout() = default;

	//! Assignment operator.
	RadialTreeLayout& operator=(const RadialTreeLayout& tl);

	//! Calls the algorithm for graph attributes \p GA.
	/**
	 * The algorithm preserve the order of children which is given by
	 * the adjacency lists.
	 *
	 * \pre The graph is a tree.
	 * @param GA represents the input graph and is assigned the computed layout.
	 */
	virtual void call(GraphAttributes& GA) override;

	// option that determines the minimal vertical distance
	// required between levels

	//! Returns the option <i>levelDistance</i>.
	double levelDistance() const { return m_levelDistance; }

	//! Sets the option <i>levelDistance</i> to \p x.
	void levelDistance(double x) { m_levelDistance = x; }

	// option that determines if the root is on the top or on the bottom

	//! Returns the option <i>rootSelection</i>.
	RootSelectionType rootSelection() const { return m_selectRoot; }

	//! Sets the option <i>rootSelection</i> to \p sel.
	void rootSelection(RootSelectionType sel) { m_selectRoot = sel; }

	const NodeArray<double>& diameter() const { return m_diameter; }

private:
	void FindRoot(const Graph& G);
	void ComputeLevels(const Graph& G);
	void ComputeDiameters(GraphAttributes& AG);
	void ComputeAngles(const Graph& G);
	void ComputeCoordinates(GraphAttributes& AG);
	void ComputeGroupings(const Graph& G);

	OGDF_NEW_DELETE
};

}
