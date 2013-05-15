/*
 * $Revision: 3433 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-22 13:44:53 +0200 (Mo, 22. Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration of class PlanarizationLayout.
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


#ifdef _MSC_VER
#pragma once
#endif

#ifndef OGDF_PLANARIZATION_LAYOUT_LIGHT_H
#define OGDF_PLANARIZATION_LAYOUT_LIGHT_H


#include <ogdf/module/LayoutModule.h>
#include <ogdf/module/CrossingMinimizationModule.h>
#include <ogdf/module/EmbedderModule.h>
#include <ogdf/module/LayoutPlanRepModule.h>
#include <ogdf/module/CCLayoutPackModule.h>
#include <ogdf/basic/ModuleOption.h>



namespace ogdf {

	class CliqueReplacer;


//! The planarization approach for drawing graphs.
class OGDF_EXPORT PlanarizationLayout : public LayoutModule
{
public:
	//! Creates an instance of planarization layout and sets options to default values.
	PlanarizationLayout();

	//! Destructor.
	~PlanarizationLayout() { }

	//! Calls planarization layout for GraphAttributes \a ga.
	/**
	 * \pre The graph has no self-loops.
	 * @param ga is the input graph and will also be assigned the layout information.
	 */
	void call(GraphAttributes &ga);

	void call(GraphAttributes &ga, GraphConstraints & gc) { call(ga); }

	//! Calls planarization layout with clique handling for GraphAttributes \a ga with associated graph \a g.
	/**
	 * \pre \a g is the graph associated with graph attributes \a ga.
	 *
	 * This call perfoms a special handling for cliques, which are temporarily replaced by a star graph.
	 * In the final drawing, the clique edges are drawn straight-line.
	 */
	void call(GraphAttributes &ga, Graph &g);

	void callSimDraw(GraphAttributes &ga);

	/** @}
	 *  @name Optional parameters
	 *  @{
	 */

	//! Returns the current setting of option pageRatio.
	/**
	 * This option specifies the desired ration width / height of the computed
	 * layout. It is currently only used for packing connected components.
	 */
	double pageRatio() const {
		return m_pageRatio;
	}

	//! Sets the option pageRatio to \a ratio.
	void pageRatio(double ratio) {
		m_pageRatio = ratio;
	}

	//! Returns the current setting of option minCliqueSize.
	/**
	 * If preprocessing of cliques is considered, this option determines the
	 * minimal size of cliques to search for.
	 */
	int minCliqueSize() const {
		return m_cliqueSize;
	}

	//! Set the option minCliqueSize to \a i.
	void minCliqueSize(int i) {
		m_cliqueSize = max(i, 3);
	}


	/** @}
	 *  @name Module options
	 *  @{
	 */

	//! Sets the module option for crossing minimization.
	void setCrossMin(CrossingMinimizationModule *pCrossMin) {
		m_crossMin.set(pCrossMin);
	}

	//! Sets the module option for the graph embedding algorithm.
	/**
	 * The result of the crossing minimization step is a planar graph,
	 * in which crossings are replaced by dummy nodes. The embedding
	 * module then computes a planar embedding of this planar graph.
	 */
	void setEmbedder(EmbedderModule *pEmbedder) {
		m_embedder.set(pEmbedder);
	}

	//! Sets the module option for the planar layout algorithm.
	/**
	 * The planar layout algorithm is used to compute a planar layout
	 * of the planarized representation resulting from the crossing
	 * minimization step. Planarized representation means that edge crossings
	 * are replaced by dummy nodes of degree four, so the actual layout
	 * algorithm obtains a planar graph as input. By default, the planar
	 * layout algorithm produces an orthogonal drawing.
	 */
	void setPlanarLayouter(LayoutPlanRepModule *pPlanarLayouter) {
		m_planarLayouter.set(pPlanarLayouter);
	}

	//! Sets the module option for the arrangement of connected components.
	/**
	 * The planarization layout algorithm draws each connected component of
	 * the input graph seperately, and then arranges the resulting drawings
	 * using a packing algorithm.
	 */
	void setPacker(CCLayoutPackModule *pPacker) {
		m_packer.set(pPacker);
	}

	/** @}
	 *  @name Further information
	 *  @{
	 */

	//! Returns the number of crossings in the computed layout.
	int numberOfCrossings() const {
		return m_nCrossings;
	}

	//! @}

private:
	void arrangeCCs(PlanRep &PG, GraphAttributes &GA, Array<DPoint> &boundingBox) const;
	void preprocessCliques(Graph &G, CliqueReplacer &cliqueReplacer);
	void fillAdjNodes(List<node>& adjNodes,
		PlanRep& PG,
		node centerNode,
		NodeArray<bool>& isClique,
		Layout& drawing);

	//! The module for computing a planar subgraph.
	ModuleOption<CrossingMinimizationModule> m_crossMin;

	//! The module for planar embedding.
	ModuleOption<EmbedderModule> m_embedder;

	//! The module for computing a planar layout.
	ModuleOption<LayoutPlanRepModule> m_planarLayouter;

	//! The module for arranging connected components.
	ModuleOption<CCLayoutPackModule> m_packer;

	double m_pageRatio;    //!< The desired page ratio.
	int m_nCrossings;      //!< The number of crossings in the computed layout.

	int m_cliqueSize;      //!< The minimum size of cliques to search for.
};

} // end namespace ogdf


#endif
