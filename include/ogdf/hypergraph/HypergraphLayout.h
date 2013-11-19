/*
 * $Revision: 3830 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 09:55:21 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Layout algorithms for hypergraph based on edge standard
 *        representations (clique / star / tree) - HypergraphLayoutES
 *        and subset standard representation - HypergraphLayoutSS.
 *
 * ... edge standard is based partly on Section 7.2 of PhD Thesis
 * by Dr. Chimani, subset standard is based on the following paper:
 *
 * Bertault, François and Eades, Peter.:Drawing Hypergraphs in the
 * Subset Standard (Short Demo Paper) Graph Drawing Springer Berlin /
 * Heidelberg 2001. pp.45-76. ISBN 978-3-540-41554-1
 *
 * \author Ondrej Moris
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

#ifndef HYPERGRAPH_LAYOUT_H
#define HYPERGRAPH_LAYOUT_H

#include <ogdf/hypergraph/Hypergraph.h>
#include <ogdf/hypergraph/EdgeStandardRep.h>
#include <ogdf/hypergraph/HypergraphAttributes.h>
#include <ogdf/hypergraph/HypergraphLayoutModule.h>

#include <ogdf/basic/exceptions.h>
#include <ogdf/basic/ModuleOption.h>
#include <ogdf/module/EmbedderModule.h>
#include <ogdf/module/CrossingMinimizationModule.h>
#include <ogdf/module/LayoutPlanRepModule.h>

#include <ogdf/planarity/PlanRep.h>

namespace ogdf {

class OGDF_EXPORT HypergraphLayoutSS : public HypergraphLayoutModule {

public:

	enum Method {
		dummyNode    = 0x000001,
		spanningTree = 0x000002,
		steinerTree  = 0x000003
	};

private:

	//!< Determines whether polygons should be convex (if possible).
	bool m_convex;

	//!< Defines the minimum distance between polygons and hypernodes.
	int m_separation;

	//!< Defines the number of algorithm iterations.
	int m_iterations;

	//!< Defines what method is used to represent hyperedges.
	/**
	 * The following representation methods are available:
	 *
	 * 1. Dummy - star based representation of hyperedges such that
	 *            all newly introduced dummy nodes are placed in the
	 *            barycenter of their relevant hypernodes.
	 *
	 * 2. Spanning Tree - each hyperedge is represented by a minumum
	 *            euclidean spanning tree of their hypenodes, no new
	 *            dummy nodes are created.
	 *
	 * 3. Steiner Tree - each hyperedge is represented by a Steiner
	 *            tree such that its leaves are hypernodes incident
	 *            with the hyperedge, steiner vertices are represented
	 *            by newly created nodes.
	 */
	Method representationHelper;

public:

	//! Creates an instance of subset standard hypergraph layout.
	HypergraphLayoutSS();

	//! Copy constructor.
	HypergraphLayoutSS(const HypergraphLayoutSS &hl);

	//! Destructor.
	~HypergraphLayoutSS();

	//! Assignment operator.
	HypergraphLayoutSS &operator=(const HypergraphLayoutSS &hl);

	/**
	 * \brief Calls subset standard hypergraph layout.
	 *
	 * @param HA is the input hypergraph and will also be assigned the
	 *           layout information.
	 */
	void call(HypergraphAttributes &HA)
	{
	  layout(HA);
	}

private:

	/**
	 * Let a hypergraph H be given. The algorithm works as follows:
	 *
	 * 1. Hypernodes as assigned random positions (eg. on a grid).
	 *
	 * 2. Iterate the following (m_iterations):
	 *
	 *    a) Transform H into a simple graph H based on a chosen
	 *       representation method (dummy, spanning tree or steiner trees).
	 *       Make sure hypernodes positions in G are preserved. See below
	 *       for more details about representation methods.
	 *    b) Apply any energy-based layout algorithm to get a layout of G.
	 *    c) Set positions of hypernodes of H according to positions of their
	 *       corresponding nodes of G.
	 *
	 * 3. Again, transform H into a simple graph H based on a chosen
	 *    representation method preserving all hypernodes positions.
	 *
	 * 4. For each hyperedge e, draw a countour around edges of G representing
	 *    e, make sure m_separation is kept between a contour and edges.
	 *
	 * 5. We say that a convex polygon representing a hyperedge e is valid
	 *    when it does contain hypernodes incident with e only. If it m_convex
	 *    is set then transform all contours into convex polygons unless they
	 *    are valid (i.e. compute their convex hulls).
	 *
	 */
	void layout(HypergraphAttributes &HA)
	{
		OGDF_THROW_PARAM(LibraryNotSupportedException, lnscFunctionNotImplemented);
	}
};

class OGDF_EXPORT HypergraphLayoutES : public HypergraphLayoutModule {

public:

	//!< Final appearance is driven by given profile.
	enum Profile {
		Normal          = 0x000001,
		ElectricCircuit = 0x000002
	};

private:

	//!< The ration between width and height of a drawing.
	double m_ratio;

	//!< The number of crossings in the layout.
	int m_crossings;

	//!< Defines whether a drawing IO constraint is desired or not.
	bool m_constraintIO;

	//!< Defines whether inputs and outputs are placed on different "sides".
	// TODO: This might require some tweaks in Hypergraph class.
	bool m_constraintPorts;

	//!< Defines the profile of the layout (eg. Electric Circuit).
	Profile m_profile;

	//!< The module for computing the final layout.
	ModuleOption<LayoutPlanRepModule>  m_planarLayoutModule;

	//!< The module for crossing minimization.
	ModuleOption<CrossingMinimizationModule> m_crossingMinimizationModule;

	//!< The module for embedding planarization.
	ModuleOption<EmbedderModule>  m_embeddingModule;

public:

	// constructor
	HypergraphLayoutES();

	// destructor
	virtual ~HypergraphLayoutES() { }

	// Dynamic casting is currently not working as desired and hence we left
	// the following call inherited from superclass empty.
	virtual void call(HypergraphAttributes &HA);

	//void call(HypergraphAttributesES &HA);

	//! Assignment operator.
	HypergraphLayoutES &operator=(const HypergraphLayoutES &hl);

	//! Returns the number of crossings in computed layout.
	int crossings() const
	{
		return m_crossings;
	}

	//! Returns the ratio  between width and height of a drawing.
	double ratio() const
	{
		return m_ratio;
	}

	//! Sets the Input / Output drawing requirement.
	void setConstraintIO(bool pConstraintIO)
	{
		m_constraintIO = pConstraintIO;
	}

	//! Sets the layout profile.
	void setProfile(Profile pProfile)
	{
		m_profile = pProfile;
	}

	/** @}
	 *  @name Modules
	 *  @{
	 */

	/**
	 * \brief Sets the module option for the planar layout.
	 *
	 * Crossing minimization produces a planar representation of a hypergraph
	 * such that all crossings are replaced by additional dummy nodes.
	 * This is in fact a planar graph and hence it can be drawn quite
	 * easily by any planar layout algorithm.
	 */
	void setPlanarLayoutModule
		(LayoutPlanRepModule *pPlanarLayoutModule)
	{
		m_planarLayoutModule.set(pPlanarLayoutModule);
	}


	/**
	 * \brief Sets the module option for crossing minimization.
	 *
	 * The crossing minimization module minimizes the crossings of a hypergraph
	 * in an edge standard  representation.
	 */
	void setCrossingMinimizationModule
		(CrossingMinimizationModule *pCrossingMinimizationModule)
	{
		m_crossingMinimizationModule.set(pCrossingMinimizationModule);
	}

	/**
	 * \brief Sets the module option for embedding.
	 *
	 * When a planarized edge representation of a hypergraph in computed,
	 * we have to found a way how to embed it (ie. find faces).
	 */
	void setEmbeddingModule
		(EmbedderModule *pEmbeddingModule)
	{
		m_embeddingModule.set(pEmbeddingModule);
	}

private:

	void layout(HypergraphAttributesES &pHA);

	//void planarizeCC(PlanRep &ccPlanarRep, List<edge> &fixedShell);

	void packAllCC(PlanRep &planarRep,
				   HypergraphAttributesES &pHA,
				   Array<DPoint> &bounding);

	std::pair<node,node> * insertShell(GraphCopySimple &planarRep,
									   List<node> &src,
									   List<node> &tgt,
									   List<edge> &fixedShell);

	void removeShell(PlanRep &planarRep, std::pair<node,node> &st);

	void applyProfile(HypergraphAttributesES &HA);

};

} // end namespace ogdf

#endif
