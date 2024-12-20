/** \file
 * \brief Declaration of interface for grid layout algorithms.
 *
 * \author Carsten Gutwenger
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

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/LayoutModule.h>
#include <ogdf/basic/LayoutStandards.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/geometry.h>

namespace ogdf {
class GraphAttributes;
class GridLayout;
class PlanRep;

/**
 * \brief Base class for grid layout algorithms.
 *
 * A grid layout algorithm computes a grid layout of a graph.
 * Such a grid layout does not take real node sizes into account
 * and places a node simply on a grid point; edges may be routed
 * via bend points on grid points.
 *
 * The class GridLayoutModule provides the infrastructure
 * to transform a grid layout into a (usual) layout of a graph,
 * turning a grid layout algorithm automatically into a
 * LayoutModule.
 */
class OGDF_EXPORT GridLayoutModule : public LayoutModule {
	friend class GridLayoutPlanRepModule;
	friend class PlanarGridLayoutModule;

public:
	//! Initializes a grid layout module.
	GridLayoutModule() : LayoutModule(), m_separation(LayoutStandards::defaultNodeSeparation()) { }

	virtual ~GridLayoutModule() { }

	/**
	 * \brief Calls the grid layout algorithm (general call).
	 *
	 * This method implements the call function of the base class LayoutModule.
	 * A derived algorithm implements the call by implementing doCall().
	 *
	 * @param GA is the input graph; the new layout is also stored in \p GA.
	 */
	virtual void call(GraphAttributes& GA) override final;

	/**
	 * \brief Calls the grid layout algorithm (call for GridLayout).
	 *
	 * A derived algorithm implements the call by implementing doCall().
	 *
	 * @param G is the input graph.
	 * @param gridLayout is assigned the computed grid layout.
	 */
	void callGrid(const Graph& G, GridLayout& gridLayout);

	//! Returns the current setting of the minimum distance between nodes.
	/**
	 * This minimum distance is used for mapping grid coordinates to double coordinates as stored
	 * in GraphAttributes. This mapping occurs automatically when the grid layout algorithm is
	 * called with LayoutModule's call method.
	 */
	double separation() const { return m_separation; }

	//! Sets the minimum distance between nodes.
	/**
	 * This minimum distance is used for mapping grid coordinates to double coordinates as stored
	 * in GraphAttributes. This mapping occurs automatically when the grid layout algorithm is
	 * called with LayoutModule's call method.
	 */
	void separation(double sep) { m_separation = sep; }

	const IPoint& gridBoundingBox() const { return m_gridBoundingBox; }

protected:
	/**
	 * \brief Implements the algorithm call.
	 *
	 * A derived algorithm must implement this method and return the computed grid
	 * layout in \p gridLayout.
	 *
	 * @param G is the input graph.
	 * @param gridLayout is assigned the computed grid layout.
	 * @param boundingBox returns the bounding box of the grid layout. The lower left
	 *        corner of the bounding box is always (0,0), thus this IPoint defines the
	 *        upper right corner as well as the width and height of the grid layout.
	 */
	virtual void doCall(const Graph& G, GridLayout& gridLayout, IPoint& boundingBox) = 0;

	IPoint m_gridBoundingBox; //!< The computed bounding box of the grid layout.

private:
	double m_separation; //!< The minimum distance between nodes.

	//! Internal transformation of grid coordinates to real coordinates.
	void mapGridLayout(const Graph& G, GridLayout& gridLayout, GraphAttributes& AG);
};

/**
 * \brief Base class for planar grid layout algorithms.
 *
 * A planar grid layout algorithm is a grid layout algorithm
 * that produces a crossing-free grid layout of a planar
 * graph. It provides an additional call method for producing
 * a planar layout with a predefined planar embedding.
 */
class OGDF_EXPORT PlanarGridLayoutModule : public GridLayoutModule {
public:
	//! Initializes a planar grid layout module.
	PlanarGridLayoutModule() : GridLayoutModule() { }

	virtual ~PlanarGridLayoutModule() { }

	/**
	 * \brief Calls the grid layout algorithm with a fixed planar embedding (general call).
	 *
	 * A derived algorithm implements the call by implementing doCall().
	 *
	 * @param AG is the input graph; the new layout is also stored in \p AG.
	 * @param adjExternal specifies an adjacency entry on the external face,
	 *        or is set to 0 if no particular external face shall be specified.
	 */
	void callFixEmbed(GraphAttributes& AG, adjEntry adjExternal = nullptr);

	/**
	 * \brief Calls the grid layout algorithm with a fixed planar embedding (call for GridLayout).
	 *
	 * A derived algorithm implements the call by implementing doCall().
	 *
	 * @param G is the input graph.
	 * @param gridLayout is assigned the computed grid layout.
	 * @param adjExternal specifies an adjacency entry (of \p G) on the external face,
	 *        or is set to 0 if no particular external face shall be specified.
	 */
	void callGridFixEmbed(const Graph& G, GridLayout& gridLayout, adjEntry adjExternal = nullptr);

protected:
	/**
	 * \brief Implements the algorithm call.
	 *
	 * A derived algorithm must implement this method and return the computed grid
	 * layout in \p gridLayout.
	 *
	 * @param G is the input graph.
	 * @param adjExternal is an adjacency entry on the external face, or 0 if no external
	 *        face is specified.
	 * @param gridLayout is assigned the computed grid layout.
	 * @param boundingBox returns the bounding box of the grid layout. The lower left
	 *        corner of the bounding box is always (0,0), thus this IPoint defines the
	 *        upper right corner as well as the width and height of the grid layout.
	 * @param fixEmbedding determines if the input graph is embedded and that embedding
	 *        has to be preserved (true), or if an embedding needs to be computed (false).
	 */

	virtual void doCall(const Graph& G, adjEntry adjExternal, GridLayout& gridLayout,
			IPoint& boundingBox, bool fixEmbedding) = 0;

	//! Implements the GridLayoutModule::doCall().
	virtual void doCall(const Graph& G, GridLayout& gridLayout, IPoint& boundingBox) override {
		doCall(G, nullptr, gridLayout, boundingBox, false);
	}

	/**
	 * Handles the special cases of graphs with less than 3 nodes.
	 *
	 * @return true iff the instance is handled
	 */
	bool handleTrivial(const Graph& G, GridLayout& gridLayout, IPoint& boundingBox);
};

/**
 * \brief Base class for grid layout algorithms operating on a PlanRep.
 *
 * A GridLayoutPlanRepModule is a special class of a grid layout module
 * that produces a planar layout of a planar graph, and that has a
 * special call method (taking a PlanRep as input) for using the
 * layout module within the planarization approach.
 *
 * \see PlanarizationGridLayout
 */
class OGDF_EXPORT GridLayoutPlanRepModule : public PlanarGridLayoutModule {
public:
	//! Initializes a plan-rep grid layout module.
	GridLayoutPlanRepModule() : PlanarGridLayoutModule() { }

	virtual ~GridLayoutPlanRepModule() { }

	/**
	 * \brief Calls the grid layout algorithm (call for GridLayout).
	 *
	 * The implementation of this call method temporarily constructs a
	 * PlanRep and copies the resulting grid layout to the grid layout for
	 * the input graph.
	 *
	 * @param G is the input graph.
	 * @param gridLayout is assigned the computed grid layout.
	 */
	void callGrid(const Graph& G, GridLayout& gridLayout) {
		PlanarGridLayoutModule::callGrid(G, gridLayout);
	}

	/**
	 * \brief Calls the grid layout algorithm (call for GridLayout of a PlanRep).
	 *
	 * @param PG is the input graph which may be modified by the algorithm.
	 * @param gridLayout is assigned the computed grid layout of \p PG.
	 */
	void callGrid(PlanRep& PG, GridLayout& gridLayout);

	/**
	 * \brief Calls the grid layout algorithm with a fixed planar embedding (call for GridLayout).
	 *
	 * A derived algorithm implements the call by implementing doCall().
	 *
	 * @param G is the input graph.
	 * @param gridLayout is assigned the computed grid layout.
	 * @param adjExternal specifies an adjacency entry (of \p G) on the external face,
	 *        or is set to 0 if no particular external face shall be specified.
	 */
	void callGridFixEmbed(const Graph& G, GridLayout& gridLayout, adjEntry adjExternal = nullptr) {
		PlanarGridLayoutModule::callGridFixEmbed(G, gridLayout, adjExternal);
	}

	/**
	 * \brief Calls the grid layout algorithm with a fixed planar embedding (call for GridLayout of a PlanRep).
	 *
	 * A derived algorithm implements the call by implementing doCall().
	 *
	 * @param PG is the input graph which may be modified by the algorithm.
	 * @param gridLayout is assigned the computed grid layout.
	 * @param adjExternal specifies an adjacency entry (of \p PG) on the external face,
	 *        or is set to 0 if no particular external face shall be specified.
	 */
	void callGridFixEmbed(PlanRep& PG, GridLayout& gridLayout, adjEntry adjExternal = nullptr);

protected:
	/**
	 * \brief Implements the algorithm call.
	 *
	 * A derived algorithm must implement this method and return the computed grid
	 * layout of \p PG in \p gridLayout.
	 *
	 * @param PG is the input graph which may be modified by the algorithm.
	 * @param adjExternal is an adjacency entry on the external face, or 0 if no external
	 *        face is specified.
	 * @param gridLayout is assigned the computed grid layout.
	 * @param boundingBox is assigned the bounding box of the computed layout.
	 * @param fixEmbedding determines if the input graph is embedded and that embedding
	 *        has to be preserved (true), or if an embedding needs to be computed (false).
	 */
	virtual void doCall(PlanRep& PG, adjEntry adjExternal, GridLayout& gridLayout,
			IPoint& boundingBox, bool fixEmbedding) = 0;


	//! Implements PlanarGridLayoutModule::doCall().
	void doCall(const Graph& G, adjEntry adjExternal, GridLayout& gridLayout, IPoint& boundingBox,
			bool fixEmbedding) override;

	using PlanarGridLayoutModule::doCall;
};

}
