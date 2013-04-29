/*
 * $Revision: 3472 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-29 15:52:12 +0200 (Mo, 29. Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration of interface for planar subgraph algorithms.
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

#ifndef OGDF_PLANAR_SUBGRAPH_MODULE_H
#define OGDF_PLANAR_SUBGRAPH_MODULE_H



#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/Module.h>
#include <ogdf/basic/Logger.h>
#include <ogdf/basic/Timeouter.h>


namespace ogdf {

/**
 * \brief Interface for planar subgraph algorithms.
 *
 * \see PlanarizationLayout, PlanarizationGridLayout
 */
class OGDF_EXPORT PlanarSubgraphModule : public Module, public Timeouter {

	int m_maxThreads;	//!< The maximal number of used threads.

public:
	//! Initializes a planar subgraph module (default constructor).
	PlanarSubgraphModule() {
#ifdef OGDF_MEMORY_POOL_NTS
		m_maxThreads = 1;
#else
		m_maxThreads = System::numberOfProcessors();
#endif
	}

	//! Initializes a planar subgraph module (copy constructor).
	PlanarSubgraphModule(const PlanarSubgraphModule &psm) : Timeouter(psm) {
		m_maxThreads = psm.m_maxThreads;
	}

	// destruction
	virtual ~PlanarSubgraphModule() { }


	/**
	 * \brief Returns the set of edges \a delEdges which have to be deleted to obtain the planar subgraph.
	 * @param G is the input graph.
	 * @param cost are the costs of edges.
	 * @param preferedEdges are edges that should be contained in the planar subgraph.
	 * @param delEdges is the set of edges that need to be deleted to obtain the planar subgraph.
	 * @param preferedImplyPlanar indicates that the edges \a preferedEdges induce a planar graph.
	 */
	ReturnType call(const Graph &G,
		const EdgeArray<int> &cost,
		const List<edge> &preferedEdges,
		List<edge> &delEdges,
		bool preferedImplyPlanar = false)
	{
		return doCall(G,preferedEdges,delEdges,&cost,preferedImplyPlanar);
	}


	/**
	 * \brief Returns the set of edges \a delEdges which have to be deleted to obtain the planar subgraph.
	 * @param G is the input graph.
	 * @param preferedEdges are edges that should be contained in the planar subgraph.
	 * @param delEdges is the set of edges that need to be deleted to obtain the planar subgraph.
	 * @param preferedImplyPlanar indicates that the edges \a preferedEdges induce a planar graph.
	 */
	ReturnType call(const Graph &G,
		const List<edge> &preferedEdges,
		List<edge> &delEdges,
		bool preferedImplyPlanar = false)
	{
		return doCall(G,preferedEdges,delEdges,0,preferedImplyPlanar);
	}


	/**
	 * \brief Returns the set of edges \a delEdges which have to be deleted to obtain the planar subgraph.
	 * @param G is the input graph.
	 * @param cost are the costs of edges.
	 * @param delEdges is the set of edges that need to be deleted to obtain the planar subgraph.
	 */
	ReturnType call(const Graph &G, const EdgeArray<int> &cost, List<edge> &delEdges) {
		List<edge> preferedEdges;
		return doCall(G,preferedEdges,delEdges, &cost);
	}

	/**
	 * \brief Returns the set of edges \a delEdges which have to be deleted to obtain the planar subgraph.
	 * @param G is the input graph.
	 * @param delEdges is the set of edges that need to be deleted to obtain the planar subgraph.
	 */
	ReturnType call(const Graph &G, List<edge> &delEdges) {
		List<edge> preferedEdges;
		return doCall(G,preferedEdges,delEdges);
	}


	//! Returns the set of edges \a delEdges which have to be deleted to obtain the planar subgraph.
	ReturnType operator()(const Graph &G,
		const List<edge> &preferedEdges,
		List<edge> &delEdges,
		bool preferedImplyPlanar = false)
	{
		return call(G,preferedEdges,delEdges,preferedImplyPlanar);
	}

	//! Returns the set of edges \a delEdges which have to be deleted to obtain the planar subgraph.
	ReturnType operator()(const Graph &G, List<edge> &delEdges) {
		return call(G,delEdges);
	}


	/**
	 * \brief Makes \a G planar by deleting edges.
	 * @param GC is a copy of the input graph.
	 * @param preferedEdges are edges in \a GC that should be contained in the planar subgraph.
	 * @param delOrigEdges is the set of original edges whose copy has been deleted in \a GC.
	 * @param preferedImplyPlanar indicates that the edges \a preferedEdges induce a planar graph.
	 */
	ReturnType callAndDelete(GraphCopy &GC,
		const List<edge> &preferedEdges,
		List<edge> &delOrigEdges,
		bool preferedImplyPlanar = false);

	/**
	 * \brief Makes \a G planar by deleting edges.
	 * @param GC is a copy of the input graph.
	 * @param delOrigEdges is the set of original edges whose copy has been deleted in \a GC.
	 */
	ReturnType callAndDelete(GraphCopy &GC, List<edge> &delOrigEdges) {
		List<edge> preferedEdges;
		return callAndDelete(GC,preferedEdges,delOrigEdges);
	}

	//! Returns a new instance of the planar subgraph module with the same option settings.
	virtual PlanarSubgraphModule *clone() const = 0;

	//! Returns the maximal number of used threads.
	int maxThreads() const { return m_maxThreads; }

	//! Sets the maximal number of used threads to \a n.
	void maxThreads(int n) {
#ifndef OGDF_MEMORY_POOL_NTS
		m_maxThreads = n;
#endif
	}

protected:
	// computes set of edges delEdges, which have to be deleted
	// in order to get a planar subgraph; edges in preferedEdges
	// should be contained in planar subgraph
	// must be implemented by derived classes!
	/**
	 * \brief Computes the set of edges \a delEdges which have to be deleted to obtain the planar subgraph.
	 *
	 * This is the actual algorithm call and needs to be implemented
	 * by derived classes.
	 * @param G is the input graph.
	 * @param preferedEdges are edges that should be contained in the planar subgraph.
	 * @param delEdges is the set of edges that need to be deleted to obtain the planar subgraph.
	 * @param pCost is apointer to an edge array containing the edge costs; this pointer
	 *        can be 0 if no costs are given (all edges have cost 1).
	 * @param preferedImplyPlanar indicates that the edges \a preferedEdges induce a planar graph.
	 */
	virtual ReturnType doCall(const Graph &G,
		const List<edge> &preferedEdges,
		List<edge> &delEdges,
		const EdgeArray<int>  *pCost = 0,
		bool preferedImplyPlanar = false) = 0;



	OGDF_MALLOC_NEW_DELETE
};

} // end namespace ogdf

#endif
