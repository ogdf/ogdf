/*
 * $Revision: 3367 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-04 16:29:19 +0200 (Do, 04. Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration of interface for edge insertion algorithms
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

#ifndef OGDF_EDGE_INSERTION_MODULE_H
#define OGDF_EDGE_INSERTION_MODULE_H


#include <ogdf/basic/Logger.h>
#include <ogdf/basic/Module.h>
#include <ogdf/basic/Timeouter.h>
#include <ogdf/planarity/PlanRepLight.h>


namespace ogdf {

	//! Interface for edge insertion algorithms.
	/**
	 * \see SubgraphPlanarizer
	 */
	class OGDF_EXPORT EdgeInsertionModule : public Module, public Timeouter {
	public:
		//! Initializes an edge insertion module (default constructor).
		EdgeInsertionModule() { }

		//! Initializes an edge insertion module (copy constructor).
		EdgeInsertionModule(const EdgeInsertionModule &eim) : Timeouter(eim) { }

		//! Destructor.
		virtual ~EdgeInsertionModule() { }

		//! Returns a new instance of the edge insertion module with the same option settings.
		virtual EdgeInsertionModule *clone() const = 0;

		//! Inserts all edges in \a origEdges into \a pr.
		/**
		 * @param pr        is the input planarized representation and will also receive the result.
		 * @param origEdges is the array of original edges (edges in the original graph of \a pr)
		 *                  that have to be inserted.
		 * @return the status of the result.
		 */
		ReturnType call(PlanRepLight &pr, const Array<edge> &origEdges) {
			return doCall(pr, origEdges, 0, 0, 0);
		}

		//! Inserts all edges in \a origEdges with given costs into \a pr.
		/**
		 * @param pr        is the input planarized representation and will also receive the result.
		 * @param costOrig  is an edge array containing the costs of original edges; edges in
		 *                  \a pr without an original edge have zero costs.
		 * @param origEdges is the array of original edges (edges in the original graph of \a pr)
		 *                  that have to be inserted.
		 * @return the status of the result.
		 */
		ReturnType call(PlanRepLight &pr,
			const EdgeArray<int> &costOrig,
			const Array<edge>    &origEdges)
		{
			return doCall(pr, origEdges, &costOrig, 0, 0);
		}


		//! Inserts all edges in \a origEdges with given costs and subgraphs (for simultaneous drawing) into \a pr.
		/**
		 * @param pr            is the input planarized representation and will also receive the result.
		 * @param costOrig      is an edge array containing the costs of original edges; edges in
		 *                      \a pr without an original edge have zero costs.
		 * @param origEdges     is the array of original edges (edges in the original graph of \a pr)
		 *                      that have to be inserted.
		 * @param edgeSubGraphs is an edge array specifying to which subgraph an edge belongs.
		 * @return the status of the result.
		 */
		ReturnType call(PlanRepLight &pr,
			const EdgeArray<int>      &costOrig,
			const Array<edge>         &origEdges,
			const EdgeArray<__uint32> &edgeSubGraphs)
		{
			return doCall(pr, origEdges, &costOrig, 0, &edgeSubGraphs);
		}


		//! Inserts all edges in \a origEdges with given forbidden edges into \a pr.
		/**
		 * \pre No forbidden edge may be in \a origEdges.
		 *
		 * @param pr            is the input planarized representation and will also receive the result.
		 * @param forbiddenOrig is an edge array indicating if an original edge is forbidden to be crossed.
		 * @param origEdges     is the array of original edges (edges in the original graph of \a pr)
		 *                      that have to be inserted.
		 * @return the status of the result.
		 */
		ReturnType call(PlanRepLight &pr,
			const EdgeArray<bool> &forbiddenOrig,
			const Array<edge>     &origEdges)
		{
			return doCall(pr, origEdges, 0, &forbiddenOrig, 0);
		}

		//! Inserts all edges in \a origEdges with given costs and forbidden edges into \a pr.
		/**
		 * \pre No forbidden edge may be in \a origEdges.
		 *
		 * @param pr            is the input planarized representation and will also receive the result.
		 * @param costOrig      is an edge array containing the costs of original edges; edges in
		 *                      \a pr without an original edge have zero costs.
		 * @param forbiddenOrig is an edge array indicating if an original edge is forbidden to be crossed.
		 * @param origEdges     is the array of original edges (edges in the original graph of \a pr)
		 *                      that have to be inserted.
		 * @return the status of the result.
		 */
		ReturnType call(PlanRepLight &pr,
			const EdgeArray<int>  &costOrig,
			const EdgeArray<bool> &forbiddenOrig,
			const Array<edge>     &origEdges)
		{
			return doCall(pr, origEdges, &costOrig, &forbiddenOrig, 0);
		}


		//! Inserts all edges in \a origEdges with given costs, forbidden edges, and subgraphs (for simultaneous drawing) into \a pr.
		/**
		 * \pre No forbidden edge may be in \a origEdges.
		 *
		 * @param pr            is the input planarized representation and will also receive the result.
		 * @param costOrig      is an edge array containing the costs of original edges; edges in
		 *                      \a pr without an original edge have zero costs.
		 * @param forbiddenOrig is an edge array indicating if an original edge is forbidden to be crossed.
		 * @param origEdges     is the array of original edges (edges in the original graph of \a pr)
		 *                      that have to be inserted.
		 * @param edgeSubGraphs is an edge array specifying to which subgraph an edge belongs.
		 * @return the status of the result.
		 */
		ReturnType call(PlanRepLight &pr,
			const EdgeArray<int>      &costOrig,
			const EdgeArray<bool>     &forbiddenOrig,
			const Array<edge>         &origEdges,
			const EdgeArray<__uint32> &edgeSubGraphs)
		{
			return doCall(pr, origEdges, &costOrig, &forbiddenOrig, &edgeSubGraphs);
		}


		//! Inserts all edges in \a origEdges into \a pr, optionally costs, forbidden edges, and subgraphs (for simultaneous drawing) may be given.
		/**
		 * @param pr             is the input planarized representation and will also receive the result.
		 * @param origEdges      is the array of original edges (edges in the original graph of \a pr)
		 *                       that have to be inserted.
		 * @param pCostOrig      points to an edge array containing the costs of original edges; edges in
		 *                       \a pr without an original edge have zero costs. May be a 0-pointer, in which case all edges have cost 1.
		 * @param pForbiddenOrig points to an edge array indicating whether an original edge is forbidden to be crossed.
		 *                       May be a 0-pointer, in which case no edges are forbidden.
		 * @param pEdgeSubGraphs points to an edge array specifying to which subgraph an edge belongs.
		 *                       May be a 0-poiner, in which case no subgraphs / simultaneous embedding is used.
		 * @return the status of the result.
		 */
		ReturnType callEx(
			PlanRepLight              &pr,
			const Array<edge>         &origEdges,
			const EdgeArray<int>      *pCostOrig = 0,
			const EdgeArray<bool>     *pForbiddenOrig = 0,
			const EdgeArray<__uint32> *pEdgeSubGraphs = 0)
		{
			return doCall(pr, origEdges, pCostOrig, pForbiddenOrig, pEdgeSubGraphs);
		}


	protected:
		//! Actual algorithm call that has to be implemented by derived classes.
		/**
		 * @param pr             is the input planarized representation and will also receive the result.
		 * @param origEdges      is the array of original edges (edges in the original graph of \a pr)
		 *                       that have to be inserted.
		 * @param pCostOrig      points to an edge array containing the costs of original edges; edges in
		 *                       \a pr without an original edge have zero costs.
		 * @param pForbiddenOrig points to an edge array indicating whether an original edge is forbidden to be crossed.
		 * @param pEdgeSubGraphs points to an edge array specifying to which subgraph an edge belongs.
		 * @return the status of the result.
		 */
		virtual ReturnType doCall(
			PlanRepLight              &pr,
			const Array<edge>         &origEdges,
			const EdgeArray<int>      *pCostOrig,
			const EdgeArray<bool>     *pForbiddenOrig,
			const EdgeArray<__uint32> *pEdgeSubGraphs) = 0;


		OGDF_MALLOC_NEW_DELETE
	};

} // end namespace ogdf

#endif
