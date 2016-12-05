/** \file
 * \brief Declares class LayoutStatistics which provides various
 *        functions for computing statistical measures of a layout.
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

#include <ogdf/basic/GraphAttributes.h>


namespace ogdf {


//! Utility class providing functions for computing statistical information about a layout.
/**
 * @ingroup graph-drawing
 */
class OGDF_EXPORT LayoutStatistics
{
public:
	//! Computes the total edge lengths in the layout given by \a ga.
	/**
	 * \param ga            is the input layout.
	 * \param pMinLength    points to a variable in which the minimum lengths of an edge is stored;
	 *                      \a pMinLength may be a \c nullptr in which case nothing is stored.
	 * \param pMaxLength    points to a variable in which the maximum lengths of an edge is stored;
	 *                      \a pMaxLength may be a \c nullptr in which case nothing is stored.
	 * \param pAvgLength    points to a variable in which the average lengths of an edge is stored;
	 *                      \a pAvgLength may be a \c nullptr in which case nothing is stored.
	 * \param pStdDeviation points to a variable in which the standard deviation of edge lengths is stored;
	 *                      \a pStdDeviation may be a \c nullptr in which case nothing is stored.
	 * \param considerSelfLoops determines whether the lengths of self-loops shall be considered (true) or not (false).
	 * \return the sum of all edge lengths.
	 */
	static double edgeLengths(
		const GraphAttributes &ga,
		double *pMinLength        = nullptr,
		double *pMaxLength        = nullptr,
		double *pAvgLength        = nullptr,
		double *pStdDeviation     = nullptr,
		bool    considerSelfLoops = false);


	//! Computes the total number of bends in the layout given by \a ga.
	/**
	 * \param ga               is the input layout.
	 * \param pMinBendsPerEdge points to a variable in which the minimum lengths of an edge is stored;
	 *                         \a pMinBendsPerEdge may be a \c nullptr in which case nothing is stored.
	 * \param pMaxBendsPerEdge points to a variable in which the maximum lengths of an edge is stored;
	 *                         \a pMaxBendsPerEdge may be a \c nullptr in which case nothing is stored.
	 * \param pAvgBendsPerEdge points to a variable in which the average lengths of an edge is stored;
	 *                         \a pAvgBendsPerEdge may be a \c nullptr in which case nothing is stored.
	 * \param pStdDeviation    points to a variable in which the standard deviation of edge lengths is stored;
	 *                         \a pStdDeviation may be a \c nullptr in which case nothing is stored.
	 * \param considerSelfLoops determines whether the bends of self-loops shall be considered (true) or not (false).
	 * \return the total number of bend points on all edges.
	 */
	static int numberOfBends(
		const GraphAttributes &ga,
		int    *pMinBendsPerEdge  = nullptr,
		int    *pMaxBendsPerEdge  = nullptr,
		double *pAvgBendsPerEdge  = nullptr,
		double *pStdDeviation     = nullptr,
		bool    considerSelfLoops = false);


	//! Computes the angular resolution of the layout given by \a ga.
	/**
	 * The angular resolution of a layout is the smallest angle formed by two edge segments that meet in a common endpoint or bend point.
	 * Angles are given in radians.
	 *
	 * \param ga               is the input layout.
	 * \param pMaxAngle points to a variable in which the largest angle in the layout is stored;
	 *                         \a pMaxAngle may be a \c nullptr in which case nothing is stored.
	 * \param pAvgAngle points to a variable in which the average angle in the layout is stored;
	 *                         \a pAvgAngle may be a \c nullptr in which case nothing is stored.
	 * \param pStdDeviation    points to a variable in which the standard deviation of angles is stored;
	 *                         \a pStdDeviation may be a \c nullptr in which case nothing is stored.
	 * \param considerBends    determines whether bend points of edges shall be considered (true) or not (false).
	 * \return the angular resolution (smallest angle) of the layout.
	 */
	static double angularResolution(
		const GraphAttributes &ga,
		double *pMaxAngle     = nullptr,
		double *pAvgAngle     = nullptr,
		double *pStdDeviation = nullptr,
		bool    considerBends = true);


	//! Computes the number of edge crossings in the layout given by \a ga.
	/**
	 * If several edge segments cross in the same point, this is counted as if all of these segments
	 * would cross pairwise. E.g., if three edge segments cross in a common points, this counts as
	 * three crossings.
	 *
	 * \param ga is the input layout; if it contains bend points, edge segement of an edge's polyline is considered as line segment.
	 * \return the number of line crossings.
	 */
	static int numberOfCrossings(const GraphAttributes &ga);


	//! Computes the intersection graph \a H of the line segments in the layout given by \a ga.
	/**
	 * The nodes of the intersection graph are all the endpoints of segments in the layout and intersection points.
	 * The edges corrsepond to edges in the input layout: If an edge connecting points \a p and \a q in \a H corresponds
	 * to an edge \a e in the input graph, then \a e contains a line segment \a s such that both \a p and \a q are
	 * endpoints or intersection points of \a s.
	 *
	 * To put it more simple, we obtain graph \a H from the input layout \a ga by putting a dummy vertex on each crossing
	 * and bend point, and joining all nodes representing the same point in the plane.
	 *
	 * \param ga        is the input layout; if it contains bend points, edge segement of an edge's polyline is considered as line segment.
	 * \param H         is assigned the intersection graph.
	 * \param points    maps nodes in \a H to their geometric position in the layout.
	 * \param origNode  maps nodes in \a H to nodes in the input graph (given by \a ga); points which are only intersection points of segments are mapped to \c nullptr.
	 * \param origEdge  maps edges in \a H to the corresponding edges in the input graph (given by \a ga).
	 */
	static void intersectionGraph(const GraphAttributes &ga, Graph &H, NodeArray<DPoint> &points, NodeArray<node> &origNode, EdgeArray<edge> &origEdge);
};


} // end namespace ogdf
