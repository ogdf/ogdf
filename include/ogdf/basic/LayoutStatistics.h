/** \file
 * \brief Declares class LayoutStatistics which provides various
 *        functions for computing statistical measures of a layout.
 *
 * \author Carsten Gutwenger, Sharif Wurz
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

#include <ogdf/basic/ArrayBuffer.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/geometry.h>

namespace ogdf {
class GraphAttributes;

//! Computes statistical information about a layout.
/**
 * @ingroup graph-drawing
 */
class OGDF_EXPORT LayoutStatistics {
public:
	//! Computes graph height, in coordinate length metric.
	/**
	 * 
	 * The graph height is defined as the difference between the
	 * maximum and minimum y-coordinates of the nodes in the layout.
	 * 
	 * Returns graph height double value.
	 * Returns 0 if graph has one or less nodes.
	 * 
	 * \param G                 Input graph.
	 * \param ga                Input layout.
	 * \return                  The graph height.
	 */
	static inline double graphHeight(const Graph& G, const GraphAttributes& ga);

	//! Computes graph width, in coordinate length metric.
	/**
	 * 
	 * The graph width is defined as the difference between the
	 * maximum and minimum x-coordinates of the nodes in the layout.
	 * 
	 * Returns graph width double value.
	 * Returns 0 if graph has one or less nodes.
	 * 
	 * \param G                 Input graph.
	 * \param ga                Input layout.
	 * \return                  The graph width.
	 */
	static inline double graphWidth(const Graph& G, const GraphAttributes& ga);

	//! Computes graph area, in coordinate length metric.
	/**
	 * 
	 * Uses graphHeight() and graphWidth() as helper functions,
	 * to compute area of graph.
	 * In effect graphArea() = graphHeight() * graphWidth().
	 * 
	 * Returns graph area double value
	 * Returns 0 if graph has one or less nodes.
	 * 
	 * \param G                 Input graph.
	 * \param ga                Input layout.
	 * \return                  The graph width.
	 */
	static inline double graphArea(const Graph& G, const GraphAttributes& ga);

	//! Computes the edge length for each edge in the layout \p ga.
	/**
	 * \param ga                Input layout.
	 * \param considerSelfLoops Determines whether the lengths of self-loops are considered.
	 * \return                  The edge length for each edge.
	 */
	static ArrayBuffer<double> edgeLengths(const GraphAttributes& ga, bool considerSelfLoops = false);


	//! Computes the number of bends (i.e. bend-points) for each edge in the layout \p ga.
	/**
	 * \param ga                Input layout.
	 * \param considerSelfLoops Determines whether the bends of self-loops are considered.
	 * \return                  The number of bends for each edge.
	 */
	static ArrayBuffer<int> numberOfBends(const GraphAttributes& ga, bool considerSelfLoops = false);


	//! Computes the angle for each pair of adjacent edge segments of the layout \p ga.
	/**
	 * Angles are given in radians.
	 *
	 * \param ga            Input layout.
	 * \param considerBends Determines whether bend points of edges shall be considered.
	 * \return              The angle for each two adjacent edge segments.
	 */
	static ArrayBuffer<double> angles(const GraphAttributes& ga, bool considerBends = true);


	//! Computes the number of edge crossings for each edge in the layout \p ga.
	/**
	 * If several edge segments cross in the same point, this is counted as if
	 * all of these segments would cross pairwise. E.g., if three edge segments
	 * cross in a common points, this counts as two crossings for each of the
	 * edges.
	 *
	 * \warning The same warning as for #intersectionGraph applies.
	 * \warning The sum of all returned values is twice the number of crossings
	 * as each crossing involves two edges.
	 *
	 * \param ga Input layout. If it contains bend points, each segment of an edge's polyline is considered as a line segment.
	 *           Otherwise, a straight-line drawing is assumed.
	 * \return   The number of crossings for each edge.
	 */
	static ArrayBuffer<int> numberOfCrossings(const GraphAttributes& ga);


	//! Computes the number of crossings through a non-incident node for each
	//! edge in the layout \p ga.
	/**
	 * If several edge segments cross a node in the same point, one crossing per
	 * edge segment is counted. E.g., if three edge segments cross a node in a
	 * common point, this counts as three node crossings.
	 * Each node is treated as if it had the shape of the rectangle with the
	 * corresponding width and height given by \p ga.
	 *
	 * \param ga Input layout. If it contains bend points, each segment of an edge's polyline is considered as a line segment.
	 *           Otherwise, a straight-line drawing is assumed.
	 * \return   The number of node crossings for each edge.
	 */
	static ArrayBuffer<int> numberOfNodeCrossings(const GraphAttributes& ga);


	//! Computes the percentage of crossings for each edge in the layout \p ga
	//! compared to the maximum number of crossings.
	/**
	 * The maximum number of crossings is defined as the number of crossings
	 * that would occur if all edges were straight lines and crossed at most once.
	 * This is computed by the numberOfCrossings() function.
	 *
	 * \warning sum of all returned values is twice the number of crossings
	 * as each crossing involves two edges (the edge itself and the crossing edge),
	 * therefore we divide the sum by two to get the correct amount.
	 *
	 *
	 * \return The percentage, of crossings for each edge compared to the maximum amount of crossings.
	 *
	 * Source:
	 * https://drops.dagstuhl.de/storage/00lipics/lipics-vol320-gd2024/LIPIcs.GD.2024.45/LIPIcs.GD.2024.45.pdf
	 * 
	 */
	static double LayoutStatistics::percentageCrossingVsMaxCrossings(const GraphAttributes& ga);

	//! Computes the number of node overlaps for each node in the layout \p ga.
	/**
	 * Each node is treated as if it had the shape of the rectangle with the
	 * corresponding width and height given by \p ga.
	 *
	 * \warning The sum of all returned values is twice the number of node
	 * overlaps as each node overlap involves two nodes.
	 *
	 * \param ga Input layout.
	 * \return   The number of node overlaps for each node.
	 */
	static ArrayBuffer<int> numberOfNodeOverlaps(const GraphAttributes& ga);


	//! Computes the intersection graph \p H of the line segments in the layout given by \p ga.
	/**
	 * The nodes of the intersection graph are all endpoints of segments in \p ga plus all intersection points.
	 * The edges correspond to edges in the input layout: If an edge connecting points \a v and \a w in \p H corresponds
	 * to an edge \a e in the input graph, then \a e contains a line segment \a s such that both \a v and \a w are
	 * endpoints or intersection points of \a s.
	 *
	 * To put it more simple, we obtain graph \p H from \p ga by putting a dummy vertex on each crossing
	 * and bend point, and joining all nodes representing the same point in the plane.
	 *
	 * \warning Do not call this algorithm on drawings with arbitrarily close curves (e.g., curves overlapping on an interval).
	 *
	 * \param ga        Input layout. If it contains bend points, each segment of an edge's polyline is considered as a line segment.
	 *                  Otherwise, a straight-line drawing is assumed.
	 * \param H         Is assigned the intersection graph.
	 * \param points    Maps nodes in \p H to their geometric position in the layout.
	 * \param origNode  Maps nodes in \p H to nodes in \p ga's graph.
	 * Points that are only intersection points of segments are mapped to \c nullptr.
	 * \param origEdge  Maps edges in \p H to the corresponding edges in \p ga's graph.
	 */
	static void intersectionGraph(const GraphAttributes& ga, Graph& H, NodeArray<DPoint>& points,
			NodeArray<node>& origNode, EdgeArray<edge>& origEdge);


	//! Computes the edge length deviation \p H of the edges in the graph given, in \p ga.
	/**
 	 * Edge length deviation def.: The deviation of the actual edge length from the average edge length of the graph.
	 * Returns an array of edge deviations (of doubles) of all edges of the graph.
	 * Also generates a private variable with an edge length deviation average.
	 * 
	 * Source: https://www2.cs.arizona.edu/people/kobourov/gd-metrics2024.pdf
	 * 
	 * Returns an array of deviation values, of each edge (deviation of avg edge length) in the Graph
	*/
	static ArrayBuffer<double> edgeLengthDeviation(const GraphAttributes& ga);


	//! Computes the distances between each pair of nodes in the graph given, in \p ga.
	/**
	 * Returns an array of distances, a pair where the first entry is the node pair, and the second entry is the distance between them as doubles.
	 * Nodes are added in the order they are in the graph container "mainGraph.nodes".
	 * Per default, each edge is added uniquely -> ( \a u, \a v ) or ( \a v , \a u ),
	 * else (bidirectional = true) both edges ( \a u, \a v ) and ( \a v , \a u ) are added.
	 * handles size allocation, so no need to preallocate specific vector size.
	 * 
	*/
	static void LayoutStatistics::distancesBetweenAllNodes(const Graph& mainGraph,
			const GraphAttributes& ga,
			ArrayBuffer<std::pair<std::pair<node, node>, double>>& allDistances,
			bool edgesTwice = false);

	//! Computes neighborhood preservation \p H of the edges in the graph given, in \p ga.
	/**
	 * Neighborhood preservation def.: How many nodes in the graph are connected to the actually closest (least coordinative distance) nodes.
	 * Returns an array of neighbourhood preservations in percentile (of doubles from \p 0.0 for lowest preservation, to \p 1.0 for highest preservation) of all nodes of the graph.
	 * Neighborhood preservation def.: How many nodes in the graph are connected to the actually closest nodes (least coordinative distance.
	 * Also generates a private variable with an overall graph neighbourhood preservation average (also from \p 0.0 to \p 1.0 ).
	 * 
	 * Source: https://www2.cs.arizona.edu/people/kobourov/gd-metrics2024.pdf
	 * 
	 * Returns an array of preservation values, of each node in Graph \p mainGraph
	 * 
	*/
	static ArrayBuffer<double> neighbourhoodPreservation(const Graph& mainGraph,
			const GraphAttributes& ga);

	//! Computes Gabriel Ratio \p H of the edges in the graph given, in \p ga and also returns new set of edges that fulfill Gabriel criteria.
	/**
	 * Gabriel Ratio def.: 	The ratio of the distance between two nodes and the distance to the closest node to the edge between the two nodes,
	 * 						that is only allowed to contain the two nodes itself, not any other ones.
	 * 
	 * Source: https://www2.cs.arizona.edu/people/kobourov/gd-metrics2024.pdf
	 * 
	 * Returns an array of Gabriel Ratios (of doubles) of all nodes of the graph
	 * Also generates a private variable with an overall graph Gabriel Ratio average
	 */
	static ArrayBuffer<double> gabrielRatio(Graph& mainGraph, const GraphAttributes& ga);

	//! Computes Node Ratio \p H of the nodes in the graph given, in \p ga.
	/**
	 * Node Ratio def.: 	Node Resolution (NR) metric is the ratio between
	 * 						the smallest distance between two nodes and the
	 * 						largest distance between two nodes.
	 * Source: https://www2.cs.arizona.edu/people/kobourov/gd-metrics2024.pdf
	 * 
	 * Returns a double that gives the Node Ratio as output, which is: smallest_distance_between_2_nodes / biggest_distance_between_2_nodes
	 * 
	 */
	static double LayoutStatistics::nodeResolution(Graph& mainGraph, const GraphAttributes& ga);


	//! Computes Angular Resolution (AR) \p H of the nodes in the graph given, in \p ga.
	/**
	 * def. Angular Resolution (AR):
	 * "After calculating the angular deviation of each edge from the ideal angle
	 * (based on degree), the AR metric calculates the mean over all nodes
	 * (excluding degree-1 nodes)." - https://www2.cs.arizona.edu/people/kobourov/gd-metrics2024.pdf
	 * 
	 * Returns a double of the angular resolution of the graph
	 * Returns 0.0 if there are no nodes with degree greater than 2, or if the graph has less than 3 nodes.
	 * 
	 */
	static double LayoutStatistics::angularResolution(const Graph& mainGraph,
			const GraphAttributes& ga);


	//! Computes Aspect Ratio (Asp) \p H of the layout/graph \p g.
	/**
	 * def. Aspect Ratio (Asp):
	 * "The Aspect Ratio is the ratio of the height of
	 * the drawing’s bounding box to its width (or vice versa,
	 * depending on which is greater)." - https://www2.cs.arizona.edu/people/kobourov/gd-metrics2024.pdf
	 * 
	 * Returns aspect ration double 
	 * Returns 0.0 if the graph has less than 2 nodes, or either height or width of the bounding box are equal or smaller than 0.
	 * 
	 */
	static double LayoutStatistics::aspectRatio(const Graph& mainGraph, const GraphAttributes& ga);

	//! Computes Node Uniformity (NU) \p H of the graph \p g.
	/**
	 * def. Node Uniformity (NU):
	 * Calculates how well the node distribution is in the bounding box of the graph.
	 * "The Node Uniformity metric measures the distribution of nodes in the bounding
	 * box, by splitting it into grid cells based on the number of nodes
	 * in the graph, counting the number of nodes in each cell, and
	 * comparing them with an ideal distribution." - https://www2.cs.arizona.edu/people/kobourov/gd-metrics2024.pdf
	 * 
	 * Can also receive grid width and height as params, defaults to 10x10 grid.
	 * 
	 * Returns node uniformity metric between 1.0 and 0.0, where 1.0 is a perfect uniformity
	 * and 0.0 is the worst uniformity.
	 * Returns 0.0 if the graph has less than 2 nodes, or number of grid cells is 0 (e.g. gridWidth or/and gridHeight is 0).
	 * 
	 */
	static double LayoutStatistics::nodeUniformity(const Graph& mainGraph,
			const GraphAttributes& ga, size_t gridWidth = 10, size_t gridHeight = 10);

	//! Computes Edge Orthogonality (EO) \p H of the graph \p g.
	/**
	 * def. Edge Orthogonality (EO):
	 * Calculates mean of edge orthogonality deviation.
	 * "The Edge Orthogonality metric takes the mean of the angular deviation
	 * of all edges from the horizontal or vertical axis (whichever is
	 * closest)." - https://www2.cs.arizona.edu/people/kobourov/gd-metrics2024.pdf
	 * 
	 * Returns mean Edge Orthogonality of all edges
	 */
	static double LayoutStatistics::edgeOrthogonality(const Graph& mainGraph,
			const GraphAttributes& ga);

	//! Computes center of mass, where most nodes are.
	/**
	 * def. Center Of Mass:
	 * Calculates mean position of all nodes.
	 * 
	 * Source:
	 * https://drops.dagstuhl.de/storage/00lipics/lipics-vol320-gd2024/LIPIcs.GD.2024.45/LIPIcs.GD.2024.45.pdf
	 * 
	 * Returns double pair containing center of mass coordinates.
	 * Returns a pair with (0.0, 0.0) if graph is empty.
	 */
	static std::pair<double, double> LayoutStatistics::centerOfMass(const Graph& mainGraph,
			const GraphAttributes& ga);


	//! Computes Closest pair of points.
	/**
	 * def. Closest pair of points:
	 * Finds the two nodes in the graph, that have the smallest euclidean distance to each other,
	 * and calculates their distance.
	 * 
	 * Source:
	 * https://drops.dagstuhl.de/storage/00lipics/lipics-vol320-gd2024/LIPIcs.GD.2024.45/LIPIcs.GD.2024.45.pdf
	 *
	 * Returns euclidean distance double of two closest nodes.
	 * Returns -1.0 if graph is empty.
	 */
	static double LayoutStatistics::closestPairOfPoints(const Graph& mainGraph,
			const GraphAttributes& ga);

	//! Computes horizontal node balance.
	/**
	 * def. Horizontal/Vertical node balance:
	 * Split nodes into two groups, left and right (or above or below) of the center,
	 * and divide it by total number of nodes.
	 * 
	 * Source:
	 * https://drops.dagstuhl.de/storage/00lipics/lipics-vol320-gd2024/LIPIcs.GD.2024.45/LIPIcs.GD.2024.45.pdf
	 *
	 * 
	 * Returns balance double value.
	 * Returns -1.0 if graph is empty.
	 */
	static double LayoutStatistics::horizontalVerticalBalance(const Graph& mainGraph,
			const GraphAttributes& ga, const bool vertical = false);

	//! Retrieves min and max, x- and y-coordinates.
	/**
	 * 
	 * def. Border coordinates:
	 * Finds the min and max x- and y-coordinates of all nodes in the graph,
	 * and returns them as pairs. They make up the bounding box of the graph.
	 * 
	 * Source:
	 * https://drops.dagstuhl.de/storage/00lipics/lipics-vol320-gd2024/LIPIcs.GD.2024.45/LIPIcs.GD.2024.45.pdf
	 *
	 * 
	 * Returns a pair of pairs, where the first pair contains minX & maxX coordinates,
	 * and the other pair contains minY & maxY coordinates.
	 * Returns a pair of pairs where each double is -1.0 when graph is empty.
	 */
	static std::pair<std::pair<double, double>, std::pair<double, double>>
	LayoutStatistics::borderCoordinates(const Graph& mainGraph, const GraphAttributes& ga);

	//! Calculating percentage of nodes with integer coordinates.
	/**
	 * 
	 * For each node, checks if its x- and y-coordinates are integers,
	 * and counts how many of them are.
	 * Can also receives epsilon value for tolerance of closeness to integer (default is 1e-9).
	 * 
	 * Source:
	 * https://drops.dagstuhl.de/storage/00lipics/lipics-vol320-gd2024/LIPIcs.GD.2024.45/LIPIcs.GD.2024.45.pdf
	 *
	 * Returns percentage of nodes with integer coordinates (or within epsilon range).
	 * Returns -1.0 when graph is empty.
	 */
	static double LayoutStatistics::nodeOrthogonality(const Graph& mainGraph,
			const GraphAttributes& ga, const double epsilon = 1e-9);

	//! Calculates mean edge direction (vector) angle.
	/**
	 * 
	 * For each edge, calculates angle of its directional vector,
	 * and returns the mean angle of all edges.
	 * 
	 * Source:
	 * https://drops.dagstuhl.de/storage/00lipics/lipics-vol320-gd2024/LIPIcs.GD.2024.45/LIPIcs.GD.2024.45.pdf
	 *
	 * Returns degree of mean edge direction angle (as double).
	 * Returns -1.0 when graph is empty.
	 */
	double LayoutStatistics::averageFlow(const Graph& mainGraph, const GraphAttributes& ga);

	//! Calculates percentage of edges that point upwards.
	/**
	 * 
	 * For each edge, checks if its directional vector has a positive y-value,
	 * meaning that its direction points upwards.
	 * 
	 * Can also receives epsilon value for tolerance of closeness to integer (default is 1e-9).
	 * 
	 * Source:
	 * https://drops.dagstuhl.de/storage/00lipics/lipics-vol320-gd2024/LIPIcs.GD.2024.45/LIPIcs.GD.2024.45.pdf
	 *
	 * Returns percentage of edges pointing upwards.
	 * Returns -1.0 when graph is empty, or undirected.
	 */
	double LayoutStatistics::upwardsFlow(const Graph& mainGraph, const GraphAttributes& ga);


	//! Calculates the variance of node distances from the center of mass.
	/**
	 * 
	 * For each node, compute distance to center of mass,
	 * then calculate variance in all these distances.
	 * 
	 * Source:
	 * https://drops.dagstuhl.de/storage/00lipics/lipics-vol320-gd2024/LIPIcs.GD.2024.45/LIPIcs.GD.2024.45.pdf
	 *
	 * Returns concentration nodes in graph.
	 * Returns -1.0 when graph is empty.
	 */
	double LayoutStatistics::concentration(const Graph& mainGraph, const GraphAttributes& ga);
};
}
