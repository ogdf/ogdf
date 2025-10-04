/** \file
 * \brief Implements class LayoutStatistics which provides various
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

#include <ogdf/basic/ArrayBuffer.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/LayoutStatistics.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/Math.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/geometry.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include <math.h>

namespace ogdf {

// Forward declaration for borderCoordinates() function
std::pair<std::pair<double, double>, std::pair<double, double>> borderCoordinates(
		const GraphAttributes& ga);

inline double graphHeight(const GraphAttributes& ga) {
	const Graph& G = ga.constGraph();
	if (G.numberOfNodes() <= 1) {
		return 0.0; // height is 0, for one or less nodes
	}
	DRect boundingBoxVals = ga.GraphAttributes::boundingBox();
	return boundingBoxVals.height();
}

inline double graphWidth(const GraphAttributes& ga) {
	const Graph& G = ga.constGraph();
	if (G.numberOfNodes() <= 1) {
		return 0.0; // width is 0, for one or less nodes
	}
	DRect boundingBoxVals = ga.GraphAttributes::boundingBox();
	return boundingBoxVals.width();
}

inline double graphArea(const GraphAttributes& ga) { return graphHeight(ga) * graphWidth(ga); }

ArrayBuffer<double> LayoutStatistics::edgeLengths(const GraphAttributes& ga, bool considerSelfLoops) {
	ArrayBuffer<double> values;
	for (const edge& e : ga.constGraph().edges) {
		if (!considerSelfLoops && e->isSelfLoop()) {
			continue;
		}

		const DPolyline& dpl = ga.bends(e);
		DPoint pv = ga.point(e->source());
		DPoint pw = ga.point(e->target());

		double len = 0;
		if (!dpl.empty()) {
			len = dpl.length();
			len += pv.distance(dpl.front());
			len += pw.distance(dpl.back());
		} else {
			len = pv.distance(pw);
		}

		values.push(len);
	}

	return values;
}

ArrayBuffer<int> LayoutStatistics::numberOfBends(const GraphAttributes& ga, bool considerSelfLoops) {
	ArrayBuffer<int> values;
	for (const edge& e : ga.constGraph().edges) {
		if (considerSelfLoops || !e->isSelfLoop()) {
			values.push(ga.bends(e).size());
		}
	}

	return values;
}

ArrayBuffer<double> LayoutStatistics::angles(const GraphAttributes& ga, bool considerBends) {
	ArrayBuffer<double> values;
	const Graph& G = ga.constGraph();

	for (const node& v : G.nodes) {
		double vx = ga.x(v);
		double vy = ga.y(v);

		// Get angles for edge segments incident to v.
		List<double> angles;
		for (const adjEntry& adj : v->adjEntries) {
			const DPolyline& dpl = ga.bends(adj->theEdge());
			double ex, ey;
			if (dpl.empty()) {
				ex = ga.x(adj->twinNode());
				ey = ga.y(adj->twinNode());
			} else {
				ex = dpl.front().m_x;
				ey = dpl.front().m_y;
			}

			angles.pushBack(atan2(ex - vx, ey - vy));
		}

		if (angles.size() < 2) {
			continue;
		}

		angles.quicksort();
		double lastAngle = angles.back();
		for (const double& psi : angles) {
			double alpha = psi - lastAngle;

			// happens in the first iteration only
			if (alpha < 0) {
				OGDF_ASSERT(psi == angles.front());
				alpha += 2 * Math::pi;
			}

			values.push(alpha);
			lastAngle = psi;
		}
	}

	if (considerBends) {
		for (const edge& e : G.edges) {
			DPolyline dpl = ga.bends(e);

			dpl.pushFront(ga.point(e->source()));
			dpl.pushBack(ga.point(e->target()));
			dpl.normalize();

			if (dpl.size() < 3) {
				continue;
			}

			for (ListConstIterator<DPoint> it = dpl.begin().succ(); it.succ().valid(); ++it) {
				double bx = (*it).m_x, by = (*it).m_y;

				const DPoint& p1 = *it.pred();
				double psi1 = atan2(p1.m_x - bx, p1.m_y - by);

				const DPoint& p2 = *it.succ();
				double psi2 = atan2(p2.m_x - bx, p2.m_y - by);

				double alpha = fabs(psi1 - psi2);
				if (alpha > Math::pi) {
					alpha -= Math::pi;
				}

				values.push(alpha);
				values.push(alpha * Math::pi);
			}
		}
	}

	return values;
}

ArrayBuffer<int> LayoutStatistics::numberOfCrossings(const GraphAttributes& ga) {
	ArrayBuffer<int> values;
	const Graph& G = ga.constGraph();
	EdgeArray<int> crossings(G, 0);

	Graph H;
	NodeArray<DPoint> points;
	NodeArray<node> origNode;
	EdgeArray<edge> origEdge;
	intersectionGraph(ga, H, points, origNode, origEdge);

	for (const node& v : H.nodes) {
		node vOrig = origNode[v];
		int d = (vOrig != nullptr) ? vOrig->degree() : 0;
		int k = (v->degree() - d) / 2;

		// If there are two or more intersecting edges:
		if (k >= 2) {
			// For every original edge involved in the crossing:
			for (const adjEntry& adj : v->adjEntries) {
				if (adj->isSource()) {
					edge e = adj->theEdge();
					edge eOrig = origEdge[e];

					// Ignore original edges incident to vOrig.
					if (eOrig->source() != e->source() || eOrig->target() != e->target()) {
						crossings[eOrig] += (k - 1);
					}
				}
			}
		}
	}

	for (const edge& e : G.edges) {
		values.push(crossings[e]);
	}

	return values;
}

ArrayBuffer<int> LayoutStatistics::numberOfNodeCrossings(const GraphAttributes& ga) {
	ArrayBuffer<int> values;
	const Graph& G = ga.constGraph();
	DPoint inter;

	// Get bounding rectangle of every node.
	NodeArray<DRect> nodeRects(G);
	ga.nodeBoundingBoxes<DRect>(nodeRects);

	// For all edges, get the target point of each of their edge segments.
	for (const edge& e : G.edges) {
		int nCrossingsE = 0;
		node src = e->source();
		node tgt = e->target();
		DPoint vPoint = ga.point(src);

		DPolyline edgeSegmentTargets = ga.bends(e);
		edgeSegmentTargets.pushBack(ga.point(tgt));

		int i = 0;
		int last = edgeSegmentTargets.size() - 1;

		// For all edge segments from vPoint to wPoint:
		for (const DPoint& wPoint : edgeSegmentTargets) {
			DSegment segment = DSegment(vPoint, wPoint);

			// Count crossing of segment with nodes u, but do not count
			// "crossing" of source/target node with first/last edge segment.
			for (const node& u : G.nodes) {
				if ((u != src || i != 0) && (u != tgt || i != last)
						&& nodeRects[u].intersection(segment)) {
					nCrossingsE++;
				}
			}
			vPoint = wPoint;
			i++;
		}
		values.push(nCrossingsE);
	}
	return values;
}

double LayoutStatistics::percentageCrossingVsMaxCrossings(const GraphAttributes& ga) {
	const Graph& G = ga.constGraph();
	size_t m = G.numberOfEdges();
	if (m < 2) { // No crossings possible
		return 0.0; // Avoid division by zero
	}
	// Get crossings for each edge, contains number of crossings for each edge twice
	ArrayBuffer<int> crossings = LayoutStatistics::numberOfCrossings(ga);
	// ArrayBuffer<int> maxCrossings; // Initialize
	size_t sumOfCrossings = 0;

	// Calculate sum of actual crossings for all edges
	for (const int& val : crossings) {
		sumOfCrossings += val;
	}
	sumOfCrossings /= 2; // Each crossing is counted twice, so we divide by 2

	size_t sumMaxCrossings = 0;
	// Calculate maximum crossings for all edges (non-incident edges)
	// Explanation: The number of all edges which are not adjacent to the edge, and thus able to cross
	// each edge, can get crossed by all the m other edges
	for (const edge& e : G.edges) {
		const size_t degreeU = (e->source())->degree(); // degree of node u
		const size_t degreeV = (e->target())->degree(); // degree of node v

		// Number of edges incident to u and v
		size_t incident = (degreeU + degreeV > 2) ? (degreeU + degreeV - 1) : 0;

		// Max number of possible crossings to edge e
		size_t maxCrossing = (m > incident) ? (m - incident) : 0;

		// sum up all max crossings
		sumMaxCrossings += maxCrossing;
	}
	if (sumMaxCrossings == 0) { // if no crossings possible
		return 0.0;
	}
	return (static_cast<double>(sumOfCrossings) / static_cast<double>(sumMaxCrossings))
			* 100.0; // Return percentage
}

ArrayBuffer<int> LayoutStatistics::numberOfNodeOverlaps(const GraphAttributes& ga) {
	ArrayBuffer<int> values;
	const Graph& G = ga.constGraph();

	// Get bounding rectangle of every node.
	NodeArray<DIntersectableRect> nodeRects(G);
	ga.nodeBoundingBoxes<DIntersectableRect>(nodeRects);

	// For all pairs of nodes, test whether they overlap.
	for (const node& v : G.nodes) {
		int nOverlapsV = 0;
		for (const node& w : G.nodes) {
			if (v != w && nodeRects[v].intersects(nodeRects[w])) {
				nOverlapsV++;
			}
		}
		values.push(nOverlapsV);
	}

	return values;
}

ArrayBuffer<double> LayoutStatistics::edgeLengthDeviation(const GraphAttributes& ga) {
	ArrayBuffer<double> edgeLengths = LayoutStatistics::edgeLengths(ga, true);
	double edgeSum = 0.0;
	double delta = 0.0;
	size_t edgesNum = edgeLengths.size();

	for (const double& len : edgeLengths) {
		edgeSum += len;
	}

	// if edgesNum (amount of edges) > 0, then edgeSum / edgeCount, else avgEdgeLen = 0.0
	double avgEdgeLen = (edgesNum > 0) ? (edgeSum / edgesNum) : 0.0;

	ArrayBuffer<double> edgeLenDev; // array of deviation values, of each edge
	double edgeLenDevVal = 0.0; // for one edge deviation value, the lower the better

	// For all edges, subtract all deltas of edge lengths to avg.
	for (const double& len : edgeLengths) {
		// calculating delta for each edge using std::fabs, for floating numbers
		delta = std::fabs(len - avgEdgeLen);
		// adding delta to array of edge length deviations
		edgeLenDev.push(delta);
		edgeLenDevVal += (delta);
	}

	edgeLenDevVal /= edgesNum; // calculating average edge length deviation
	return edgeLenDev;
}

void LayoutStatistics::distancesBetweenAllNodes(const GraphAttributes& ga,
		ArrayBuffer<std::pair<std::pair<node, node>, double>>& allDistances, bool bidirectional) {
	const Graph& mainGraph = ga.constGraph();
	const size_t n = mainGraph.numberOfNodes(); // for efficiency
	if (n < 2) {
		return; // No pairs possible
	}

	// Sizes: bidirectional {(u_i,v_i),(v_i,u_i)} = n * n, unique {(u_i,v_i)} = n * (n - 1) / 2
	const size_t reqSize = (bidirectional ? n * n : (n * (n - 1)) / 2);

	// Correct memory handling, to avoid back and forth memory allocation if array is too small
	if (bidirectional && static_cast<size_t>(allDistances.size()) < reqSize) {
		ArrayBuffer<std::pair<std::pair<node, node>, double>> tmp(reqSize);
		allDistances = std::move(tmp); // transfer reserved memory (size n * n)
		tmp.clear(); // clear temp ArrayBuffer
	} else if (!bidirectional && static_cast<size_t>(allDistances.size()) < reqSize) {
		ArrayBuffer<std::pair<std::pair<node, node>, double>> tmp(reqSize);
		allDistances = std::move(tmp); // transfer reserved memory (size n * (n - 1) / 2)
		tmp.clear(); // clear temp ArrayBuffer
	}

	double dx = 0.0; // x distance between two nodes
	double dy = 0.0; // y distance between two nodes
	double distance = 0.0; // distance between two nodes

	// for each node, calculate distance to all other nodes
	for (auto it_i = mainGraph.nodes.begin(); it_i != mainGraph.nodes.end(); ++it_i) {
		const node& u = *it_i;

		/* start from 0 for visiting each edge twice from each direction if edgesTwice is true
		 else start from i + 1 for only visiting each edge only once */
		for (auto it_j = (bidirectional ? mainGraph.nodes.begin() : std::next(it_i));
				it_j != mainGraph.nodes.end(); ++it_j) {
			if (bidirectional && u == *it_j) {
				continue; // skip same nodes
			}

			const node& v = *it_j; // second node

			// calculating euclidean distance between two nodes
			dx = ga.x(u) - ga.x(v);
			dy = ga.y(u) - ga.y(v);
			// using hypot for safety { hypot(arg1^2, arg2^2) }
			distance = std::hypot(dx, dy); // distance calculation

			allDistances.push(std::make_pair(std::make_pair(u, v), distance));
		}
	}
}

ArrayBuffer<double> LayoutStatistics::neighbourhoodPreservation(const GraphAttributes& ga) {
	const Graph& mainGraph = ga.constGraph();
	double nodePreservationValue = 0.0; // preservation value for each node, higher = better
	double graphPreservationValue = 0.0; // preservation value for whole graph, higher = better
	double preservationValue = 0.0;
	ArrayBuffer<double> nodePreservations; // array of preservation values, of each node/vertex

	size_t currentNodeDeg = 0; // current node degree
	size_t numOfNodes = mainGraph.numberOfNodes(); // for efficiency

	// Init array, where each entry is a pair that contains the edge as pair of nodes, and the distance between the node pair
	ArrayBuffer<std::pair<std::pair<node, node>, double>> allDistances;
	// calculating distances between all nodes, edges are added twice {(u, v), (v, u)}
	distancesBetweenAllNodes(ga, allDistances, true);

	size_t i = 0;
	auto startIt = allDistances.begin();
	auto endIt = allDistances.begin() + numOfNodes;

	// for each node, sorting all by distance
	for (; i < numOfNodes; ++i) {
		// sorting all distances ascendingly for current node within spectrum
		std::sort(startIt, endIt,
				[](const std::pair<std::pair<node, node>, double>& a,
						const std::pair<std::pair<node, node>, double>& b) {
					return a.second < b.second;
				});
		startIt = endIt;
		endIt += numOfNodes;
	}

	node currentNode; // current node
	node neighbor; // neighbor node
	i = 0; // counter for current node index

	// for each node, calculate preservation value
	for (auto it_i = mainGraph.nodes.begin(); it_i != mainGraph.nodes.end(); ++it_i) {
		nodePreservationValue = 0.0; // reseting value for each node/round

		currentNode = *it_i; //current node
		currentNodeDeg = currentNode->degree(); // set current node degree

		if (currentNodeDeg == 0) {
			nodePreservations.push(0.0); // if has no neighbors, skip and add 0.0
			continue;
		}

		size_t start = i * numOfNodes; // current start index within allDistances array

		// for |deg(currentNode)|, calculate preservation value for current node
		for (size_t j = 0; j < currentNodeDeg; ++j) {
			neighbor = allDistances[start + j].first.second; // neighbor of current node

			if (mainGraph.searchEdge(currentNode, neighbor) != nullptr) {
				nodePreservationValue += 1;
				continue;
			}
		}
		preservationValue = nodePreservationValue / currentNodeDeg;
		nodePreservations.push(preservationValue);
		graphPreservationValue += preservationValue;
		++i;
	}

	// final graph preservation value calculation
	graphPreservationValue /= numOfNodes;
	return nodePreservations;
}

ArrayBuffer<double> LayoutStatistics::gabrielRatio(const GraphAttributes& ga,
		Graph& gabrielGraphReference) {
	const Graph& mainGraph = ga.constGraph();
	// array for per-node Gabriel ratios
	ArrayBuffer<double> nodeGabrielRatios;

	size_t numOfNodes = mainGraph.numberOfNodes(); // for efficiency

	ArrayBuffer<double> gabrielCount(numOfNodes, 0.0); // array for counting Gabriel edges for each node
	ArrayBuffer<double> incidentCount(numOfNodes,
			0.0); // array for counting incident edges for each node

	// all nodes distances array, ((u, v), distance)
	ArrayBuffer<std::pair<std::pair<node, node>, double>> allDistances;

	// Calculating distances between all nodes, edges are added uniquely, allDistances size: n*(n-1)/2
	distancesBetweenAllNodes(ga, allDistances);

	// creating new graph for Gabriel edges & adding all nodes
	Graph gabrielGraph;
	for (node n : mainGraph.nodes) {
		gabrielGraph.newNode(n->index());
	}

	bool isGabrielEdge = true; // Gabriel edge true by default

	// for each edge, calculate Gabriel edge
	for (const auto& [edgeNodes, distance] : allDistances) {
		isGabrielEdge = true; // set to true for each loop

		node u = edgeNodes.first; // first node of edge
		node v = edgeNodes.second; // second node of edge

		// incrementing incident edges for both nodes (u, v)
		incidentCount[u->index()] += 1.0;
		incidentCount[v->index()] += 1.0;

		// calculating midpoints (between the two nodes) both of x and y, and the radius range of the distance
		double midPointX = (ga.x(u) + ga.x(v)) / 2.0;
		double midPointY = (ga.y(u) + ga.y(v)) / 2.0;
		double radius = distance / 2.0;

		for (const node& n : mainGraph.nodes) {
			if (n == v || n == u) { // skip same nodes
				continue;
			}
			// calculating distance to midpoint
			double distToMidPoint = std::hypot(ga.x(n) - midPointX, ga.y(n) - midPointY);

			// if distance to midpoint is less than radius => not a Gabriel edge
			if (distToMidPoint < radius) {
				isGabrielEdge = false;
				break;
			}
		}

		if (isGabrielEdge) {
			gabrielGraph.newEdge(u, v); // adding Gabriel edge to output graph
			// incrementing Gabriel count for both nodes (u, v)
			gabrielCount[u->index()] += 1.0;
			gabrielCount[v->index()] += 1.0;
		}
	}

	// calculating Gabriel ratio for each node
	for (size_t i = 0; i < numOfNodes; ++i) {
		double ratio = (incidentCount[i] > 0) ? (gabrielCount[i] / incidentCount[i]) : 0.0;
		nodeGabrielRatios.push(ratio);
	}

	// calculating Gabriel ratio for whole graph
	double graphGabrielRatio = 0.0;
	for (size_t i = 0; i < numOfNodes; ++i) {
		graphGabrielRatio += nodeGabrielRatios[i];
	}
	graphGabrielRatio /= numOfNodes; // final division for Gabriel Ratio of graph

	gabrielGraphReference = gabrielGraph; // setting mainGraph to output graph (gabrielGraph)

	// returning array of Gabriel Ratios for each node
	return nodeGabrielRatios;
}

double LayoutStatistics::nodeResolution(const GraphAttributes& ga) {
	const Graph& mainGraph = ga.constGraph();
	size_t numOfNodes = mainGraph.numberOfNodes(); // for efficiency

	// if graph has less than 3 nodes, no sensible resolution exists as the longest and shortest distance
	// between nodes would be the same
	if (numOfNodes < 3) {
		return 0.0;
	}
	ArrayBuffer<std::pair<std::pair<node, node>, double>> allDistances;
	distancesBetweenAllNodes(ga, allDistances); // each edge once, allDistances size: n*(n-1)/2

	double minDist = std::numeric_limits<double>::max(); // for min. distance between two nodes
	double maxDist = 0.0; // for max distance between two nodes

	for (const auto& pairDist : allDistances) // for each edge/node-pair
	{
		const double& dist = pairDist.second; // distance between nodes
		if (dist < minDist) {
			minDist = dist; // update min. distance
		}
		if (dist > maxDist) {
			maxDist = dist; // update max. distance
		}
	}

	// prevents division by zero, else calculates and returns Node Resolution (NR)
	return minDist > 0.0 ? (minDist / maxDist) : 0.0;
}

double LayoutStatistics::angularResolution(const GraphAttributes& ga) {
	const Graph& mainGraph = ga.constGraph();

	// if graph has less than 3 nodes, there cannot exist two edges
	if (mainGraph.numberOfNodes() < 3) {
		return 0.0;
	}

	double perfectAngle = 0.0; // optimal angle var
	double totalAngleDev = 0.0; // total angle deviation for all nodes
	size_t nodesCount = 0; // count for nodes with degree greater than 1

	// Iterate through each node in graph, and calculate angle deviations for each
	for (const node& u : mainGraph.nodes) {
		const size_t degree = u->degree(); // degree of the source node
		if (degree < 2) {
			continue; // skip nodes with degree less than 2
		}

		List<double> angles; // list of angles for the current node
		// Calculate angles for each adjacent node
		for (const adjEntry& adj : u->adjEntries) {
			node v = adj->twinNode(); // v becomes a neighbor of u
			double x = ga.x(u) - ga.x(v); // x coordinate difference
			double y = ga.y(u) - ga.y(v); // y coordinate difference
			angles.pushBack((atan2(y, x) * 180 / M_PI)); // convert to degrees
		}

		// Sort angles in ascending order
		angles.quicksort();

		angles.pushBack(angles.front() + 2 * Math::pi); // for wrapping around to first angle at the end

		// Ideal angle calculation
		perfectAngle = 2 * Math::pi / degree; // optimal angle in degrees
		double nodeDev = 0.0;

		// Calculate the angle deviation for current node
		for (ListConstIterator<double> it = angles.begin(); it.succ().valid();
				++it) { // checks if next angle exists (last angle is the first angle we added)

			// Calc angle deviation: angle2 - angle1
			double deviationGap = *(it.succ()) - *it;
			nodeDev += std::fabs(deviationGap - perfectAngle); // add absolute deviation to nodeDev
		}

		totalAngleDev += nodeDev / degree; // average angle deviation for current node
		++nodesCount; // increment count for nodes with degree greater than 1
	}

	// Return mean angle deviation ratio (or 0.0 if no nodes with degree greater than 1)
	return (nodesCount > 0) ? (totalAngleDev / nodesCount) : 0.0;
}

double LayoutStatistics::aspectRatio(const GraphAttributes& ga) {
	const Graph& mainGraph = ga.constGraph();
	size_t numOfNodes = mainGraph.numberOfNodes(); // for efficiency

	// if graph has less than 2 nodes, aspect ratio is trivial
	if (numOfNodes < 2) {
		if (numOfNodes == 1) {
			return 1.0; // aspect ratio is 1.0 for single node
		}

		return 0.0; // aspect ratio is 0.0 for no node
	}

	auto maxVals = borderCoordinates(ga);
	const double minX = maxVals.first.first; // min. x coordinate
	const double maxX = maxVals.first.second; // max. x coordinate
	const double minY = maxVals.second.first; // min. y coordinate
	const double maxY = maxVals.second.second; // max. y coordinate

	// width and height of the bounding box
	const double width = maxX - minX; // width of bounding box
	const double height = maxY - minY; // height of bounding box

	// if width or height is <=0,
	// return 0.0 to avoid division by zero and negative values
	if (width <= 0.0 || height <= 0.0) {
		return 0.0;
	}

	return width / height; // aspect ratio
}

double LayoutStatistics::nodeUniformity(const GraphAttributes& ga, size_t gridWidth,
		size_t gridHeight) {
	const Graph& mainGraph = ga.constGraph();
	size_t numOfNodes = mainGraph.numberOfNodes(); // for efficiency

	// if graph has less than 2 nodes, node uniformity is trivial
	if (numOfNodes < 2) {
		return 0.0;
	}

	auto maxVals = borderCoordinates(ga);
	const double minX = maxVals.first.first; // min. x coordinate
	const double maxX = maxVals.first.second; // max. x coordinate
	const double minY = maxVals.second.first; // min. y coordinate
	const double maxY = maxVals.second.second; // max. y coordinate

	// width and height of the bounding box
	// adding 1 for minimum cell size, in case minX == maxX or minY == maxY (e.g. points could be on the same line)
	double width = maxX - minX + 1; // width of bounding box
	double height = maxY - minY + 1; // height of bounding box

	// if width or height is <=0,
	// return 0.0 to avoid division by zero and negative values
	if (width <= 0.0 || height <= 0.0) {
		return 0.0;
	}
	// size of each cell in the grid
	double cellWidth = width / gridWidth; // width of each cell
	double cellHeight = height / gridHeight; // height of each cell

	size_t gridCount = gridWidth * gridHeight; // total number of cells in the grid

	if (gridCount == 0) // if gridCount = 0, return 0.0
	{
		return 0.0;
	}

	std::vector<size_t> nodeCount(gridCount, 0); // array for counting nodes in each cell

	// Iterate through each node, and count nodes in each cell
	for (const node& n : mainGraph.nodes) {
		double ux = ga.x(n);
		double uy = ga.y(n);
		// calculate cell index for each node, and round down
		size_t cellX = static_cast<size_t>((ux - minX) / cellWidth);
		size_t cellY = static_cast<size_t>((uy - minY) / cellHeight);

		// Clamp cell indices to valid range, in case some node is exactly on border
		if (cellX >= gridWidth) {
			cellX = gridWidth - 1;
		}
		if (cellY >= gridHeight) {
			cellY = gridHeight - 1;
		}

		// increment node count for grid-cell
		nodeCount[cellX + cellY * gridWidth]++; // cellY * gridWidth because cells are stored in row-major order
	}

	// calculate ideal uniformity
	size_t idealUniformity = ceil(static_cast<double>(numOfNodes) / static_cast<double>(gridCount));

	double totalUniformityDeviation = 0.0; // uniformity value
	size_t deviationCheck = 0;
	size_t devCount = 0; // count of deviations

	for (const size_t& cellCount : nodeCount) {
		// if cell has nodes, calculate deviation
		if (cellCount > 0) {
			// calculate deviation from ideal uniformity
			deviationCheck = fabs(cellCount - idealUniformity);
			// if deviation is 0 or 1, skip it
			// as it is considered uniform enough
			if (deviationCheck == 0 || deviationCheck == 1) {
				continue;
			} else { //
				totalUniformityDeviation += deviationCheck;
				devCount++;
			}
		}
	}

	double worstUniformityDeviation =
			gridCount * (numOfNodes - idealUniformity); // every cell has nodes off by all nodes
	totalUniformityDeviation /= static_cast<double>(gridCount);
	double nodeUniformityRatio = 1.0 - (totalUniformityDeviation / worstUniformityDeviation);

	return nodeUniformityRatio; // return node uniformity ratio
}

double LayoutStatistics::edgeOrthogonality(const GraphAttributes& ga) {
	const Graph& mainGraph = ga.constGraph();
	size_t numOfEdges = mainGraph.numberOfEdges(); // for efficiency

	// if no edge exists, return 0.0
	if (numOfEdges < 1) {
		return 0.0;
	}

	double meanEdgeOrthogonality = 0.0; // edge orthogonality value
	double angleTemp = 0.0;

	for (const edge& e : mainGraph.edges) {
		// Calculate angles for each edge to the x-axis
		double x = ga.x(e->source()) - ga.x(e->target()); // x coordinate difference
		double y = ga.y(e->source()) - ga.y(e->target()); // y coordinate difference

		// calculate angle relative to x-axis
		angleTemp = atan2(y, x) * 180 / M_PI; // convert radians to degrees

		// applying mod 90 to get the angle in the range [0, 90)
		angleTemp = fmod(angleTemp, 90.0);

		// checking if angle is closer to x-axis
		if (angleTemp > 45) {
			angleTemp = 90.0 - angleTemp; // adjust angle to closer axis
		}
		meanEdgeOrthogonality += angleTemp; // add absolute angle to mean edge orthogonality
	}

	// mean edge orthogonality deviation: sum of all angles divided by number of edges
	return meanEdgeOrthogonality /= static_cast<double>(numOfEdges);
}

std::pair<double, double> LayoutStatistics::centerOfMass(const GraphAttributes& ga) {
	const Graph& mainGraph = ga.constGraph();
	size_t n = mainGraph.numberOfNodes();
	if (n < 2) {
		if (n == 1) {
			return std::make_pair(ga.x(*mainGraph.nodes.begin()),
					ga.x(*mainGraph.nodes.begin())); // return only node as coordinate
		} else { // no nodes, no center of mass
			return std::make_pair(0.0, 0.0); // return (0.0, 0.0) as center of mass
		}
	}

	double sumX = 0.0; // sum of x coordinates
	double sumY = 0.0; // sum of y coordinates

	// Summing up x and y coordinates of all nodes
	for (const node& u : mainGraph.nodes) {
		sumX += ga.x(u);
		sumY += ga.y(u);
	}
	sumX /= n; // average x coordinate
	sumY /= n; // average y coordinate

	return std::make_pair(sumX, sumY); // return center of mass node
}

double LayoutStatistics::closestPairOfPoints(const GraphAttributes& ga) {
	const Graph& mainGraph = ga.constGraph();
	size_t n = mainGraph.numberOfNodes();

	if (n < 2) {
		return -1.0;
	}

	ArrayBuffer<std::pair<std::pair<node, node>, double>> allDistances;
	distancesBetweenAllNodes(ga, allDistances);
	double smallestDist = std::numeric_limits<double>::max(); // for smallest distance

	// Extracting smallest distance from allDistances, e.g. smallest distance between two nodes
	for (const auto& pair : allDistances) {
		if (smallestDist > pair.second) {
			smallestDist = pair.second; // update smallest distance
		}
	}

	return smallestDist; // return smallest distance
}

std::pair<std::pair<double, double>, std::pair<double, double>> LayoutStatistics::borderCoordinates(
		const GraphAttributes& ga) {
	const Graph& mainGraph = ga.constGraph();
	if (mainGraph.numberOfNodes() == 0) {
		return std::make_pair(std::make_pair(-1.0, -1.0),
				std::make_pair(-1.0, -1.0)); // no nodes, no border coordinates
	}

	// Get min and max coordinates via DRect
	DRect borderCoords = ga.GraphAttributes::boundingBox();
	double minX = borderCoords.p1().m_x;
	double maxX = borderCoords.p2().m_x;
	double minY = borderCoords.p1().m_y;
	double maxY = borderCoords.p2().m_y;

	return std::make_pair(std::make_pair(minX, maxX), std::make_pair(minY, maxY));
}

double LayoutStatistics::horizontalVerticalBalance(const GraphAttributes& ga, const bool vertical) {
	const Graph& mainGraph = ga.constGraph();
	size_t n = mainGraph.numberOfNodes();

	if (n == 0) {
		return -1.0; // no nodes, no horizontal balance
	}
	if (n <= 2) {
		return 0.0; // one or two nodes are always balanced
	}

	// Initialize min and max variables
	double min = 0.0, max = 0.0;

	auto borderCoords = borderCoordinates(ga);
	if (!vertical) { // horizontal balance
		min = borderCoords.first.first; // min x-coordinate
		max = borderCoords.first.second; // max x-coordinate
	} else { // vertical balance
		min = borderCoords.second.first; // min y-coordinate
		max = borderCoords.second.second; // max y-coordinate
	}

	const double center = (min + max) / 2.0; // center x-coordinate
	double leftTopCount = 0.0; // sum of coordinates of nodes on the left side
	double rightBottomCount = 0.0; // sum of x-coordinates of nodes on the right side

	// Counting nodes on left and right side of center
	for (const node& u : mainGraph.nodes) {
		// init maxVal based on either vertical or horizontal balance
		double maxVal = vertical ? ga.y(u) : ga.x(u);

		// Counting nodes on each side.
		// Excluding nodes exactly on the center line, as they don't contribute to balance
		if (maxVal < center) {
			leftTopCount += maxVal; // add to left side sum
		} else if (maxVal > center) {
			rightBottomCount += maxVal; // add to right side sum
		}
	}
	// Calculate and return balance ratio
	if (leftTopCount > rightBottomCount) {
		return rightBottomCount /= leftTopCount;
	} else {
		return leftTopCount /= rightBottomCount;
	}
}

double LayoutStatistics::nodeOrthogonality(const GraphAttributes& ga, const double epsilon) {
	const Graph& mainGraph = ga.constGraph();
	size_t n = mainGraph.numberOfNodes();

	if (n == 0) {
		return -1.0; // no nodes, no orthogonality
	}
	size_t nodeOrthogonalityCount = 0; // count of orthogonal nodes

	double x = 0.0, y = 0.0; // coordinates of first node
	for (const node& u : mainGraph.nodes) {
		x = ga.x(u);
		y = ga.y(u);

		if (std::fabs(std::fmod(x + y, 1.0)) < epsilon) {
			// if x + y is not an integer, node is not orthogonal
			nodeOrthogonalityCount++; // increment orthogonal node count
		}
	}
	return (static_cast<double>(nodeOrthogonalityCount) / n)
			* 100.0; // return node orthogonality percentage
}

double LayoutStatistics::averageFlow(const GraphAttributes& ga) {
	const Graph& mainGraph = ga.constGraph();
	size_t m = mainGraph.numberOfEdges(); // for efficiency

	if (m == 0 || !ga.directed()) {
		return -1.0; // if graph has no edge, or isn't directed, return -1.0
	}
	double totalAngle = 0.0;
	for (const edge& e : mainGraph.edges) {
		node source = e->source(); // source node of edge
		node target = e->target(); // target node of edge

		// calculate edge angle
		double x = ga.x(source) - ga.x(target); // x coordinate difference
		double y = ga.y(source) - ga.y(target); // y coordinate difference
		totalAngle += atan2(y, x) * 180 / M_PI; // convert to degrees
	}
	return totalAngle /= m; // return average flow
}

double LayoutStatistics::upwardsFlow(const GraphAttributes& ga) {
	const Graph& mainGraph = ga.constGraph();
	size_t m = mainGraph.numberOfEdges(); // for efficiency

	if (m == 0 || !ga.directed()) {
		return -1.0; // if graph has no edge, or isn't directed, return -1.0
	}
	size_t upwardsFlowCount = 0;
	for (const edge& e : mainGraph.edges) {
		// If y-coordinate of source node is smaller than of target node,
		// it points upwards
		if (ga.y(e->source()) < ga.y(e->target())) {
			upwardsFlowCount++;
		}
	}
	return (static_cast<double>(upwardsFlowCount) / m) * 100.0; // return upwards flow percentage
}

double LayoutStatistics::concentration(const GraphAttributes& ga) {
	const Graph& mainGraph = ga.constGraph();
	size_t n = mainGraph.numberOfNodes(); // for efficiency

	if (n == 0) { // No concentration possible
		return -1.0;
	}
	if (n == 1) {
		return 0.0; // Single node has zero variance
	}

	const std::pair<double, double> center = centerOfMass(ga);
	const double centerX = center.first; // x coordinate of center of mass
	const double centerY = center.second; // y coordinate of center of mass
	double sumOfDistancesToCenter = 0.0;
	std::vector<double> distances;
	distances.reserve(n); // size n

	// Calculate sum of distances from all nodes to center of mass
	for (const node& u : mainGraph.nodes) {
		// Calculating distances for each node to center of mass
		const double dx = ga.x(u) - centerX;
		const double dy = ga.y(u) - centerY;
		const double dist = std::hypot(dx, dy); // hypot = sqrt(dx^2 + dy^2)

		distances.push_back(dist); // adding dist to vector
		sumOfDistancesToCenter += dist; // add distance to sum
	}

	// mean dist to center of mass of all nodes
	const double meanDist = sumOfDistancesToCenter / static_cast<double>(n);

	// Calculating variance of distances to center of mass
	double sumSqrdDiff = 0.0; // sum of squared differences from mean dist
	for (const double& dist : distances) {
		const double diff = dist - meanDist;
		sumSqrdDiff += diff * diff;
	}

	// returning variance center of mass
	return sumSqrdDiff / static_cast<double>(n);
}
}
