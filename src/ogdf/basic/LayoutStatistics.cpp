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
#include <cstddef>
#include <iterator>
#include <limits>
#include <utility>
#include <vector>

namespace ogdf {

double graphArea(const GraphAttributes& ga) {
	if (ga.constGraph().numberOfNodes() < 2) {
		return 0.0;
	}

	const DRect bBox = ga.GraphAttributes::boundingBox();
	return bBox.height() * bBox.width();
}

EdgeArray<double> LayoutStatistics::edgeLengths(const GraphAttributes& ga, bool considerSelfLoops) {
	EdgeArray<double> values(ga.constGraph(), 0.0);
	if (ga.constGraph().numberOfEdges() < 1) {
		return values;
	}

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
		values[e] = len;
	}

	return values;
}

EdgeArray<size_t> LayoutStatistics::numberOfBends(const GraphAttributes& ga, bool considerSelfLoops) {
	const Graph& G = ga.constGraph();
	EdgeArray<size_t> values(G, 0);

	if (G.numberOfEdges() == 0) {
		return values;
	}

	for (const edge& e : G.edges) {
		if (considerSelfLoops || !e->isSelfLoop()) {
			values[e] = ga.bends(e).size();
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

			angles.pushBack(atan2(ey - vy, ex - vx));
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
				double psi1 = atan2(p1.m_y - by, p1.m_x - bx);

				const DPoint& p2 = *it.succ();
				double psi2 = atan2(p2.m_y - by, p2.m_x - bx);

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

EdgeArray<size_t> LayoutStatistics::numberOfCrossings(const GraphAttributes& ga) {
	const Graph& G = ga.constGraph();
	EdgeArray<size_t> crossings(G, 0);

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
				if (!adj->isSource()) {
					continue;
				}
				edge e = adj->theEdge();
				edge eOrig = origEdge[e];

				// Ignore original edges incident to vOrig.
				if (eOrig->source() != vOrig && eOrig->target() != vOrig) {
					crossings[eOrig] += (k - 1);
				}
			}
		}
	}
	return crossings;
}

EdgeArray<size_t> LayoutStatistics::numberOfNodeCrossings(const GraphAttributes& ga) {
	const Graph& G = ga.constGraph();
	if (G.numberOfEdges() == 0) {
		EdgeArray<size_t> values(G, 0);
		return values;
	}

	EdgeArray<size_t> values(G, 0);
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
		values[e] = nCrossingsE;
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
	EdgeArray<size_t> crossings = LayoutStatistics::numberOfCrossings(ga);
	// ArrayBuffer<int> maxCrossings; // Initialize
	size_t sumOfCrossings = 0;

	// Calculate sum of actual crossings for all edges
	for (auto e : G.edges) {
		sumOfCrossings += crossings[e];
	}

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
	double ratio = static_cast<double>(sumOfCrossings) / static_cast<double>(sumMaxCrossings);
	return ratio * 100.0; // Return percentage
}

NodeArray<size_t> LayoutStatistics::numberOfNodeOverlaps(const GraphAttributes& ga) {
	const Graph& G = ga.constGraph();
	NodeArray<size_t> values(G, 0);

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
		values[v] = nOverlapsV;
	}

	return values;
}

double LayoutStatistics::edgeLengthDeviation(const GraphAttributes& ga, EdgeArray<double>& out) {
	const Graph& mainGraph = ga.constGraph();
	const size_t m = mainGraph.numberOfEdges();
	if (m < 2) {
		out.init(mainGraph, 0.0);
		return 0.0; // no deviation with less than two edges
	}

	out.init(mainGraph, 0.0); // init out to label the correct graph
	EdgeArray<double> lengths(mainGraph, true);
	for (const edge& e : mainGraph.edges) {
		node u = e->source();
		node v = e->target();
		double length = hypot(ga.x(u) - ga.x(v), ga.y(u) - ga.y(v));
		lengths[e] = (length > 0.0) ? length : 0.0;
	}

	double edgeSum = 0.0;

	double delta = 0.0;

	for (edge e : mainGraph.edges) {
		const double len = lengths[e];
		out[e] = len;
		edgeSum += len;
	}

	// if edgesNum (amount of edges) > 0, then edgeSum / edgeCount, else avgEdgeLen = 0.0
	double avgEdgeLen = edgeSum / static_cast<double>(m);

	double edgeLenDevVal = 0.0; // for one edge deviation value, the lower the better

	// For all edges, subtract all deltas of edge lengths to avg.
	for (const edge& e : ga.constGraph().edges) {
		// calculating delta for each edge using std::fabs, for floating numbers
		delta = std::fabs(out[e] - avgEdgeLen);
		// adding delta to edge length deviation sum
		edgeLenDevVal += (delta);
		out[e] = delta;
	}

	// calculating average edge length deviation
	return edgeLenDevVal /= static_cast<double>(m);
}

NodeArray<double> LayoutStatistics::neighbourhoodPreservation(const GraphAttributes& ga) {
	const Graph& mainGraph = ga.constGraph();
	size_t n = mainGraph.numberOfNodes(); // for efficiency

	// array of preservation values, of each node/vertex
	NodeArray<double> nodePreservations(mainGraph, 0.0);

	if (n < 2) {
		if (n == 1) {
			nodePreservations[mainGraph.firstNode()] = 1.0;
		}
		return nodePreservations;
	}

	double nodePreservationValue = 0.0; // preservation value for each node, higher = better
	// Uncomment 1. | double graphPreservationValue = 0.0; // preservation value for whole graph, higher = better


	// Distance NodeArray-Matrix n x n: vector of (distance(u, v), v)
	NodeArray<std::vector<std::pair<double, node>>> distM(mainGraph);

	// calculating distances between all nodes {(u, v), (v, u)}
	for (const node u : mainGraph.nodes) {
		auto& rowU = distM[u];
		rowU.reserve(n);

		for (const node v : mainGraph.nodes) {
			if (u == v) {
				distM[u].emplace_back(0.0, v);
				continue;
			}
			const double dist = std::hypot(ga.x(u) - ga.x(v), ga.y(u) - ga.y(v));
			distM[u].emplace_back(dist, v);
		}

		// sorting all distances for each node within the distance NodeArray-Matrix
		std::sort(rowU.begin(), rowU.end(),
				[](const std::pair<double, node>& a, const std::pair<double, node>& b) {
					return a.first < b.first;
				});
	}

	size_t nodeDegree = 0; // current node degree

	// for each node, calculate preservation value
	for (node u : mainGraph.nodes) {
		nodeDegree = u->degree(); // set current node degree
		const auto& rowU = distM[u];

		if (nodeDegree == 0 || rowU.size() < 2) {
			// if has no neighbors, skip and leave value at 0.0
			continue;
		}

		// limiting max degree if input graph has self-loops or multi-edges
		const size_t maxDegree = std::min(nodeDegree, rowU.size() - 1);
		if (maxDegree == 0) {
			continue;
		}


		nodePreservationValue = 0.0; // reseting value for each node/round

		// for |deg(currentNode)|, calculate preservation value for current node
		for (size_t i = 1; i <= maxDegree; ++i) {
			const node& neighbor = distM[u][i].second; // closest neighbor node to current node n

			// checking if neighbor is connected in graph
			if (mainGraph.searchEdge(u, neighbor) != nullptr) {
				nodePreservationValue += 1;
				continue;
			}
		}
		nodePreservations[u] = nodePreservationValue / static_cast<double>(maxDegree);
		// Uncomment 1. | graphPreservationValue += nodePreservationValue / static_cast<double>(nodeDegree);
	}

	// final graph preservation value calculation
	// Uncomment 1. | graphPreservationValue /= n;
	return nodePreservations;
}

NodeArray<double> LayoutStatistics::gabrielRatio(const GraphAttributes& ga,
		Graph& gabrielGraphReference) {
	const Graph& mainGraph = ga.constGraph();
	gabrielGraphReference.clear();

	NodeArray<double> nodeGabrielRatios(mainGraph, 0.0);

	if (mainGraph.numberOfNodes() < 3) {
		return nodeGabrielRatios;
	}


	// size_t numOfNodes = mainGraph.numberOfNodes(); // for efficiency

	// array for per-node Gabriel ratios
	NodeArray<size_t> gabrielCount(mainGraph, 0); // count Gabriel-edges for each node
	NodeArray<size_t> degreeCount(mainGraph, 0); // count degree of each node

	// create map from original nodes to new gabrielGraph nodes
	NodeArray<node> nodeMap(mainGraph, nullptr);
	for (const node& n : mainGraph.nodes) {
		nodeMap[n] = gabrielGraphReference.newNode(); // assigns node to node map
	}

	// for each edge, calculate Gabriel edge
	for (edge e : mainGraph.edges) {
		node u = e->source(); // first node of edge
		node v = e->target(); // second node of edge

		degreeCount[u]++;
		degreeCount[v]++;

		// calculating midpoints (between the two nodes) both of x and y, and the radius range of the distance
		double midPointX = (ga.x(u) + ga.x(v)) * 0.5;
		double midPointY = (ga.y(u) + ga.y(v)) * 0.5;

		double dx = ga.x(u) - ga.x(v);
		double dy = ga.y(u) - ga.y(v);

		// instead of sqrt (hypot for distToMidPoint), we square both sides (distToMidpoint^2 < radius^2) to avoid sqrt
		// which turns radius = distance / 2 into: radius^2 = (distance^2) / 4; * 0.25 for efficiency
		double radiusSquared = (dx * dx + dy * dy) * 0.25;


		bool isGabrielEdge = true;
		for (const node& w : mainGraph.nodes) {
			if (w == v || w == u) { // skip same nodes
				continue;
			}
			// calculating distance to midpoint
			double dx2 = ga.x(w) - midPointX;
			double dy2 = ga.y(w) - midPointY;
			// avoiding sqrt by squaring both sides (radius and distToMidPoint)
			double distToMidPointSquared = dx2 * dx2 + dy2 * dy2; // dx^2 + dy^2

			// if distance to midpoint is less than radius => not a Gabriel edge
			if (distToMidPointSquared < radiusSquared) {
				isGabrielEdge = false;
				break;
			}
		}

		if (isGabrielEdge) {
			gabrielGraphReference.newEdge(nodeMap[u], nodeMap[v]);
			// incrementing Gabriel count for both nodes (u, v)
			gabrielCount[u]++;
			gabrielCount[v]++;
		}
	}

	// calculating Gabriel ratio for each node
	for (const node& n : mainGraph.nodes) {
		int incident = (n->degree());
		if (incident) {
			nodeGabrielRatios[n] =
					static_cast<double>(gabrielCount[n]) / static_cast<double>(incident);
		} else {
			nodeGabrielRatios[n] = 0.0;
		}
	}

	// calculating Gabriel ratio for whole graph
	/* Uncomment 3. | double graphGabrielRatio = 0.0;
	for (size_t i = 0; i < numOfNodes; ++i) {
		graphGabrielRatio += nodeGabrielRatios[i];
	}
	graphGabrielRatio /= numOfNodes; // final division for Gabriel Ratio of graph
	*/

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

	double minDist = std::numeric_limits<double>::max(); // for min. distance between two nodes
	double maxDist = 0.0; // for max distance between two nodes

	// for (auto u : mainGraph.nodes) {
	// 	for (auto v : mainGraph.nodes) {
	for (node u : mainGraph.nodes) {
		for (node v = u->succ(); v; v = v->succ()) {
			if (u == v) {
				continue; // skip same nodes
			}

			const double dist = std::hypot(ga.x(u) - ga.x(v), ga.y(u) - ga.y(v));
			if (dist < minDist) {
				minDist = dist; // update min. distance
			}
			if (dist > maxDist) {
				maxDist = dist; // update max. distance
			}
		}
	}

	// prevents division by zero, else calculates and returns Node Resolution (NR)
	return maxDist > 0.0 ? (minDist / maxDist) : 0.0;
}

double LayoutStatistics::angularResolution(const GraphAttributes& ga) {
	const Graph& mainGraph = ga.constGraph();

	// if graph has less than 3 nodes, there cannot exist two edges
	if (mainGraph.numberOfNodes() < 3) {
		return 0.0;
	}

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
			angles.pushBack(Math::radiansToDegrees(atan2(y, x))); // convert to degrees
		}

		// Sort angles in ascending order
		angles.quicksort();

		angles.pushBack(angles.front() + 2 * Math::pi); // for wrapping around to first angle at the end

		// Ideal angle calculation
		double perfectAngle = 2 * Math::pi / static_cast<double>(degree); // optimal angle in degrees
		double nodeDev = 0.0;

		// Calculate the angle deviation for current node
		for (ListConstIterator<double> it = angles.begin(); it.succ().valid();
				++it) { // checks if next angle exists (last angle is the first angle we added)

			// Calc angle deviation: angle2 - angle1
			double deviationGap = *(it.succ()) - *it;
			nodeDev += std::fabs(deviationGap - perfectAngle); // add absolute deviation to nodeDev
		}

		totalAngleDev +=
				nodeDev / static_cast<double>(degree); // average angle deviation for current node
		++nodesCount; // increment count for nodes with degree greater than 1
	}

	// Return mean angle deviation ratio (or 0.0 if no nodes with degree greater than 1)
	return (nodesCount > 0) ? (totalAngleDev / static_cast<double>(nodesCount)) : 0.0;
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

	const DRect bBox = ga.GraphAttributes::boundingBox();

	// width and height of the bounding box
	const double width = bBox.width(); // width of bounding box
	const double height = bBox.height(); // height of bounding box

	// check for infinite values
	if (!std::isfinite(width) || !std::isfinite(height)) {
		return 0.0;
	}

	const double epsilon = std::numeric_limits<double>::epsilon() * 1e3;
	// if width or height is <=0,
	// return 0.0 to avoid division by zero, negative values, or extremely small values
	if (width <= epsilon || height <= epsilon) {
		return 0.0;
	}

	return width / height; // aspect ratio
}

double LayoutStatistics::nodeUniformity(const GraphAttributes& ga, size_t gridWidth,
		size_t gridHeight) {
	const Graph& mainGraph = ga.constGraph();
	size_t numOfNodes = mainGraph.numberOfNodes(); // for efficiency

	if (gridWidth == 0 || gridHeight == 0) {
		return 0.0;
	}

	// if graph has less than 2 nodes, node uniformity is trivial
	if (numOfNodes < 2) {
		return 0.0;
	}

	auto bBox = ga.GraphAttributes::boundingBox();
	// Get min x- and y-coordinates via DRect
	double minX = bBox.p1().m_x;
	double minY = bBox.p1().m_y;

	// width and height of the bounding box
	// adding 1 for minimum cell size, in case minX == maxX or minY == maxY (e.g. points could be on the same line)
	double width = bBox.width() + 1.0; // width of bounding box
	double height = bBox.height() + 1.0; // height of bounding box

	// if width or height is <=0,
	// return 0.0 to avoid division by zero and negative values
	if (!std::isfinite(width) || !std::isfinite(height) || width <= 0.0 || height <= 0.0) {
		return 0.0;
	}

	// size of each cell in the grid
	double cellWidth = width / static_cast<double>(gridWidth); // width of each cell
	double cellHeight = height / static_cast<double>(gridHeight); // height of each cell

	size_t gridCount = gridWidth * gridHeight; // total number of cells in the grid

	if (gridCount == 0 || !std::isfinite(gridCount)) // if gridCount = 0, return 0.0
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
		size_t cellIndex = cellX + cellY * gridWidth;
		if (!std::isfinite(cellIndex) || cellIndex < 0) {
			continue;
		}
		if (cellIndex >= nodeCount.size()) {
			continue;
		}

		nodeCount[cellIndex]++; // cellY * gridWidth because cells are stored in row-major order
	}

	// calculate ideal uniformity
	size_t idealUniformity = ceil(static_cast<double>(numOfNodes) / static_cast<double>(gridCount));

	double totalUniformityDeviation = 0.0; // uniformity value
	size_t deviationCheck = 0;

	for (const size_t& cellCount : nodeCount) {
		// if cell has nodes, calculate deviation
		if (cellCount > 0) {
			// calculate deviation from ideal uniformity
			deviationCheck = cellCount > idealUniformity ? cellCount - idealUniformity
														 : idealUniformity - cellCount;
			// if deviation is 0 or 1, skip it
			// as it is considered uniform enough
			if (deviationCheck == 0 || deviationCheck == 1) {
				continue;
			} else {
				totalUniformityDeviation += deviationCheck;
			}
		}
	}

	// every cell has nodes off by all nodes
	double worstUniformityDeviation = 0.0;
	long long diffWorst = numOfNodes - idealUniformity;
	if (diffWorst > 0) {
		worstUniformityDeviation = static_cast<double>(gridCount) * static_cast<double>(diffWorst);
	} else {
		return 0.0;
	}


	totalUniformityDeviation /= static_cast<double>(gridCount);

	if (worstUniformityDeviation <= 0.0) {
		return 0.0;
	}

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
		// store source and target nodes of edge
		const node u = e->source();
		const node v = e->target();

		// Calculate angles for each edge to the x-axis
		double x = ga.x(u) - ga.x(v); // x coordinate difference
		double y = ga.y(u) - ga.y(v); // y coordinate difference

		// calculate angle relative to x-axis
		angleTemp = Math::radiansToDegrees(atan2(y, x)); // convert radians to degrees

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

DPoint LayoutStatistics::centerOfMass(const GraphAttributes& ga) {
	const Graph& mainGraph = ga.constGraph();
	size_t n = mainGraph.numberOfNodes();

	if (n == 0) { // no nodes, no center
		return DPoint(0.0, 0.0); // return (0.0, 0.0) as center of mass
	}
	if (n == 1) {
		const node singleNode = mainGraph.firstNode();
		const double x = ga.x(singleNode);
		const double y = ga.y(singleNode);
		if (!std::isfinite(x) || !std::isfinite(y)) {
			std::cout << "LayoutStatistics::centerOfMass: single node has non-finite cords.\n";
			return DPoint(0.0, 0.0);
		}

		return DPoint(x, y); // return only node as coordinate
	}

	double sumX = 0.0; // sum of x coordinates
	double sumY = 0.0; // sum of y coordinates
	size_t invalidCordsCount = 0;

	// Summing up x and y coordinates of all nodes
	for (const node& u : mainGraph.nodes) {
		const double x = ga.x(u);
		const double y = ga.y(u);

		if (!std::isfinite(x) || !std::isfinite(y)) {
			std::cout << "LayoutStatistics::centerOfMass: non-finite cords node found.\n";
			invalidCordsCount++;
			continue;
		}

		sumX += x;
		sumY += y;
	}
	if (invalidCordsCount == 0) {
		return DPoint(-1.0, -1.0);
	}

	n -= invalidCordsCount;

	sumX /= n; // average x coordinate
	sumY /= n; // average y coordinate

	return DPoint(sumX, sumY); // return center of mass node
}

double LayoutStatistics::closestPairOfPoints(const GraphAttributes& ga) {
	const Graph& mainGraph = ga.constGraph();
	std::vector<node> validNodes;

	for (const node& u : mainGraph.nodes) {
		if (std::isfinite(ga.x(u)) && std::isfinite(ga.y(u))) {
			validNodes.push_back(u);
		}
	}
	size_t validNodeCount = validNodes.size();
	if (validNodeCount < mainGraph.numberOfNodes()) {
		std::cout << "LayoutStatistics::closestPairOfPoints: Ignoring "
				  << (mainGraph.numberOfNodes() - validNodeCount)
				  << " nodes with non-finite coordinates.\n";
	}

	if (validNodeCount < 2) {
		return -1.0;
	}


	double smallestDist = std::numeric_limits<double>::infinity(); // for smallest distance
	double distance = 0.0;

	// calculating smallest distance from between two nodes of graph
	for (size_t i = 0; i < validNodeCount; ++i) {
		const node u = validNodes[i];
		for (size_t j = i + 1; j < validNodeCount; ++j) {
			const node v = validNodes[j];
			// calculating euclidean distance between two nodes
			distance = std::hypot(ga.x(u) - ga.x(v), ga.y(u) - ga.y(v));
			if (distance < smallestDist) {
				smallestDist = distance; // update smallest distance
			}
		}
	}
	return smallestDist; // return smallest distance
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

	const DRect bBox = ga.GraphAttributes::boundingBox();

	if (!vertical) { // horizontal balance
		min = bBox.p1().m_x; // min x-coordinate
		max = bBox.p2().m_x; // max x-coordinate
	} else { // vertical balance
		min = bBox.p1().m_y; // min y-coordinate
		max = bBox.p2().m_y; // max y-coordinate
	}

	const double center = (min + max) * 0.5; // center x-coordinate
	double leftTopCount = 0.0; // sum of coordinates of nodes on the left side
	double rightBottomCount = 0.0; // sum of x-coordinates of nodes on the right side

	// Counting nodes on left and right side of center
	for (const node& u : mainGraph.nodes) {
		// init maxVal based on either vertical or horizontal balance
		double maxVal = vertical ? ga.y(u) : ga.x(u);

		// Counting nodes on each side.
		// Excluding nodes exactly on the center line, as they don't contribute to balance
		if (maxVal < center) {
			leftTopCount += std::fabs(center - maxVal); // add to left side sum
		} else if (maxVal > center) {
			rightBottomCount += std::fabs(maxVal - center); // add to right side sum
		}
	}
	// If both sides have zero sum (all on center)
	if (leftTopCount == 0.0 && rightBottomCount == 0.0) {
		return 0.0;
	}

	// Calculate and return balance ratio
	const double maxSide = std::max(leftTopCount, rightBottomCount);
	const double minSide = std::min(leftTopCount, rightBottomCount);

	// If all on center, as latter case differs slightly
	if (maxSide + minSide == 0.0) {
		return 0.0;
	}

	// 0 = perfectly balanced, 1 = maximally unbalanced
	return (maxSide - minSide) / (maxSide + minSide);
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
		totalAngle += Math::radiansToDegrees(atan2(y, x)); // convert to degrees
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

	const DPoint center = centerOfMass(ga);
	const double centerX = center.m_x; // x coordinate of center of mass
	const double centerY = center.m_y; // y coordinate of center of mass
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
