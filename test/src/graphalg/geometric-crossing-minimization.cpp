/** \file
 * \brief Basic test for geometric crossing minimization
 *
 * \author Kassian Köck
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

#include <ogdf/basic/LayoutStatistics.h>
#include <ogdf/energybased/StressMinimization.h>
#include <ogdf/geometric/CrossingMinimalPosition.h>
#include <ogdf/geometric/GeometricEdgeInsertion.h>
#include <ogdf/geometric/GeometricVertexInsertion.h>
#include <ogdf/geometric/VertexMovement.h>
#include <ogdf/geometric/VertexOrder.h>
#include <ogdf/planarity/MaximalPlanarSubgraphSimple.h>

#include <resources.h>

#ifdef OGDF_INCLUDE_CGAL

#	include <CGAL/Random.h>

#endif

enum method { EDGE_INSERTION = 1, VERTEX_INSERTION = 2, VERTEX_MOVEMENT = 3, NONE = 0 };

void readGML(const ResourceFile* file, Graph& g, GraphAttributes& ga) {
	std::stringstream is {file->data()};
	AssertThat(GraphIO::readGML(ga, g, is), IsTrue());
	AssertThat(g.empty(), IsFalse());
	AssertThat(g.numberOfEdges(), IsGreaterThan(0));
}

void removeDegreeOne(Graph& g) {
	bool had_degree_one;
	do {
		had_degree_one = false;
		safeForEach(g.nodes, [&](node v) {
			if (v->degree() <= 1) {
				had_degree_one = true;
				g.delNode(v);
			}
		});

	} while (had_degree_one);
}

void edgeInsertion(Graph& g, GraphAttributes& ga) {
#ifdef OGDF_INCLUDE_CGAL
	//compute planar subgraph
	List<edge> edges_to_hide;
	MaximalPlanarSubgraphSimple<int> planar_sg;
	planar_sg(g, edges_to_hide);

	//setup and perform geometric edge insertion
	StressMinimization sm;
	sm.hasInitialLayout(false);

	CrossingMinimalPositionFast pos;
	pos.setExactComputation();
	int dim = g.numberOfNodes() * 200;
	pos.setBoundingBox(-dim, -dim, dim, dim);

	GeometricEdgeInsertion gei(g);
	gei.setInitialLayouter(&sm);
	gei.setVertexPosition(&pos);
	gei.setHiddenEdgeSet(&edges_to_hide);
	gei(ga);
#endif
}

void vertexInsertion(Graph& g, GraphAttributes& ga) {
#ifdef OGDF_INCLUDE_CGAL
	StressMinimization sm;
	sm(ga);

	//compute vertex order
	auto first = g.firstNode();
	double x_min = ga.x(first), x_max = ga.x(first), y_min = ga.y(first), y_max = ga.y(first);

	//compute bounding box
	std::vector<node> vertices;
	for (auto v = g.firstNode(); v; v = v->succ()) {
		x_min = std::min(x_min, ga.x(v));
		x_max = std::max(x_max, ga.x(v));
		y_min = std::min(y_min, ga.y(v));
		y_max = std::max(y_max, ga.y(v));
		ga.width(v) = 0;
		ga.height(v) = 0;
		vertices.push_back(v);
	}

	std::sort(vertices.begin(), vertices.end(),
			[&](node a, node b) { return a->degree() > b->degree(); });

	List<node> vertex_order;
	for (auto v : vertices) {
		vertex_order.pushBack(v);
	}

	//setup and perform vertex insertion
	CrossingMinimalPositionFast pos;
	pos.setExactComputation();
	pos.setBoundingBox(x_min - 1000, y_min - 1000, x_max + 1000, y_max + 1000);

	GeometricVertexInsertion vi(g);
	vi.setVertexPosition(&pos);
	vi.setVertexOrder(&vertex_order);
	vi(ga);
#endif
}

void vertexMovement(Graph& g, GraphAttributes& ga) {
#ifdef OGDF_INCLUDE_CGAL
	StressMinimization sm;
	sm(ga);

	CrossingMinimalPositionPrecise pos;
	pos.setExactComputation();

	VertexMovement vm;

	CrossingVertexOrder vo(ga, OrderEnum::asc, MeasureEnum::squared);
	auto vertex_order = vo.get_vertex_order();

	auto first = vertex_order.front();
	double x_min = ga.x(first), x_max = ga.x(first), y_min = ga.y(first), y_max = ga.y(first);

	for (auto v : vertex_order) {
		x_min = std::min(x_min, ga.x(v));
		x_max = std::max(x_max, ga.x(v));
		y_min = std::min(y_min, ga.y(v));
		y_max = std::max(y_max, ga.y(v));

		ga.width(v) = 0;
		ga.height(v) = 0;
	}

	pos.setBoundingBox(x_min - 1000, y_min - 1000, x_max + 1000, y_max + 1000);
	vm.setVertexPosition(&pos);
	vm.setVertexOrder(&vertex_order);
	vm(ga);
#endif
}

int expectedNumberOfCrossings(method m, const ResourceFile* file) {
	std::stringstream is {file->data()};
	std::string line;
	do {
		is >> line;
	} while (line.find("Crossing_Number:") == std::string::npos);
	for (int i = 0; i < m; ++i) {
		is >> line;
	}
	int crossingNumber;
	is >> crossingNumber;
	return crossingNumber;
}

int numberOfCrossings(GraphAttributes& ga) {
	int crossingNumber = 0;
	auto buffer = LayoutStatistics::numberOfCrossings(ga);
	for (const auto& item : buffer) {
		crossingNumber += item;
	}
	return crossingNumber / 2;
}

void testGeometricCrossingMinimization(method m) {
	std::string description = "Geometric Crossing Minimization via ";
	switch (m) {
	case EDGE_INSERTION:
		description += "Edge Insertion";
		break;
	case VERTEX_INSERTION:
		description += "Vertex Insertion";
		break;
	case VERTEX_MOVEMENT:
		description += "Vertex Movement";
		break;
	case NONE:
		break;
	}
#ifdef OGDF_INCLUDE_CGAL
	describe(description, [m]() {
#else
	describe_skip(description, [m]() {
#endif
		for_each_file("geometry", [&](const ResourceFile* file) {
			it("works on " + file->fullPath(), [&] {
				setSeed(0);
#ifdef OGDF_INCLUDE_CGAL
				CGAL::get_default_random() = CGAL::Random(randomSeed());
#endif
				Graph g;
				GraphAttributes ga(g, GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics);
				readGML(file, g, ga);
				AssertThat(numberOfCrossings(ga), Equals(expectedNumberOfCrossings(NONE, file)));


				removeDegreeOne(g);
				switch (m) {
				case EDGE_INSERTION:
					edgeInsertion(g, ga);
					break;
				case VERTEX_INSERTION:
					vertexInsertion(g, ga);
					break;
				case VERTEX_MOVEMENT:
					vertexMovement(g, ga);
					break;
				case NONE:
					break;
				}

#ifdef OGDF_GEOMETRIC_CR_MIN_DEBUG
				string filePrefix =
						file->name().substr(0, file->name().size() - 4) + "_" + to_string(m);
				GraphIO::write(ga, filePrefix + ".gml");
				GraphIO::write(ga, filePrefix + ".svg");
#endif
				AssertThat(numberOfCrossings(ga),
						IsLessThanOrEqualTo(expectedNumberOfCrossings(m, file) * 1.5));
			});
		});
	});
}

go_bandit([] {
	testGeometricCrossingMinimization(EDGE_INSERTION);
	testGeometricCrossingMinimization(VERTEX_INSERTION);
	testGeometricCrossingMinimization(VERTEX_MOVEMENT);
});
