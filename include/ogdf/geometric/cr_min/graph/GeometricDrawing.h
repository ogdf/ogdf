/** \file
 *
 * \author Marcel Radermacher
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

#ifdef OGDF_INCLUDE_CGAL

#	include <ogdf/basic/EdgeArray.h>
#	include <ogdf/basic/GraphAttributes.h>
#	include <ogdf/fileformats/GraphIO.h>
#	include <ogdf/fileformats/SvgPrinter.h>
#	include <ogdf/geometric/cr_min/datastructure/OGDFVector.h>
#	include <ogdf/geometric/cr_min/geometry/objects/LineSegment.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Point.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Rectangle.h>
#	include <ogdf/geometric/cr_min/graph/Drawing.h>
#	include <ogdf/geometric/cr_min/graph/OGDFFaceWrapper.h>
#	include <ogdf/geometric/cr_min/graph/OGDFGraphWrapper.h>
#	include <ogdf/geometric/cr_min/tools/ogdf/Converter.h>
#	include <ogdf/geometric/cr_min/tools/ogdf/Universal.h>

#	include <map>
#	include <tuple>
#	include <vector>

namespace ogdf {
namespace internal {
namespace gcm {
namespace graph {

/*! represents a geometrical graph, i.e., each vertex has a coordinate/ point
 */

template<typename Kernel_, typename Graph_>
class GeometricDrawing : public Drawing<Kernel_, Graph_> {
public:
	using Kernel = Kernel_;
	using Graph = Graph_;
	using Node = typename Graph::Node;
	using Edge = typename Graph::Edge;

	using Point = geometry::Point_t<Kernel>;
	using LineSegment = geometry::LineSegment_t<Kernel>;

private:
	Graph& g;

	using PointToNodeMap = std::map<Point, Node>;
	using PointToNodeMapIterator = typename PointToNodeMap::iterator;

	std::vector<Edge> visited_edges;

	PointToNodeMap point_to_node;

	datastructure::NodeVector<Point, Graph> node_to_point;

	OGDFFaceWrapper m_outer_face;
	bool outer_face_initialized {false};

public:
	explicit GeometricDrawing(Graph& _g) : Drawing<Kernel_, Graph>(_g), g(_g), node_to_point(g) {
		/*nothing to do*/
	}

	inline void set_outer_face(OGDFFaceWrapper ef) {
		m_outer_face = ef;
		outer_face_initialized = true;
	}

	inline OGDFFaceWrapper& outer_face() {
		OGDF_ASSERT(outer_face_initialized);
		return m_outer_face;
	}

	inline const OGDFFaceWrapper& outer_face() const {
		OGDF_ASSERT(outer_face_initialized);
		return m_outer_face;
	}

	inline void clear() {
		point_to_node.clear();
		node_to_point.clear();
		outer_face_initialized = false;
	}

	size_t number_of_points() const { return node_to_point.size(); }

	const std::vector<Point>& points() const { return node_to_point; }

	// does not keep the consitency with add node!
	inline void set_point(const Node& v, const Point& p) {
		OGDF_ASSERT((size_t)v->index() <= g.max_node_index());
		if (g.max_node_index() >= node_to_point.size()) {
			node_to_point.resize(g.max_node_index() + 1);
		}
		node_to_point[v] = p;
	}

	inline Point get_point(const Node& v) const {
		OGDF_ASSERT((size_t)v->index() < node_to_point.size());
		return node_to_point[v];
	}

	inline Point& operator[](const Node& v) {
		OGDF_ASSERT((size_t)v->index() < node_to_point.size());
		return node_to_point[v];
	}

	inline Point operator[](const Node& v) const { return get_point(v); }

	inline LineSegment get_segment(const Edge& e) const {
		return {get_point(e->source()), get_point(e->target())};
	}

	inline geometry::Bbox bbox() const {
		double xmin = std::numeric_limits<double>::infinity();
		double ymin = std::numeric_limits<double>::infinity();
		double xmax = -xmin;
		double ymax = -ymin;

		geometry::Bbox bb(xmin, ymin, xmax, ymax);
		for (Node v : g.nodes()) {
			bb += get_point(v).bbox();
		}
		return bb;
	}
};

template<typename Kernel, typename Graph, typename EdgeList>
void generate_graph_from_list(const std::vector<geometry::Point_t<Kernel>>& nodes,
		const EdgeList& edge_list, std::vector<typename Graph::Edge>& map_edge,
		GeometricDrawing<Kernel, Graph>& d) {
	using Node = typename Graph::Node;
	using Edge = typename Graph::Edge;
	using Point = typename geometry::Point_t<Kernel>;

	std::vector<Node> node_map;
	std::vector<Point> map;
	Graph& g = d.get_graph();
	for (const Point& p : nodes) {
		node_map.push_back(g.add_node());
		d.set_point(node_map.back(), p); //TODO check ?
	}
	for (const std::tuple<unsigned int, unsigned int, unsigned int>& edge : edge_list) {
		if (node_map[std::get<0>(edge)] != node_map[std::get<1>(edge)]) {
			const Edge e = g.add_edge(node_map[std::get<0>(edge)], node_map[std::get<1>(edge)]);
			map_edge.push_back(e);
		} else {
			map_edge.push_back(nullptr);
		}
	}
#	ifdef OGDF_DEBUG
	for (const auto e : map_edge) {
		OGDF_ASSERT(!e || e->graphOf() == &g.get_ogdf_graph());
	}
#	endif
}

template<typename Kernel, typename Graph, typename EdgeList>
void generate_graph_from_list(std::map<unsigned int, geometry::Point_t<Kernel>>& nodes,
		const EdgeList& edge_list, std::vector<typename Graph::Edge>& map_edge,
		std::map<unsigned int, typename Graph::Node>& node_map, GeometricDrawing<Kernel, Graph>& d) {
	using Edge = typename Graph::Edge;

	Graph& g = d.get_graph();
	for (const std::tuple<unsigned int, unsigned int, unsigned int>& edge : edge_list) {
		unsigned int source = std::get<0>(edge);
		unsigned int target = std::get<1>(edge);


		if (node_map.find(source) == node_map.end()) {
			auto u = g.add_node();
			d.set_point(u, nodes[source]);
			node_map.insert(std::make_pair(source, u));
		}
		if (node_map.find(target) == node_map.end()) {
			auto u = g.add_node();
			d.set_point(u, nodes[target]);
			node_map.insert(std::make_pair(target, u));
		}
	}

	for (const std::tuple<unsigned int, unsigned int, unsigned int>& edge : edge_list) {
		if (node_map[std::get<0>(edge)] != node_map[std::get<1>(edge)]) {
			const Edge e = g.add_edge(node_map[std::get<0>(edge)], node_map[std::get<1>(edge)]);
			OGDF_ASSERT(e != nullptr);
			map_edge.push_back(e);
		} else {
			map_edge.push_back(nullptr);
		}
	}

#	ifdef OGDF_DEBUG
	for (const auto e : map_edge) {
		OGDF_ASSERT(!e || e->graphOf() == &g.get_ogdf_graph());
	}
#	endif
}

template<typename Kernel, typename Graph, typename EdgeList>
void generate_graph_from_list(std::map<unsigned int, geometry::Point_t<Kernel>>& nodes,
		const EdgeList& edge_list, std::vector<typename Graph::Edge>& map_edge,
		GeometricDrawing<Kernel, Graph>& g) {
	std::map<unsigned int, typename Graph::Node> node_map;
	generate_graph_from_list(nodes, edge_list, map_edge, node_map, g);
}

template<typename Kernel, typename Graph>
std::tuple<datastructure::NodeVector<typename Graph::Node, Graph>,
		datastructure::EdgeVector<typename Graph::Edge, Graph>>
copy(GeometricDrawing<Kernel, Graph>& original, GeometricDrawing<Kernel, Graph>& copy) {
	using Node = typename Graph::Node;
	using Edge = typename Graph::Edge;

	datastructure::NodeVector<Node, Graph> node_map(original);
	datastructure::EdgeVector<Edge, Graph> edge_map(original);

	for (auto v : original.nodes()) {
		node_map[v] = copy.add_node(original.get_point(v));
	}

	for (auto e : original.edges()) {
		edge_map[e] = copy.add_edge(node_map[e->source()], node_map[e->target()]);
	}

	return std::make_tuple(node_map, edge_map);
}

template<typename Kernel, typename Graph>
void ogdf_attributes_to_geometric_drawing(const GraphAttributes& ga,
		GeometricDrawing<Kernel, Graph>& d) {
	OGDF_ASSERT(&d.get_graph().get_ogdf_graph() == &ga.constGraph());
	for (auto v : d.get_graph().nodes()) {
		geometry::Point_t<Kernel> p(ga.x(v), ga.y(v));
		d.set_point(v, p);
	}
}

template<typename Kernel, typename Graph>
void geometric_drawing_to_ogdf_attributes(const GeometricDrawing<Kernel, Graph>& d,
		GraphAttributes& ga) {
	OGDF_ASSERT(&d.get_graph().get_ogdf_graph() == &ga.constGraph());
	for (auto v : d.get_graph().nodes()) {
		auto p = d.get_point(v);
		ga.x(v) = CGAL::to_double(p.x());
		ga.y(v) = CGAL::to_double(p.y());
	}
}

}
}
}
}

#endif
