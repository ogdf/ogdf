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

#	include <ogdf/geometric/cr_min/geometry/objects/Polyline.h>
#	include <ogdf/geometric/cr_min/graph/GeometricDrawing.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace graph {

template<typename Kernel_, typename Graph>
class PolylineDrawing
	: public GeometricDrawing<Kernel_, Graph> { //TODO a polyline drawing is not really a geometric drawing :-/
public:
	using Kernel = Kernel_;
	using Node = node;
	using Edge = edge;

	using Point = geometry::Point_t<Kernel>;
	using LineSegment = geometry::LineSegment_t<Kernel>;
	using Polyline = geometry::Polyline_t<Kernel>;


private:
	using parent = GeometricDrawing<Kernel, Graph>;
	datastructure::EdgeVector<Polyline, Graph> edge_shape;

public:
	PolylineDrawing(Graph& g) : GeometricDrawing<Kernel, Graph>(g) { }

	void clear() {
		edge_shape.clear();
		parent::clear();
	}

	inline const Polyline& get_polyline(const Edge& e) const { return edge_shape[e]; }

	inline void set_polyline(const Edge& e, Polyline& p) {
		OGDF_ASSERT(!p.is_degenerate());
		OGDF_ASSERT(p.front() == this->get_point(e->source()));
		OGDF_ASSERT(p.back() == this->get_point(e->target()));
		edge_shape[e] = p;
	}

	inline void set_polyline(const Edge& e, Polyline&& p) {
		OGDF_ASSERT(!p.is_degenerate());
		OGDF_ASSERT(p.front() == this->get_point(e->source()));
		OGDF_ASSERT(p.back() == this->get_point(e->target()));

		edge_shape[e] = std::move(p);
	}

	inline geometry::Bbox bbox() const {
		double xmin = std::numeric_limits<double>::infinity();
		double ymin = std::numeric_limits<double>::infinity();
		double xmax = -xmin;
		double ymax = -ymin;

		geometry::Bbox bb(xmin, ymin, xmax, ymax);
		for (Node v : parent::nodes()) {
			bb += parent::get_point(v).bbox();
		}

		for (Edge e : parent::edges()) {
			bb += edge_shape[e].bbox();
		}

		return bb;
	}
};

}
}
}
}

#endif
