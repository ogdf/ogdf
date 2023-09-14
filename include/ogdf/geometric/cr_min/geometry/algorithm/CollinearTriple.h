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

#	include <ogdf/geometric/cr_min/geometry/objects/Point.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Vector.h>
#	include <ogdf/geometric/cr_min/graph/GeometricDrawing.h>

#	include <vector>

namespace ogdf {
namespace internal {
namespace gcm {
namespace geometry {

template<typename Kernel>
bool has_collinear_pair(const std::vector<Point_t<Kernel>>& points, const unsigned int v) {
	std::vector<typename Kernel::FT> angles;
	geometry::Vector_t<Kernel> x(1, 0);
	auto p_v = points[v];
	for (unsigned int w = 0; w < points.size(); w++) {
		if (w != v) {
			auto p_w = points[w];
			if (p_w == p_v) {
				return true;
			}
			for (unsigned int u = w + 1; u < points.size(); u++) {
				if (u != v) {
					auto p_u = points[u];
					if (CGAL::collinear(p_u, p_v, p_w) || p_u == p_v) {
						return true;
					}
				}
			}
		}
	}
	return false;
}

template<typename Kernel, typename Graph>
bool has_collinear_pair(const graph::GeometricDrawing<Kernel, Graph>& d,
		const typename Graph::Node v) {
	std::vector<typename Kernel::FT> angles;
	geometry::Vector_t<Kernel> x(1, 0);
	auto p_v = d.get_point(v);
	for (auto w : d.get_graph().nodes()) {
		if (w != v) {
			auto p_w = d.get_point(w);
			if (p_w == p_v) {
				return true;
			}
			for (auto u : d.get_graph().nodes()) {
				if (u != v && u != w) {
					auto p_u = d.get_point(u);
					if (CGAL::collinear(p_u, p_v, p_w) || p_u == p_v) {
						return true;
					}
				}
			}
		}
	}
	return false;
}

template<typename Kernel, typename Graph>
bool has_collinear_triple(const graph::GeometricDrawing<Kernel, Graph>& d) {
	for (auto v : d.get_graph().nodes()) {
		if (has_collinear_pair(d, v)) {
			return true;
		}
	}
	return false;
}

template<typename Kernel>
bool has_collinear_triple(const std::vector<Point_t<Kernel>>& points) {
	for (unsigned int v = 0; v < points.size(); ++v) {
		if (has_collinear_pair(points, v)) {
			return true;
		}
	}
	return false;
}

template<typename Kernel>
void resolve_collinearity(std::vector<Point_t<Kernel>>& points) {
	if (!geometry::has_collinear_triple(points)) {
		return;
	}

	std::mt19937_64 gen(ogdf::randomSeed());

	std::uniform_real_distribution<double> distr(-1, 1);

	for (unsigned int v = 0; v < points.size(); ++v) {
		bool col = geometry::has_collinear_pair(points, v);
		while (col) {
			auto prev = points[v];
			auto x = prev.x() + distr(gen);
			auto y = prev.y() + distr(gen);
			points[v] = Point_t<Kernel>(x, y);
			col = geometry::has_collinear_pair(points, v);
			if (col) {
				points[v] = prev;
			}
		}
		if (col) {
			break;
		}
	}
}

}
}
}
}

#endif
