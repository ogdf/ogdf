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

#	include <ogdf/geometric/cr_min/geometry/algorithm/MapBoxTriangulation.h>
#	include <ogdf/geometric/cr_min/geometry/algorithm/RestrictedTriangulation.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Point.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Polygon.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Vector.h>

#	include <random>
#	include <vector>

#	include <CGAL/Triangle_2.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace geometry {


/*! Compute a random point within a polygon
 */
template<typename Kernel>
Point_t<Kernel> random_point_in_polygon(const Polygon_t<Kernel>& polygon, std::mt19937_64& gen) {
	auto triang = triangulation(polygon);
	std::vector<double> weights; // index * 3 yields index in triang
	std::vector<double> interval;
	interval.push_back(0);

	for (unsigned int i = 0; i < triang.size(); i = i + 3) {
		unsigned int a = triang[i];
		unsigned int b = triang[i + 1];
		unsigned int c = triang[i + 2];

		CGAL::Triangle_2<Kernel> t(polygon[a], polygon[b], polygon[c]);

		weights.push_back(CGAL::to_double(CGAL::abs(t.area())));

		interval.push_back(interval.back() + 1);
	}

	// select triangle in polygon proportional to weight
	std::piecewise_constant_distribution<double> dist(interval.begin(), interval.end(),
			weights.begin());

	// barycenter coordinate at random	-> random point in triangle
	std::uniform_real_distribution<double> t_dist(0.0, 1.0);
	double w_1 = t_dist(gen);
	double w_2 = t_dist(gen);
	double w_3 = t_dist(gen);
	double n_w = w_1 + w_2 + w_3;

	using Vector = Vector_t<Kernel>;
	using Point = Point_t<Kernel>;
	Point p(0, 0);


	unsigned tri_id = std::floor(dist(gen)) * 3;

	unsigned int a = triang[tri_id];
	unsigned int b = triang[tri_id + 1];
	unsigned int c = triang[tri_id + 2];
	Vector v_1(p, polygon[a]);
	Vector v_2(p, polygon[b]);
	Vector v_3(p, polygon[c]);

	Vector v = (v_1 * w_1 + v_2 * w_2 + v_3 * w_3) * (1.0 / n_w);
	Point r(v.x(), v.y());
	return r;
}


}
}
}
}

#endif
