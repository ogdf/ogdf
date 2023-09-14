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

#	include <ogdf/geometric/cr_min/geometry/objects/Polygon.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace geometry {

/** Peturb vertex points p of polygons by a small epsilon along the bisector of the two edges incident to p
 *
 */
template<typename Kernel>
Polygon_t<Kernel> peturb(const Polygon_t<Kernel>& polygon, double eps = 1e-5) {
	using Segment = LineSegment_t<Kernel>;

	auto bisector = [&](const Segment& s1, const Segment& s2) {
		const auto& v_in = s1.to_vector();
		const auto& v_out = s2.to_vector();
		auto bisected = bisect(-v_in, v_out);
		if (geometry::right_turn(v_in, bisected)) {
			bisected = -bisected;
		}

		return bisected;
	};

	Polygon_t<Kernel> perturbed;
	auto itr = polygon.edges_circulator();
	for (unsigned int i = 0; i < polygon.size(); ++i) {
		perturbed.push_back(polygon[i] + bisector(*(itr - 1), *itr) * eps);
		++itr;
	}
	return perturbed;
}

}
}
}
}

#endif
