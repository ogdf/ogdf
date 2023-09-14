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

#	include <ogdf/lib/mapbox/mapbox_earcut.h>

#	include <array>
#	include <vector>

#	include <CGAL/number_utils.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace geometry {

template<typename Kernel>
std::vector<unsigned int> triangulation(const geometry::Polygon_t<Kernel>& polygon) {
	using MB_Point = std::array<double, 2>;

	if (polygon.size() > 3) {
		std::vector<std::vector<MB_Point>> mb_polygon;
		mb_polygon.resize(1);
		for (unsigned int i = 0; i < polygon.size(); ++i) {
			mb_polygon[0].push_back(
					{CGAL::to_double(polygon[i].x()), CGAL::to_double(polygon[i].y())});
		}

		return mapbox::earcut<unsigned int>(mb_polygon);
	} else {
		return {0, 1, 2};
	}
}

}
}
}
}

#endif
