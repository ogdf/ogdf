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

#	include <ogdf/geometric/cr_min/geometry/algorithm/CollinearTriple.h>
#	include <ogdf/geometric/cr_min/geometry/algorithm/LargestCircleInPolygon.h>
#	include <ogdf/geometric/cr_min/geometry/algorithm/RandomPointInPolygon.h>
#	include <ogdf/geometric/cr_min/graph/geometric_crossing_min/CrossingMinimalRegion.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace graph {

template<typename Kernel, typename Graph>
class CrossingMinimalPosition {
private:
	using Point = geometry::Point_t<Kernel>;
	using Segment = geometry::LineSegment_t<Kernel>;
	using Polygon = geometry::Polygon_t<Kernel>;
	using Drawing = graph::GeometricDrawing<Kernel, Graph>;
	using Node = typename Graph::Node;
	using Edge = typename Graph::Edge;

public:
	static Point compute(const Drawing& d, const Node& v, const std::vector<Edge>& sample,
			geometry::Rectangle_t<Kernel>& rect_box) {
		unsigned int dump_a, dump_b;
		auto region = CrossingMinimalRegion<Kernel, Graph>::compute(d, v, sample, rect_box, dump_a,
				dump_b);
		if (geometry::is_clockwise(region)) {
			region = geometry::reverse(region);
		}

		Point p = d.get_point(v);
		if (region.size() > 2) {
			if (!region.is_convex()) {
				p = geometry::largest_circle_in_polygon(region, 1e-5);
			} else {
				p = geometry::centroid(region);
			}
		}

		//round
		OGDF_ASSERT(rect_box.has_on_bounded_side(p));
		p = Point(CGAL::to_double(p.x()), CGAL::to_double(p.y()));
		OGDF_ASSERT(rect_box.has_on_bounded_side(p));
		return p;
	}
};

}
}
}
}

#endif
