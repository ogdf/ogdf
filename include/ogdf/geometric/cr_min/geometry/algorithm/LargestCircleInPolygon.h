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

#	include <ogdf/basic/PriorityQueue.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Polygon.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Rectangle.h>

#	include <queue>

//https://blog.mapbox.com/a-new-algorithm-for-finding-a-visual-center-of-a-polygon-7c77e6492fbc
namespace ogdf {
namespace internal {
namespace gcm {
namespace geometry {

template<typename Kernel>
Point_t<Kernel> largest_circle_in_polygon(const Polygon_t<Kernel>& polygon, double precision = 1) {
	using FT = typename Kernel::FT;

	struct Cell {
		Bbox bb;
		FT m_distance;

		Cell(Bbox& bb_, FT distance) : bb(bb_), m_distance(distance) {
			//nothing to do
		}

		FT cell_size() const {
			return std::min(bb.height(), bb.width()); //use max?
		}

		FT potential() const { return m_distance + bb.width() * bb.height() / 4; }

		bool operator<(const Cell& c) const { return potential() > c.potential(); }
	};

	Bbox bb_p(polygon.bbox());

	PriorityQueue<Cell> pq;

	Cell opt(bb_p, squared_distance(polygon, bb_p.template center<Kernel>()));
	pq.push(opt);

	while (!pq.empty()) {
		const Cell top = pq.top();
		pq.pop();

		if (top.potential() < 0) {
			continue;
		}

		if (top.m_distance > opt.m_distance) {
			opt = top;
		}

		if (top.potential() - opt.m_distance < precision) {
			continue;
		}

		auto h = top.bb.height() / 2;
		auto w = top.bb.width() / 2;

		Bbox b1(top.bb.xmin(), top.bb.ymin(), top.bb.xmin() + w, top.bb.ymin() + h);
		Bbox b2(top.bb.xmin() + w, top.bb.ymin(), top.bb.xmax(), top.bb.ymin() + h);
		Bbox b3(top.bb.xmin() + w, top.bb.ymin() + h, top.bb.xmax(), top.bb.ymax());
		Bbox b4(top.bb.xmin(), top.bb.ymin() + h, top.bb.xmin() + w, top.bb.ymax());

		pq.push(Cell(b1, squared_distance(polygon, b1.template center<Kernel>())));
		pq.push(Cell(b2, squared_distance(polygon, b2.template center<Kernel>())));
		pq.push(Cell(b3, squared_distance(polygon, b3.template center<Kernel>())));
		pq.push(Cell(b4, squared_distance(polygon, b4.template center<Kernel>())));
	}

	// center point must be in polygon
	OGDF_ASSERT(contains(polygon, opt.bb.template center<Kernel>()));

	return opt.bb.template center<Kernel>();
}

}
}
}
}

#endif
