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

#	include <ogdf/geometric/cr_min/geometry/objects/Line.h>
#	include <ogdf/geometric/cr_min/geometry/objects/LineSegment.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Point.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Ray.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Vector.h>
#	include <ogdf/geometric/cr_min/tools/math.h>

#	include <vector>

#	include <CGAL/Bbox_2.h>
#	include <CGAL/Iso_rectangle_2.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace geometry {

template<typename kernel>
using Rectangle_t = CGAL::Iso_rectangle_2<kernel>;

class Bbox : public CGAL::Bbox_2 {
private:
	using Bbox_ = CGAL::Bbox_2;

public:
	using Bbox_::Bbox_;

	Bbox() {
		// nothing to do
	}

	Bbox(const Bbox_& b) : Bbox_(b) {
		//nothing to do;
	}

	inline double width() const { return Bbox::xmax() - Bbox::xmin(); }

	inline double height() const { return Bbox::ymax() - Bbox::ymin(); }

	inline double area() const { return width() * height(); }

	template<typename Kernel>
	inline Point_t<Kernel> center() const {
		return {xmin() + width() / 2, ymin() + height() / 2};
	}
};

inline Bbox equalize(const Bbox& bb) {
	double max = std::max(bb.width(), bb.height());
	return {bb.xmin(), bb.ymin(), bb.xmin() + max, bb.ymin() + max};
}

template<typename kernel>
inline Rectangle_t<kernel> scale_up(const Rectangle_t<kernel>& rect, const typename kernel::FT v) {
	const typename kernel::FT s = std::min(rect.xmax() - rect.xmin(), rect.ymax() - rect.ymin()) * v;
	const Vector_t<kernel> t(s, s);
	return {rect.min() - t, rect.max() + t};
}
} // namespace
}
}
}

#endif
