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

#	include <ogdf/geometric/cr_min/geometry/objects/Circle.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Cone.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Direction.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Line.h>
#	include <ogdf/geometric/cr_min/geometry/objects/LineSegment.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Point.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Polygon.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Polyline.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Ray.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Rectangle.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Vector.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace geometry {

template<typename _Kernel>
class Geometry_t {
public:
	using Kernel = _Kernel;
	using Point = Point_t<Kernel>;
	using LineSegment = LineSegment_t<Kernel>;
	using Ray = Ray_t<Kernel>;
	using Polygon = Polygon_t<Kernel>;
	using Line = Line_t<Kernel>;
	using Polyline = Polyline_t<Kernel>;
	using Circle = Circle_t<Kernel>;
	using Direction = Direction_t<Kernel>;
	using Cone = Cone_t<Kernel>;
	using Rectangle = Rectangle_t<Kernel>;
	using Vector = Vector_t<Kernel>;
};

}
}
}
}

#endif
