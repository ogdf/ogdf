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

#	include <ogdf/geometric/cr_min/geometry/objects/LineSegment.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace geometry {

template<typename kernel>
class Window_t : public LineSegment_t<kernel> {
private:
	using LineSegment = LineSegment_t<kernel>;

public:
	unsigned int source_id;
	unsigned int target_id;
	using LineSegment::LineSegment; //inherit constructors
};

template<typename kernel>
std::ostream& operator<<(std::ostream& os, const Window_t<kernel>& v) {
	os << "vector[" << v.source() << ", " << v.target() << "]";
	return os;
}

} //namespace
}
}
}

#endif
