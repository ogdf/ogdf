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

#	include <ogdf/geometric/cr_min/geometry/objects/Direction.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Vector.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace geometry {

/*
 * A cone through the origin defined by two vectors. Everything that is counter clockwise in between the left and right
 * is within the cone. left and right are not supposed to be parallel
 * A cone is not necessarly acute / convex
 */
template<typename Kernel>
class Cone_t {
private:
	using Vector = Vector_t<Kernel>;
	using Direction = Direction_t<Kernel>;
	Direction m_left;
	Direction m_right;

public:
	Cone_t(const Vector& left, const Vector& right) : Cone_t(left.direction, right.direction()) {
		// nothing to do
	}

	Cone_t(const Direction& left, const Direction& right) : m_left(left), m_right(right) {
		OGDF_ASSERT(left != right);
		// nothing to do
	}

	const Direction& left() const { return m_left; }

	const Direction& right() const { return m_right; }
};

/*
 * is c_1 subset of c_2
 */
template<typename Kernel>
bool is_subset(const Cone_t<Kernel>& c_1, const Cone_t<Kernel>& c_2) {
	return c_1.left().counterclockwise_in_between(c_2.left(), c_2.right())
			&& c_1.right().counterclockwise_in_between(c_2.left(), c_2.right());
}

template<typename Kernel>
bool is_in(const Direction_t<Kernel>& d, const Cone_t<Kernel>& c) {
	return d.counterclockwise_in_between(c.left(), c.right());
}

template<typename Kernel>
bool do_intersect(const Cone_t<Kernel>& c_1, const Cone_t<Kernel>& c_2) {
	return is_in(c_1.left(), c_2) || is_in(c_1.right(), c_2) || is_in(c_2.left(), c_1)
			|| is_in(c_2.right(), c_1);
}

template<typename kernel>
std::ostream& operator<<(std::ostream& os, const Cone_t<kernel>& c) {
	os << "cone[(" << c.left().vector() << ", " << c.right().vector() << ")]";
	return os;
}

}
}
}
}

#endif
