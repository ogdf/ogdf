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

#	include <ogdf/geometric/cr_min/tools/math.h>

#	include <iostream>
#	include <limits>

#	include <CGAL/Aff_transformation_2.h>
#	include <CGAL/Cartesian.h>
#	include <CGAL/Vector_2.h>
#	include <CGAL/aff_transformation_tags.h>
#	include <math.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace geometry {
using namespace tools;

template<typename kernel>
using Vector_t = CGAL::Vector_2<kernel>;

template<typename kernel>
inline bool is_the_same(const Vector_t<kernel>& v1, const Vector_t<kernel>& v2) {
	return v1.x() == v2.x() && v1.y() == v2.y();
}

template<typename kernel>
inline Vector_t<kernel> rotate(const Vector_t<kernel>& v, const double angle) {
	CGAL::Aff_transformation_2<kernel> rotation(CGAL::ROTATION, std::sin(angle), std::cos(angle));
	return std::move(rotation(v));
}

template<typename kernel>
inline typename kernel::FT dot(const Vector_t<kernel>& v1, const Vector_t<kernel>& v2) {
	return v1 * v2;
}

template<typename kernel>
inline typename kernel::FT cross(const Vector_t<kernel>& v1, const Vector_t<kernel>& v2) {
	return CGAL::determinant(v1, v2);
}

template<typename kernel>
inline bool turn(const Vector_t<kernel>& v1, const Vector_t<kernel>& v2) {
	return cross(v1, v2) >= 0;
}

template<typename kernel>
inline bool left_turn(const Vector_t<kernel>& v1, const Vector_t<kernel>& v2) {
	return CGAL::orientation(v1, v2) == CGAL::LEFT_TURN;
}

template<typename kernel>
inline bool right_turn(const Vector_t<kernel>& v1, const Vector_t<kernel>& v2) {
	return CGAL::orientation(v1, v2) == CGAL::RIGHT_TURN;
}

template<typename kernel>
inline bool parallel(const Vector_t<kernel>& v1, const Vector_t<kernel>& v2) {
	return CGAL::orientation(v1, v2) == CGAL::COLLINEAR;
}

template<typename kernel>
inline Vector_t<kernel> normalize(const Vector_t<kernel>& v) {
	OGDF_ASSERT(v.squared_length() > 0);
	return v * CGAL::inverse(tools::approx_sqrt(v.squared_length()));
}

template<typename kernel>
inline Vector_t<kernel> normal(const Vector_t<kernel>& v) {
	return v.perpendicular(CGAL::POSITIVE);
}

template<typename kernel>
inline Vector_t<kernel> bisect(const Vector_t<kernel> v1, const Vector_t<kernel>& v2) {
	if ((-v1).direction() == v2.direction()) {
		return normalize(normal(v1));
	}
	return (normalize(v1) + normalize(v2)) * 0.5;
}

template<typename kernel>
inline typename kernel::FT cos_angle(const Vector_t<kernel>& v1, const Vector_t<kernel>& v2) {
	OGDF_ASSERT(!isZero(v1.squared_length()) && !isZero(v2.squared_length()));
	typename kernel::FT v = v1 * v2;
	return v * CGAL::inverse(tools::approx_sqrt(v1.squared_length() * v2.squared_length()));
}

template<typename kernel>
inline typename kernel::FT angle(const Vector_t<kernel>& v1, const Vector_t<kernel>& v2) {
	return tools::acos(cos_angle(v1, v2));
}

//ccw
template<typename kernel>
inline typename kernel::FT full_angle(const Vector_t<kernel>& v1, const Vector_t<kernel>& v2) {
	using type = typename kernel::FT;
	const type alpha = geometry::angle(v1, v2);
	const bool left = left_turn(v1, v2);
	return (type)left * alpha + (type)(1.0 - left) * ((type)2.0 * tools::const_pi<type>() - alpha);
}

template<typename kernel>
inline bool is_valid(const Vector_t<kernel>& v) {
	return !isinf(v.x()) && !isinf(v.y()) && !isnan(v.x()) && !isnan(v.y());
}

template<typename kernel>
void serialize(const Vector_t<kernel>& v, std::ostream& os) {
	os.write(reinterpret_cast<const char*>(&v.x()), sizeof(v.x()));
	os.write(reinterpret_cast<const char*>(&v.y()), sizeof(v.y()));
}

template<typename kernel>
void deserialize(Vector_t<kernel>& v, std::istream& in) {
	in.read(reinterpret_cast<char*>(&v.x()), sizeof(v.x()));
	in.read(reinterpret_cast<char*>(&v.y()), sizeof(v.y()));
}

template<typename Vector>
struct VectorExactLess {
	bool operator()(const Vector& a, const Vector& b) const {
		return a.x() < b.x() || (a.x() == b.x() && a.y() < b.y());
	}
};

} //namespace
}
}
}

#endif
