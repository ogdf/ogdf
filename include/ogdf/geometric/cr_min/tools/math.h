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

#	include <CGAL/CORE_Expr.h>
#	include <CGAL/Gmpfr.h>
#	include <CGAL/Gmpq.h>
#	include <CGAL/gmpxx.h>
#	include <mpfr.h>

#	ifndef OGDF_GEOMETRIC_MPFR_NUMBER_OF_DIGITS

#		define OGDF_GEOMETRIC_MPFR_NUMBER_OF_DIGITS 77 // equals 256 bits

#	endif

#	ifndef OGDF_GEOMETRIC_MPFR_EQUALITY_THRESHOLD
#		define OGDF_GEOMETRIC_MPFR_EQUALITY_THRESHOLD 1e-20
#	endif


namespace ogdf {
namespace internal {
namespace gcm {
namespace tools {

template<typename T>
inline double cast(const T& a) {
	//return a.to_double();
	return CGAL::to_double(a);
}

template<>
inline double cast(const double& a) {
	return a;
}

inline mpfr_rnd_t std_rnd_to_mpfr_rnd(const std::float_round_style& e) {
	switch (e) {
	case std::round_indeterminate:
		OGDF_ASSERT(false);
		return MPFR_RNDZ;
	case std::round_toward_zero:
		return MPFR_RNDZ;
	case std::round_to_nearest:
		return MPFR_RNDN;
	case std::round_toward_infinity:
		return MPFR_RNDU;
	case std::round_toward_neg_infinity:
		return MPFR_RNDD;
	default:
		OGDF_ASSERT(false);
		return MPFR_RNDZ;
	}
}

template<typename T>
inline const T const_pi() {
	return CGAL_PI;
}

template<>
inline const CGAL::Gmpfr const_pi() {
	CGAL::Gmpfr y;
	mpfr_const_pi(y.fr(), std_rnd_to_mpfr_rnd(CGAL::Gmpfr::get_default_rndmode()));
	return y;
}

template<typename T>
inline const T approx_sqrt(const T& v) {
	return (T)sqrt(CGAL::to_double(v));
}

inline const mpz_class approx_sqrt(const mpz_class& v) { return sqrt(v); }

inline const mpq_class approx_sqrt(const mpq_class& v) {
	return mpq_class(approx_sqrt(v.get_num()), approx_sqrt(v.get_den()));
}

inline const CGAL::Gmpfr approx_sqrt(const CGAL::Gmpfr& v) {
	CGAL::Gmpfr y;
	mpfr_sqrt(y.fr(), v.fr(), std_rnd_to_mpfr_rnd(CGAL::Gmpfr::get_default_rndmode()));
	return y;
}

inline const CORE::Expr approx_sqrt(const CORE::Expr& v) { return CGAL::sqrt(v); }

template<typename T>
inline const T acos(const T& v) {
	return (T)std::acos(CGAL::to_double(v));
}

template<>
inline const CGAL::Gmpfr acos(const CGAL::Gmpfr& v) {
	CGAL::Gmpfr y;
	mpfr_acos(y.fr(), v.fr(), std_rnd_to_mpfr_rnd(CGAL::Gmpfr::get_default_rndmode()));
	return y;
}

template<>
inline const CORE::Expr acos(const CORE::Expr& v) {
	return std::acos(CGAL::to_double(v));
}

// DIGITS


template<typename t>
inline bool isEqual(const t& a, const t& b) {
	return a == b;
}

template<>
inline bool isEqual(const double& a, const double& b) {
	return fabs(a - b) < OGDF_GEOMETRIC_MPFR_EQUALITY_THRESHOLD;
}

template<>
inline bool isEqual(const CGAL::Gmpfr& a, const CGAL::Gmpfr& b) {
	return abs(a - b) < OGDF_GEOMETRIC_MPFR_EQUALITY_THRESHOLD;
}

template<typename T>
inline bool isLessEqual(const T& a, const T& b) {
	return (a < b) || isEqual(a, b);
}

template<typename T>
inline bool isLess(const T& a, const T& b) {
	return (a < b) && !isEqual(a, b);
}

template<typename T>
inline bool isGreaterEqual(const T& a, const T& b) {
	return (a > b) || isEqual(a, b);
}

template<typename T>
inline bool isGreater(const T& a, const T& b) {
	return (a > b) && !isEqual(a, b);
}

template<typename T>
inline bool isZero(const T& a) {
	return isEqual(a, (T)0.0);
}

}
} // namespace

inline std::ostream& operator<<(std::ostream& os, const CGAL::Gmpq& p) {
	os.precision(20);
	os << CGAL::to_double(p);
	return os;
}

inline std::ostream& operator<<(std::ostream& os, const CORE::Expr& p) {
	os.precision(20);
	os << CGAL::to_double(p);
	return os;
}

}
}

#endif
