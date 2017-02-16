/** \file
 * \brief Tests for Math.h
 *
 * \author Ivo Hedtke
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

#include <bandit/bandit.h>
#include <ogdf/basic/Math.h>
#include <ogdf/basic/EpsilonTest.h>

using namespace bandit;
using namespace ogdf;

template<typename T>
static void testGcdAndLcm(const char *type)
{
	it(string("computes gcd of large numbers of type ") + string(type), []() {
		T big = numeric_limits<T>::max();
		AssertThat(Math::gcd(big, big), Equals(big));
	});
	it(string("computes lcm of large numbers of type ") + string(type), []() {
		T big = numeric_limits<T>::max();
		AssertThat(Math::lcm(big, big), Equals(big));
	});
}

static void testHarmonic()
{
	it("computes harmonic numbers correctly", []() {
		EpsilonTest eps;
		AssertThat(eps.equal(Math::harmonic(0), 1.0), IsTrue());
		AssertThat(eps.equal(Math::harmonic(1), 1.0), IsTrue());
		AssertThat(eps.equal(Math::harmonic(2), 1.5), IsTrue());
		AssertThat(eps.equal(Math::harmonic(3), 1.5 + 1/3.0), IsTrue());
		AssertThat(Math::harmonic(10), IsLessThan(3));
		AssertThat(Math::harmonic(11), IsGreaterThan(3));
		AssertThat(Math::harmonic(30), IsLessThan(4));
		AssertThat(Math::harmonic(31), IsGreaterThan(4));
		AssertThat(Math::harmonic(82), IsLessThan(5));
		AssertThat(Math::harmonic(83), IsGreaterThan(5));
		AssertThat(Math::harmonic(12366), IsLessThan(10));
		AssertThat(Math::harmonic(12367), IsGreaterThan(10));
	});
	it("computes huge harmonic numbers correctly", []() {
		unsigned i = 2012783313;
		double result;
		while ((result = Math::harmonic(i)) < 22.0) {
			i++;
		}
		AssertThat(i, Equals(2012783315u));
	});
}

go_bandit([]() {
	describe("Math.h", []() {
		it("computes gcd with two arguments", []() {
			AssertThat(Math::gcd(5,7), Equals(1));
			AssertThat(Math::gcd(5,15), Equals(5));
			AssertThat(Math::gcd(6,9), Equals(3));
		});
		it("computes gcd with array of arguments", []() {
			AssertThat(Math::gcd(Array<int>({5,7,11})), Equals(1));
			AssertThat(Math::gcd(Array<int>({6,12,45})), Equals(3));
		});
		testGcdAndLcm<int>("int");
		testGcdAndLcm<unsigned int>("unsigned int");
		testGcdAndLcm<long>("long");
		testGcdAndLcm<unsigned long>("unsigned long");
		testGcdAndLcm<long long>("long long");
		testGcdAndLcm<unsigned long long>("unsigned long long");

		testHarmonic();
	});
});
