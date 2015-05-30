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
 * See README.txt in the root directory of the OGDF installation for details.
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
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/

#include <bandit/bandit.h>
#include <ogdf/basic/Math.h>

using namespace bandit;
using namespace ogdf;

go_bandit([]() {
	describe("Math.h", [&]() {
		it("computes gcd with two arguments", [&]() {
            AssertThat(Math::gcd(5,7), Equals(1));
            AssertThat(Math::gcd(5,15), Equals(5));
            AssertThat(Math::gcd(6,9), Equals(3));
		});
		it("computes gcd with array of arguments", [&]() {
            AssertThat(Math::gcd(Array<int>({5,7,11})), Equals(1));
            AssertThat(Math::gcd(Array<int>({6,12,45})), Equals(3));
		});
	});
});
