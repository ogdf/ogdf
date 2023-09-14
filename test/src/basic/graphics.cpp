/** \file
 * \brief Tests for ogdf/basic/graphics
 *
 * \author Dominik Potulski.
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


#include <ogdf/basic/graphics.h>

#include <testing.h>

template<class T>
void testEnumStringConversion(T t) {
	T conversion = fromString<T>(toString(t));
	it("converts " + toString(t), [&]() { AssertThat(conversion, Equals(t)); });
}

void testStrokeTypeStringConversion() {
	describe("Conversion of StrokeType to string and vice versa", []() {
		testEnumStringConversion(StrokeType::None);
		testEnumStringConversion(StrokeType::Solid);
		testEnumStringConversion(StrokeType::Dash);
		testEnumStringConversion(StrokeType::Dot);
		testEnumStringConversion(StrokeType::Dashdot);
		testEnumStringConversion(StrokeType::Dashdotdot);
	});
}

void testFillPatternStringConversion() {
	describe("Conversion of FillPattern to string and vice versa", []() {
		testEnumStringConversion(FillPattern::None);
		testEnumStringConversion(FillPattern::Solid);
		testEnumStringConversion(FillPattern::Dense1);
		testEnumStringConversion(FillPattern::Dense2);
		testEnumStringConversion(FillPattern::Dense3);
		testEnumStringConversion(FillPattern::Dense4);
		testEnumStringConversion(FillPattern::Dense5);
		testEnumStringConversion(FillPattern::Dense6);
		testEnumStringConversion(FillPattern::Dense7);
		testEnumStringConversion(FillPattern::Horizontal);
		testEnumStringConversion(FillPattern::Vertical);
		testEnumStringConversion(FillPattern::Cross);
		testEnumStringConversion(FillPattern::BackwardDiagonal);
		testEnumStringConversion(FillPattern::ForwardDiagonal);
		testEnumStringConversion(FillPattern::DiagonalCross);
	});
}

void testShapeStringConversion() {
	describe("Conversion of Shape to string and vice versa", []() {
		testEnumStringConversion(Shape::Rect);
		testEnumStringConversion(Shape::RoundedRect);
		testEnumStringConversion(Shape::Ellipse);
		testEnumStringConversion(Shape::Triangle);
		testEnumStringConversion(Shape::Pentagon);
		testEnumStringConversion(Shape::Hexagon);
		testEnumStringConversion(Shape::Octagon);
		testEnumStringConversion(Shape::Rhomb);
		testEnumStringConversion(Shape::Trapeze);
		testEnumStringConversion(Shape::Parallelogram);
		testEnumStringConversion(Shape::InvTriangle);
		testEnumStringConversion(Shape::InvTrapeze);
		testEnumStringConversion(Shape::InvParallelogram);
		testEnumStringConversion(Shape::Image);
	});
}

go_bandit([]() {
	testStrokeTypeStringConversion();
	testFillPatternStringConversion();
	testShapeStringConversion();
});
