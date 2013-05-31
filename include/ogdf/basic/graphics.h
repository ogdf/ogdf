/*
 * $Revision: 3521 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-05-31 14:52:33 +0200 (Fr, 31. Mai 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration of basic types for graphics.
 *
 * \author Carsten Gutwenger
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

#ifdef _MSC_VER
#pragma once
#endif

#ifndef OGDF_GRAPHICS_H
#define OGDF_GRAPHICS_H

#include <ogdf/basic/basic.h>


namespace ogdf {

	//! Line types of strokes.
	enum StrokeType {
		stNone,			//!< no line
		stSolid,		//!< solid line
		stDash,			//!< dashed line
		stDot,			//!< dotted line
		stDashdot,		//!< line style "dash dot dash dot ..."
		stDashdotdot	//!< line style "dash dot dot dash dot dot ..."
	};

	//! Converts integer \a i to stroke type.
	StrokeType intToStrokeType(int i);


	//! Line cap types of strokes.
	enum StrokeLineCap {
		slcButt,
		slcRound,
		slcSquare
	};


	//! Line join types of strokes.
	enum StrokeLineJoin {
		sljMiter,
		sljRound,
		sljBevel
	};


	//! Fill patterns.
	enum FillPattern {
		fpNone,
		fpSolid,
		fpDense1,
		fpDense2,
		fpDense3,
		fpDense4,
		fpDense5,
		fpDense6,
		fpDense7,
		fpHorizontal,
		fpVertical,
		fpCross,
		fpBackwardDiagonal,
		fpForwardDiagonal,
		fpDiagonalCross
	};

	//! Converts integer \a i to fill pattern.
	FillPattern intToFillPattern(int i);


	//! Types for node shapes.
	enum Shape {
		shRect,               //!< rectangle
		shRoundedRect,        //!< rectangle with rounded corners
		shEllipse,            //!< ellipse
		shTriangle,           //!< isosceles triangle (base side down)
		shPentagon,           //!< pentagon
		shHexagon,            //!< hexagon
		shOctagon,            //!< octagon
		shRhomb,              //!< rhomb (=diamond)
		shTrapeze,            //!< trapeze (upper side shorter)
		shParallelogram,      //!< parallelogram (slanted to the right)
		shInvTriangle,        //!< isosceles triangle (base side up)
		shInvTrapeze,         //!< inverted trapeze  (upper side longer)
		shInvParallelogram,   //!< inverted parallelogram (slanted to the left)
		shImage
	};


	//! Types for edge arrows.
	enum EdgeArrow {
		eaNone,		//!< no edge arrows
		eaLast,		//!< edge arrow at target node of the edge
		eaFirst,	//!< edge arrow at source node of the edge
		eaBoth,		//!< edge arrow at target and source node of the edge
		eaUndefined
	};



	//! Colors reresented as RGBA values.
	/**
	 * The Color class represents colors with four components: R (red), G (green), B (blue), and A (alpha channel).
	 * Each component has a value between and 255. The alpha channel controls tranparency, where an opaque color
	 * has an alpha channel of 255.
	 */
	class Color {
		__uint8 m_red, m_green, m_blue, m_alpha;

	public:
		//! Named colors (same as SVG color keywords).
		enum Name {
			Aliceblue,
			Antiquewhite,
			Aqua,
			Aquamarine,
			Azure,
			Beige,
			Bisque,
			Black,
			Blanchedalmond,
			Blue,
			Blueviolet,
			Brown,
			Burlywood,
			Cadetblue,
			Chartreuse,
			Chocolate,
			Coral,
			Cornflowerblue,
			Cornsilk,
			Crimson,
			Cyan,
			Darkblue,
			Darkcyan,
			Darkgoldenrod,
			Darkgray,
			Darkgreen,
			Darkgrey,
			Darkkhaki,
			Darkmagenta,
			Darkolivegreen,
			Darkorange,
			Darkorchid,
			Darkred,
			Darksalmon,
			Darkseagreen,
			Darkslateblue,
			Darkslategray,
			Darkslategrey,
			Darkturquoise,
			Darkviolet,
			Deeppink,
			Deepskyblue,
			Dimgray,
			Dimgrey,
			Dodgerblue,
			Firebrick,
			Floralwhite,
			Forestgreen,
			Fuchsia,
			Gainsboro,
			Ghostwhite,
			Gold,
			Goldenrod,
			Gray,
			Green,
			Greenyellow,
			Grey,
			Honeydew,
			Hotpink,
			Indianred,
			Indigo,
			Ivory,
			Khaki,
			Lavender,
			Lavenderblush,
			Lawngreen,
			Lemonchiffon,
			Lightblue,
			Lightcoral,
			Lightcyan,
			Lightgoldenrodyellow,
			Lightgray,
			Lightgreen,
			Lightgrey,
			Lightpink,
			Lightsalmon,
			Lightseagreen,
			Lightskyblue,
			Lightslategray,
			Lightslategrey,
			Lightsteelblue,
			Lightyellow,
			Lime,
			Limegreen,
			Linen,
			Magenta,
			Maroon,
			Mediumaquamarine,
			Mediumblue,
			Mediumorchid,
			Mediumpurple,
			Mediumseagreen,
			Mediumslateblue,
			Mediumspringgreen,
			Mediumturquoise,
			Mediumvioletred,
			Midnightblue,
			Mintcream,
			Mistyrose,
			Moccasin,
			Navajowhite,
			Navy,
			Oldlace,
			Olive,
			Olivedrab,
			Orange,
			Orangered,
			Orchid,
			Palegoldenrod,
			Palegreen,
			Paleturquoise,
			Palevioletred,
			Papayawhip,
			Peachpuff,
			Peru,
			Pink,
			Plum,
			Powderblue,
			Purple,
			Red,
			Rosybrown,
			Royalblue,
			Saddlebrown,
			Salmon,
			Sandybrown,
			Seagreen,
			Seashell,
			Sienna,
			Silver,
			Skyblue,
			Slateblue,
			Slategray,
			Slategrey,
			Snow,
			Springgreen,
			Steelblue,
			Tan,
			Teal,
			Thistle,
			Tomato,
			Turquoise,
			Violet,
			Wheat,
			White,
			Whitesmoke,
			Yellow,
			Yellowgreen
		};

		//! Creates an opaque black color.
		Color() : m_red(0), m_green(0), m_blue(0), m_alpha(255) { }

		//! Creates a color from given RGBA-values.
		Color(__uint8 r, __uint8 g, __uint8 b, __uint8 a = 255) : m_red(r), m_green(g), m_blue(b), m_alpha(a) { }

		//! Creates a color from given RGBA-values.
		Color(int r, int g, int b, int a = 255) : m_red((__uint8)r), m_green((__uint8)g), m_blue((__uint8)b), m_alpha((__uint8)a) { }

		//! Creates a color from given color name \a name.
		Color(Color::Name name);

		//! Crates a color from string \a str.
		Color(const string &str) { fromString(str); }

		//! Crates a color from string \a str.
		Color(const char *str) { fromString(string(str)); }

		//! Returns the red component.
		__uint8 red() const { return m_red; }

		//! Returns the green component.
		__uint8 green() const { return m_green; }

		//! Returns the blue component.
		__uint8 blue() const { return m_blue; }

		//! Returns the alpha channel.
		__uint8 alpha() const { return m_alpha; }

		//! Sets the red component to \a r.
		void red(__uint8 r) { m_red = r; }

		//! Sets the green component to \a g.
		void green(__uint8 g) { m_green = g; }

		//! Sets the blue component to \a b.
		void blue(__uint8 b) { m_blue = b; }

		//! Sets the alpha channel to \a a.
		void alpha(__uint8 a) { m_alpha = a; }

		//! Converts the color to a string and returns it.
		/**
		 * Colors as represented as strings using the \#RRGGBB hex notation.
		 * Please note that in this notation the alpha channel is not represented and
		 * is assumed to be 255 (an opaque color).
		 */
		string toString() const;

		//! Sets the color the the color defined by \a str.
		bool fromString(const string &str);

		//! Returns true iff \a c and this color are equal in every component.
		bool operator==(const Color &c) const {
			return m_red == c.m_red && m_green == c.m_green && m_blue == c.m_blue && m_alpha == c.m_alpha;
		}

		//! Returns true iff \a c and this color differ in any component.
		bool operator!=(const Color &c) const {
			return !operator==(c);
		}

		//! Writes the string representation of color \a c to output stream \a os.
		friend ostream &operator<<(ostream &os, const Color &c) {
			return os << c.toString();
		}
	};


	//! Properties of strokes.
	struct Stroke {
		Color          m_color;    //!< stroke color
		float          m_width;    //!< stroke width
		StrokeType     m_type : 8; //!< stroke type (e.g. solid or dashed)
		StrokeLineCap  m_cap  : 8; //!< line-cap of the stroke
		StrokeLineJoin m_join : 8; //!< line-join of the stroke

		Stroke() : m_color(Color::Black), m_width(1.0f), m_type(stSolid), m_cap(slcButt), m_join(sljMiter) { }
		Stroke(Color c) : m_color(c), m_width(1.0f), m_type(stSolid), m_cap(slcButt), m_join(sljMiter) { }
	};


	//! Properties of fills.
	struct Fill {
		Color       m_color;   //!< fill color
		Color       m_bgColor; //!< background color of fill pattern
		FillPattern m_pattern; //!< fill pattern

		Fill() : m_color(Color::White), m_bgColor(Color::Black), m_pattern(fpSolid) { }
		Fill(Color c) : m_color(c), m_bgColor(Color::Black), m_pattern(fpSolid) { }
		Fill(Color c, FillPattern pattern) : m_color(c), m_bgColor(Color::Black), m_pattern(pattern) { }
		Fill(Color c, Color bgColor, FillPattern pattern) : m_color(c), m_bgColor(bgColor), m_pattern(pattern) { }
	};


} // end namespace ogdf


#endif
