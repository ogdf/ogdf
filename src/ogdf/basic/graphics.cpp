/*
 * $Revision: 3572 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-06-20 14:39:41 +0200 (Do, 20. Jun 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of class GraphAttributes.
 *
 * Class GraphAttributes extends a graph by graphical attributes like
 * node position, color, etc.
 *
 * \author Carsten Gutwenger, Karsten Klein
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

#include <ogdf/basic/graphics.h>


namespace ogdf {

	StrokeType intToStrokeType(int i)
	{
		switch (i) {
		case 0:
			return stNone;
		case 1:
			return stSolid;
		case 2:
			return stDash;
		case 3:
			return stDot;
		case 4:
			return stDashdot;
		case 5:
			return stDashdotdot;
		default:
			return stSolid;
		}
	}

	FillPattern intToFillPattern(int i)
	{
		switch (i) {
		case 0:
			return fpNone;
		case 1:
			return fpSolid;
		case 2:
			return fpDense1;
		case 3:
			return fpDense2;
		case 4:
			return fpDense3;
		case 5:
			return fpDense4;
		case 6:
			return fpDense5;
		case 7:
			return fpDense6;
		case 8:
			return fpDense7;
		case 9:
			return fpHorizontal;
		case 10:
			return fpVertical;
		case 11:
			return fpCross;
		case 12:
			return fpBackwardDiagonal;
		case 13:
			return fpForwardDiagonal;
		case 14:
			return fpDiagonalCross;
		default:
			return fpNone;
		}
	}


	__uint8 rgbOfColor[][3] = {
		{ 240, 248, 255 },  // Aliceblue
		{ 250, 235, 215 },  // Antiquewhite
		{   0, 255, 255 },  // Aqua
		{ 127, 255, 212 },  // Aquamarine
		{ 240, 255, 255 },  // Azure
		{ 245, 245, 220 },  // Beige
		{ 255, 228, 196 },  // Bisque
		{   0,   0,   0 },  // Black
		{ 255, 235, 205 },  // Blanchedalmond
		{   0,   0, 255 },  // Blue
		{ 138,  43, 226 },  // Blueviolet
		{ 165,  42,  42 },  // Brown
		{ 222, 184, 135 },  // Burlywood
		{  95, 158, 160 },  // Cadetblue
		{ 127, 255,   0 },  // Chartreuse
		{ 210, 105,  30 },  // Chocolate
		{ 255, 127,  80 },  // Coral
		{ 100, 149, 237 },  // Cornflowerblue
		{ 255, 248, 229 },  // Cornsilk
		{ 220,  20,  60 },  // Crimson
		{   0, 255, 255 },  // Cyan
		{   0,   0, 139 },  // Darkblue
		{   0, 139, 139 },  // Darkcyan
		{ 184, 134,  11 },  // Darkgoldenrod
		{ 169, 169, 169 },  // Darkgray
		{   0, 100,   0 },  // Darkgreen
		{ 169, 169, 169 },  // Darkgrey
		{ 189, 183, 107 },  // Darkkhaki
		{ 139,   0, 139 },  // Darkmagenta
		{  85, 107,  47 },  // Darkolivegreen
		{ 255, 140,   0 },  // Darkorange
		{ 153,  50, 204 },  // Darkorchid
		{ 139,   0,   0 },  // Darkred
		{ 233, 150, 122 },  // Darksalmon
		{ 143, 188, 143 },  // Darkseagreen
		{  72,  61, 139 },  // Darkslateblue
		{  47,  79,  79 },  // Darkslategray
		{  47,  79,  79 },  // Darkslategrey
		{   0, 206, 209 },  // Darkturquoise
		{ 148,   0, 211 },  // Darkviolet
		{ 255,  20, 147 },  // Deeppink
		{   0, 191, 255 },  // Deepskyblue
		{ 105, 105, 105 },  // Dimgray
		{ 105, 105, 105 },  // Dimgrey
		{  30, 144, 255 },  // Dodgerblue
		{ 178,  34,  34 },  // Firebrick
		{ 255, 250, 240 },  // Floralwhite
		{  34, 139,  34 },  // Forestgreen
		{ 255,   0, 255 },  // Fuchsia
		{ 220, 220, 220 },  // Gainsboro
		{ 248, 248, 255 },  // Ghostwhite
		{ 255, 215,   0 },  // Gold
		{ 218, 165,  32 },  // Goldenrod
		{ 128, 128, 128 },  // Gray
		{   0, 128,   0 },  // Green
		{ 173, 255,  47 },  // Greenyellow
		{ 128, 128, 128 },  // Grey
		{ 240, 255, 240 },  // Honeydew
		{ 255, 105, 180 },  // Hotpink
		{ 205,  92,  92 },  // Indianred
		{  75,   0, 130 },  // Indigo
		{ 255, 255, 240 },  // Ivory
		{ 240, 230, 140 },  // Khaki
		{ 230, 230, 250 },  // Lavender
		{ 255, 240, 245 },  // Lavenderblush
		{ 124, 252,   0 },  // Lawngreen
		{ 255, 250, 205 },  // Lemonchiffon
		{ 173, 216, 230 },  // Lightblue
		{ 240, 128, 128 },  // Lightcoral
		{ 224, 255, 255 },  // Lightcyan
		{ 250, 250, 210 },  // Lightgoldenrodyellow
		{ 211, 211, 211 },  // Lightgray
		{ 144, 238, 144 },  // Lightgreen
		{ 211, 211, 211 },  // Lightgrey
		{ 255, 182, 193 },  // Lightpink
		{ 255, 160, 122 },  // Lightsalmon
		{  32, 178, 170 },  // Lightseagreen
		{ 135, 206, 250 },  // Lightskyblue
		{ 119, 136, 153 },  // Lightslategray
		{ 119, 136, 153 },  // Lightslategrey
		{ 176, 196, 222 },  // Lightsteelblue
		{ 255, 255, 224 },  // Lightyellow
		{   0, 255,   0 },  // Lime
		{  50, 205,  50 },  // Limegreen
		{ 250, 240, 230 },  // Linen
		{ 255,   0, 255 },  // Magenta
		{ 128,   0,   0 },  // Maroon
		{ 102, 205, 170 },  // Mediumaquamarine
		{   0,   0, 205 },  // Mediumblue
		{ 186,  85, 211 },  // Mediumorchid
		{ 147, 112, 219 },  // Mediumpurple
		{  60, 179, 113 },  // Mediumseagreen
		{ 123, 104, 238 },  // Mediumslateblue
		{   0, 250, 154 },  // Mediumspringgreen
		{  72, 209, 204 },  // Mediumturquoise
		{ 199,  21, 133 },  // Mediumvioletred
		{  25,  25, 112 },  // Midnightblue
		{ 245, 255, 250 },  // Mintcream
		{ 255, 228, 225 },  // Mistyrose
		{ 255, 228, 181 },  // Moccasin
		{ 255, 222, 173 },  // Navajowhite
		{   0,   0, 128 },  // Navy
		{ 253, 245, 230 },  // Oldlace
		{ 128, 128,   0 },  // Olive
		{ 107, 142,  35 },  // Olivedrab
		{ 255, 165,   0 },  // Orange
		{ 255,  69,   0 },  // Orangered
		{ 218, 112, 214 },  // Orchid
		{ 238, 232, 170 },  // Palegoldenrod
		{ 152, 251, 152 },  // Palegreen
		{ 175, 238, 238 },  // Paleturquoise
		{ 219, 112, 147 },  // Palevioletred
		{ 255, 239, 213 },  // Papayawhip
		{ 255, 218, 185 },  // Peachpuff
		{ 205, 133,  63 },  // Peru
		{ 255, 192, 203 },  // Pink
		{ 221, 160, 221 },  // Plum
		{ 176, 224, 230 },  // Powderblue
		{ 128,   0, 128 },  // Purple
		{ 255,   0,   0 },  // Red
		{ 188, 143, 143 },  // Rosybrown
		{  65, 105, 225 },  // Royalblue
		{ 139,  69,  19 },  // Saddlebrown
		{ 250, 128, 114 },  // Salmon
		{ 244, 164,  96 },  // Sandybrown
		{  46, 139,  87 },  // Seagreen
		{ 255, 245, 238 },  // Seashell
		{ 160,  82,  45 },  // Sienna
		{ 192, 192, 192 },  // Silver
		{ 135, 206, 235 },  // Skyblue
		{ 106,  90, 205 },  // Slateblue
		{ 112, 128, 144 },  // Slategray
		{ 112, 128, 144 },  // Slategrey
		{ 255, 250, 250 },  // Snow
		{   0, 255, 127 },  // Springgreen
		{  70, 130, 180 },  // Steelblue
		{ 210, 180, 140 },  // Tan
		{   0, 128, 128 },  // Teal
		{ 216, 191, 216 },  // Thistle
		{ 255,  99,  71 },  // Tomato
		{  64, 224, 208 },  // Turquoise
		{ 238, 130, 238 },  // Violet
		{ 245, 222, 179 },  // Wheat
		{ 255, 255, 255 },  // White
		{ 245, 245, 245 },  // Whitesmoke
		{ 255, 255,   0 },  // Yellow
		{ 154, 205,  50 }   // Yellowgreen
	};


	Color::Color(Color::Name name) :
		m_red  (rgbOfColor[name][0]),
		m_green(rgbOfColor[name][1]),
		m_blue (rgbOfColor[name][2]),
		m_alpha(255) { }


	inline char tohex(int value)
	{
		value &= 0xf;
		return (char) ((value < 10) ? '0'+value : 'A'+value-10);
	}


	string Color::toString() const
	{
		char str[8];

		str[0] = '#';
		str[1] = tohex(m_red >> 4);
		str[2] = tohex(m_red);
		str[3] = tohex(m_green >> 4);
		str[4] = tohex(m_green);
		str[5] = tohex(m_blue >> 4);
		str[6] = tohex(m_blue);
		str[7] = 0;

		return string(str);
	}


	inline __uint8 fromHex(char c)
	{
		return (__uint8)((isdigit((int)c) ? (c - '0') : (tolower((int)c) - 'a' + 10)) & 0xf);
	}


	bool Color::fromString(const string &str)
	{
		if(str.length() != 4 && str.length() != 7)
			return false;

		if(str[0] != '#') return false;
		for(string::size_type i = 1; i < str.length(); ++i) {
			if(!isxdigit((int)str[i]))
				return false;
		}

		if(str.length() == 7) {
			m_red   = (__uint8)( (fromHex(str[1]) << 4) + fromHex(str[2]) );
			m_green = (__uint8)( (fromHex(str[3]) << 4) + fromHex(str[4]) );
			m_blue  = (__uint8)( (fromHex(str[5]) << 4) + fromHex(str[6]) );

		} else {
			__uint8 v = fromHex(str[1]);
			m_red = (__uint8)( (v << 4) + v );

			v = fromHex(str[2]);
			m_green = (__uint8)( (v << 4) + v );

			v = fromHex(str[3]);
			m_blue = (__uint8)( (v << 4) + v );
		}

		m_alpha = 255;

		return true;
	}


} // end namespace ogdf
