/** \file
 * \brief Length attribute used in EmbedderMinDepthMaxFace.
 * It contains two components (d, l) and a linear order is defined by:
 * (d, l) > (d', l') iff d > d' or (d = d' and l > l')
 *
 * \author Thorsten Kerkhof
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

#include <ogdf/basic/basic.h>

namespace ogdf {

class MDMFLengthAttribute
{
public:
	//constructors and destructor
	MDMFLengthAttribute() { d = l = 0; }
	MDMFLengthAttribute(const int& d, const int& l) : d(d), l(l) { }
	MDMFLengthAttribute(const int& d) : d(d), l(0) { }
	MDMFLengthAttribute(const MDMFLengthAttribute& x) : d(x.d), l(x.l) { }
	~MDMFLengthAttribute() { }

	MDMFLengthAttribute &operator=(const MDMFLengthAttribute& x);
	MDMFLengthAttribute &operator=(int x);
	bool operator==(const MDMFLengthAttribute& x);
	bool operator!=(const MDMFLengthAttribute& x);
	bool operator>(const MDMFLengthAttribute& x);
	bool operator<(const MDMFLengthAttribute& x);
	bool operator>=(const MDMFLengthAttribute& x);
	bool operator<=(const MDMFLengthAttribute& x);
	MDMFLengthAttribute operator+(const MDMFLengthAttribute& x);
	MDMFLengthAttribute operator-(const MDMFLengthAttribute& x);
	MDMFLengthAttribute operator+=(const MDMFLengthAttribute& x);
	MDMFLengthAttribute operator-=(const MDMFLengthAttribute& x);

public:
	//the two components:
	int d;
	int l;
};

bool operator==(const MDMFLengthAttribute& x, const MDMFLengthAttribute& y);
bool operator!=(const MDMFLengthAttribute& x, const MDMFLengthAttribute& y);
bool operator>(const MDMFLengthAttribute& x, const MDMFLengthAttribute& y);
bool operator<(const MDMFLengthAttribute& x, const MDMFLengthAttribute& y);
bool operator>=(const MDMFLengthAttribute& x, const MDMFLengthAttribute& y);
bool operator<=(const MDMFLengthAttribute& x, const MDMFLengthAttribute& y);
MDMFLengthAttribute operator+(const MDMFLengthAttribute& x, const MDMFLengthAttribute& y);
MDMFLengthAttribute operator-(const MDMFLengthAttribute& x, const MDMFLengthAttribute& y);
MDMFLengthAttribute operator+=(const MDMFLengthAttribute& x, const MDMFLengthAttribute& y);
MDMFLengthAttribute operator-=(const MDMFLengthAttribute& x, const MDMFLengthAttribute& y);
ostream& operator<<(ostream& s, const MDMFLengthAttribute& x);

} // end namespace ogdf
