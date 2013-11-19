/*
 * $Revision: 3837 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 15:19:30 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration of useful methods for processing various fileformats.
 *
 * \author Łukasz Hanuszczak
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

#ifndef OGDF_FILEFORMAT_UTILS_H
#define OGDF_FILEFORMAT_UTILS_H

#include <ogdf/basic/Hashing.h>

#include <string>

namespace ogdf {

// Provides a nicer syntax for reading formatted input through streams, e.g.
// `stream >> a >> ';' >> y`.
class TokenIgnorer {
private:
	char m_c;

public:
	TokenIgnorer(const char c): m_c(c) {};

	friend std::istream &operator >>(std::istream &is, TokenIgnorer c);
};


std::istream &operator >>(std::istream &is, TokenIgnorer token);

template <typename E>
static inline E toEnum(
	const std::string &str, // A string we want to convert.
	Hashing<std::string, E> *&map, // A map to be lazily evaluated.
	std::string toString(const E&),
	const E first, const E last, const E def) // Enum informations.
{
	if(!map) {
		map = new Hashing<std::string, E>();

		// Iterating over enums is potentially unsafe... (fixable in C++11).
		for(int it = last; it >= first; it--) {
			const E e = static_cast<E>(it);
			map->insert(toString(e), e);
		}
	}

	HashElement<std::string, E> *elem = map->lookup(str);
	return elem ? elem->info() : def;
}


} // end namespace ogdf


#endif
