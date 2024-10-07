/** \file
 * \brief Utilities for wrapping Iterators as long as we have no std::ranges. TODO should be moved to a central location and discarded once we have C++20.
 *
 * \author Simon D. Fink <ogdf@niko.fink.bayern>
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

#include <iosfwd>
#include <iterator>
#include <utility>

namespace ogdf {

//! Simple before-C++20 version for std::ranges::ref_view
template<typename IT>
class Range {
	IT m_begin;
	IT m_end;

public:
	Range() { }

	explicit Range(IT begin) : m_begin(begin) { }

	Range(IT begin, IT end) : m_begin(begin), m_end(end) { }

	IT begin() const { return m_begin; }

	IT end() const { return m_end; }
};

//! Simple before-C++20 version for std::ranges::zip_view
template<typename IT1, typename IT2>
class ZipIterator {
	IT1 m_iter1;
	IT2 m_iter2;

public:
	// iterator traits
	using iterator_category = std::input_iterator_tag;
	using value_type = typename std::pair<typename IT1::value_type, typename IT2::value_type>;
	using difference_type = std::ptrdiff_t;
	using pointer = value_type*;
	using reference = value_type&;

	explicit ZipIterator() { }

	explicit ZipIterator(const IT1& mIter1, const IT2& mIter2)
		: m_iter1(mIter1), m_iter2(mIter2) { }

	bool operator==(const ZipIterator<IT1, IT2>& rhs) const {
		return m_iter1 == rhs.m_iter1 && m_iter2 == rhs.m_iter2;
	}

	bool operator!=(const ZipIterator<IT1, IT2>& rhs) const { return !(rhs == *this); }

	value_type operator*() { return std::make_pair(*m_iter1, *m_iter2); }

	bool checkOnlyOneEnded(const IT1& mEnd1, const IT2& mEnd2) {
		bool e1 = m_iter1 == mEnd1;
		bool e2 = m_iter2 == mEnd2;
		return e1 ^ e2;
	}

	//! Increment operator (prefix, returns result).
	ZipIterator<IT1, IT2>& operator++() {
		++m_iter1;
		++m_iter2;
		return *this;
	}

	//! Increment operator (postfix, returns previous value).
	ZipIterator<IT1, IT2> operator++(int) {
		ZipIterator<IT1, IT2> before = *this;
		++m_iter1;
		++m_iter2;
		return before;
	}
};

}
