/** \file
 * \brief TODO Document
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

#include <functional>
#include <iosfwd>
#include <iterator>
#include <utility>

#pragma GCC diagnostic ignored "-Wshadow" // TODO remove

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

template<typename Wrapped>
class FilteringIterator {
	Wrapped m_iter;
	Wrapped m_end;
	std::function<bool(Wrapped)> m_filter;

public:
	// iterator traits
	using iterator_category = std::input_iterator_tag;
	using value_type = typename Wrapped::value_type;
	using difference_type = typename Wrapped::difference_type;
	using pointer = typename Wrapped::pointer;
	using reference = typename Wrapped::reference;

	explicit FilteringIterator() = default;

	explicit FilteringIterator(Wrapped mIter, Wrapped mEnd,
			const std::function<bool(Wrapped)>& mFilter)
		: m_iter(mIter), m_end(mEnd), m_filter(mFilter) {
		next(true);
	}

	bool operator==(const FilteringIterator& rhs) const { return m_iter == rhs.m_iter; }

	bool operator!=(const FilteringIterator& rhs) const { return m_iter != rhs.m_iter; }

	FilteringIterator<Wrapped> begin() const { return *this; }

	FilteringIterator<Wrapped> end() const { return FilteringIterator(m_end, m_end, m_filter); }

	value_type operator*() { return *m_iter; }

	//! Increment operator (prefix, returns result).
	FilteringIterator<Wrapped>& operator++() {
		next(false);
		return *this;
	}

	//! Increment operator (postfix, returns previous value).
	FilteringIterator<Wrapped> operator++(int) {
		FilteringIterator before = *this;
		next(false);
		return before;
	}

	void next(bool before_first) {
		if (!before_first) {
			OGDF_ASSERT(m_iter != m_end);
			++m_iter;
		}
		while (m_iter != m_end) {
			if (m_filter(m_iter)) {
				break;
			} else {
				++m_iter;
			}
		}
	}

	operator bool() const { return valid(); }

	bool valid() const { return m_iter != m_end; }
};

template<typename Wrapped>
class BoundsCheckingIterator {
	Wrapped m_iter;
	Wrapped m_end;

public:
	// iterator traits
	using iterator_category = std::input_iterator_tag;
	using value_type = typename Wrapped::value_type;
	using difference_type = typename Wrapped::difference_type;
	using pointer = typename Wrapped::pointer;
	using reference = typename Wrapped::reference;

	explicit BoundsCheckingIterator() = default;

	explicit BoundsCheckingIterator(Wrapped mIter, Wrapped mEnd) : m_iter(mIter), m_end(mEnd) { }

	bool operator==(const BoundsCheckingIterator& rhs) const { return m_iter == rhs.m_iter; }

	bool operator!=(const BoundsCheckingIterator& rhs) const { return m_iter != rhs.m_iter; }

	FilteringIterator<Wrapped> begin() const { return *this; }

	FilteringIterator<Wrapped> end() const { return FilteringIterator(m_end, m_end); }

	value_type operator*() { return *m_iter; }

	//! Increment operator (prefix, returns result).
	FilteringIterator<Wrapped>& operator++() {
		OGDF_ASSERT(m_iter != m_end);
		++m_iter;
		return *this;
	}

	//! Increment operator (postfix, returns previous value).
	FilteringIterator<Wrapped> operator++(int) {
		OGDF_ASSERT(m_iter != m_end);
		FilteringIterator before = *this;
		++m_iter;
		return before;
	}

	operator bool() const { return valid(); }

	bool valid() const { return m_iter != m_end; }
};

// TODO add TransformingIterator
