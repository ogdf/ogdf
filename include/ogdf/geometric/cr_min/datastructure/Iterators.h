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

#include <functional>
#include <utility>

namespace ogdf {
namespace internal {
namespace gcm {

namespace datastructure {
template<typename Iterator>
class FilterIterator {
protected:
	Iterator m_begin;
	Iterator m_cur;
	Iterator m_end;
	const std::function<bool(const typename Iterator::T)> m_predicate;

public:
	using T = typename Iterator::T;

	FilterIterator(Iterator begin, Iterator end,
			const std::function<bool(const typename Iterator::T)>& _predicate)
		: m_begin(std::move(begin))
		, m_cur(m_begin)
		, m_end(std::move(end))
		, m_predicate(std::move(_predicate)) {
		if (m_cur != m_end && !m_predicate(*m_cur)) {
			++(*this);
		}
		// nothing to do
	}

	FilterIterator<Iterator>& operator--() {
		do {
			--m_cur;
		} while (m_cur != m_begin && !m_predicate(*m_cur));
		return *this;
	}

	FilterIterator<Iterator>& operator++() {
		do {
			++m_cur;
		} while (m_cur != m_end && !m_predicate(*m_cur));
		return *this;
	}

	bool operator==(const FilterIterator<Iterator>& b) const { return m_cur == b.m_cur; }

	bool operator!=(const FilterIterator<Iterator>& b) const { return m_cur != b.m_cur; }

	T operator*() {
		if (m_cur != m_end && !m_predicate(*m_cur)) {
			++(*this);
		}
		return *m_cur;
	}
};

template<typename Iterator>
class IteratorRange {
private:
	Iterator m_begin;
	Iterator m_end;

public:
	IteratorRange(Iterator begin, Iterator end) : m_begin(begin), m_end(end) { }

	Iterator begin() { return m_begin; }

	Iterator end() { return m_end; }
};


}

}
}
}
