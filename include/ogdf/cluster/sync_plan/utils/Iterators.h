#pragma once

#include <utility>

template<typename T1, typename IT1, typename T2, typename IT2>
class ZipIterator {
	IT1 m_iter1;
	IT1 m_end1;
	IT2 m_iter2;
	IT2 m_end2;

public:
	explicit ZipIterator() { }

	explicit ZipIterator(const IT1& mIter1, const IT1& mEnd1, const IT2& mIter2, const IT2& mEnd2)
		: m_iter1(mIter1), m_end1(mEnd1), m_iter2(mIter2), m_end2(mEnd2) { }

	explicit ZipIterator(const IT1& mIter1, const IT1& mEnd1, const IT2& mIter2, const IT2& mEnd2,
			T1* t1, T2* t2)
		: m_iter1(mIter1), m_end1(mEnd1), m_iter2(mIter2), m_end2(mEnd2) { }

	bool operator==(const ZipIterator<T1, IT1, T2, IT2>& rhs) const {
		return m_iter1 == rhs.m_iter1 && m_end1 == rhs.m_end1 && m_iter2 == rhs.m_iter2
				&& m_end2 == rhs.m_end2;
	}

	bool operator!=(const ZipIterator<T1, IT1, T2, IT2>& rhs) const { return !(rhs == *this); }

	ZipIterator<T1, IT1, T2, IT2> begin() const { return *this; }

	ZipIterator<T1, IT1, T2, IT2> end() const {
		return ZipIterator(m_end1, m_end1, m_end2, m_end2);
	}

	std::pair<T1, T2> operator*() { return std::make_pair(*m_iter1, *m_iter2); }

	//! Increment operator (prefix, returns result).
	ZipIterator<T1, IT1, T2, IT2>& operator++() {
		OGDF_ASSERT(m_iter1 != m_end1);
		OGDF_ASSERT(m_iter2 != m_end2);
		++m_iter1;
		++m_iter2;
		return *this;
	}

	//! Increment operator (postfix, returns previous value).
	ZipIterator<T1, IT1, T2, IT2> operator++(int) {
		OGDF_ASSERT(m_iter1 != m_end1);
		OGDF_ASSERT(m_iter2 != m_end2);
		ZipIterator<T1, IT1, T2, IT2> before = *this;
		++m_iter1;
		++m_iter2;
		return before;
	}

	operator bool() const { return valid(); }

	bool valid() const { return m_iter1 != m_end1 && m_iter2 != m_end2; }

	bool bothEnded() { return m_iter1 == m_end1 && m_iter2 == m_end2; }
};

template<typename Wrapped, typename Key>
class FilteringIterator {
	Wrapped m_iter;
	Wrapped m_end;
	std::function<bool(Wrapped)> m_filter;

public:
	// iterator traits
	using iterator_category = std::input_iterator_tag;
	using value_type = Wrapped::value_type;
	using difference_type = Wrapped::difference_type;
	using pointer = Wrapped::pointer;
	using reference = Wrapped::reference;

	explicit FilteringIterator() = default;

	explicit FilteringIterator(Wrapped mIter, Wrapped mEnd,
			const std::function<bool(Wrapped)>& mFilter)
		: m_iter(mIter), m_end(mEnd), m_filter(mFilter) {
		next(true);
	}

	bool operator==(const FilteringIterator& rhs) const { return m_iter == rhs.m_iter; }

	bool operator!=(const FilteringIterator& rhs) const { return m_iter != rhs.m_iter; }

	FilteringIterator<Wrapped, Key> begin() const { return *this; }

	FilteringIterator<Wrapped, Key> end() const {
		return FilteringIterator(m_end, m_end, m_filter);
	}

	Key operator*() { return *m_iter; }

	//! Increment operator (prefix, returns result).
	FilteringIterator<Wrapped, Key>& operator++() {
		next(false);
		return *this;
	}

	//! Increment operator (postfix, returns previous value).
	FilteringIterator<Wrapped, Key> operator++(int) {
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
using StdFilteringIterator = FilteringIterator<Wrapped, typename Wrapped::value_type>;
