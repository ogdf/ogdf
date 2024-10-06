/** \file
 * \brief An intrusive list for the leaves of a PCTree. TODO should be moved to a central location; merge with GraphList?
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

#include <ogdf/basic/basic.h>
#include <ogdf/basic/internal/config_autogen.h>

#include <cstddef>
#include <iterator>

namespace ogdf::pc_tree {
template<class T>
class IntrusiveList {
	T* m_first = nullptr;
	T* m_last = nullptr;
	size_t m_count = 0;

public:
	class iterator {
		T* m_node;

	public:
		// iterator traits
		using iterator_category = std::input_iterator_tag;
		using value_type = T*;
		using difference_type = std::ptrdiff_t;
		using pointer = const T*;
		using reference = T*;

		explicit iterator(T* node) : m_node(node) { }

		iterator& operator++() {
			m_node = m_node->m_next;
			return *this;
		}

		iterator operator++(int) {
			iterator other = *this;
			++(*this);
			return other;
		}

		bool operator==(iterator other) const { return m_node == other.m_node; }

		bool operator!=(iterator other) const { return m_node != other.m_node; }

		T* operator*() const { return m_node; }
	};

	class node {
		friend class IntrusiveList;

	private:
		T* m_next;
		T* m_prev;
	};

	iterator begin() const { return iterator(m_first); }

	iterator end() const { return iterator(nullptr); }

	void clear() {
		m_first = nullptr;
		m_last = nullptr;
		m_count = 0;
	}

	[[nodiscard]] bool empty() const { return m_count == 0; }

	[[nodiscard]] size_t size() const { return m_count; }

	T* front() const {
		OGDF_ASSERT(m_first != nullptr);
		return m_first;
	}

	T* back() const {
		OGDF_ASSERT(m_last != nullptr);
		return m_last;
	}

	void push_front(T* obj) {
		OGDF_ASSERT(obj != nullptr);
		check();

		if (m_first == nullptr) {
			obj->m_next = nullptr;
			m_last = obj;
		} else {
			obj->m_next = m_first;
			m_first->m_prev = obj;
		}

		obj->m_prev = nullptr;
		m_first = obj;
		m_count++;
		check();
	}

	void push_back(T* obj) {
		OGDF_ASSERT(obj != nullptr);
		check();

		if (m_last == nullptr) {
			obj->m_prev = nullptr;
			m_first = obj;
		} else {
			obj->m_prev = m_last;
			m_last->m_next = obj;
		}

		obj->m_next = nullptr;
		m_last = obj;
		m_count++;
		check();
	}

	void pop_front() { erase(front()); }

	void pop_back() { erase(back()); }

	void erase(T* obj) {
		OGDF_ASSERT(obj != nullptr);
		OGDF_ASSERT(m_count > 0);
		check();

		if (obj == m_first) {
			m_first = obj->m_next;
		}

		if (obj == m_last) {
			m_last = obj->m_prev;
		}

		if (obj->m_prev) {
			obj->m_prev->m_next = obj->m_next;
		}

		if (obj->m_next) {
			obj->m_next->m_prev = obj->m_prev;
		}

		m_count--;
		check();
	}

	void splice(iterator at, IntrusiveList<T>& other) {
		OGDF_ASSERT(at == begin() || at == end());
		check();
		other.check();

		if (at == end()) {
			if (m_last != nullptr) {
				m_last->m_next = other.m_first;
				other.m_first->m_prev = m_last;
				m_last = other.m_last;
			} else {
				m_first = other.m_first;
				m_last = other.m_last;
			}
		} else if (at == begin()) {
			if (m_first != nullptr) {
				m_first->m_prev = other.m_last;
				other.m_last->m_next = m_first;
				m_first = other.m_first;
			} else {
				m_first = other.m_first;
				m_last = other.m_last;
			}
		}

		m_count += other.m_count;
		other.m_count = 0;
		other.m_first = nullptr;
		other.m_last = nullptr;
		check();
		other.check();
	}

private:
	void check() {
#ifdef OGDF_DEBUG
		size_t counter = 0;
		for ([[maybe_unused]] T* _ : *this) {
			counter++;
		}
		OGDF_ASSERT(counter == m_count);
#endif
	}
};
}
