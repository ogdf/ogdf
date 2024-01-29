/** \file
 * \brief Declaration and implementation of ogdf::RegisteredSet.
 *
 * \author Simon D. Fink, Matthias Pfretzschner
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

#include <ogdf/basic/List.h>
#include <ogdf/basic/RegisteredArray.h>

namespace ogdf {

//! Constant-time set operations.
/**
 * Maintains a subset of indexed keys managed by a registry.
 *
 * Provides efficient operations for testing membership,
 * iteration, insertion, and deletion of elements, as well as clearing the set.
 *
 * @tparam Registry The class which manages the registered keys. Must provide the functions defined in
 * class RegistryBase.
 * @tparam SupportFastSizeQuery Whether this set supports querying its #size in
 * constant instead of linear time (in the size).
 *
 * \sa RegisteredArray, NodeSet
 */
template<class Registry, bool SupportFastSizeQuery = true>
class RegisteredSet {
public:
	using element_type = typename Registry::key_type;
	using list_type = typename std::conditional<SupportFastSizeQuery, List<element_type>,
			ListPure<element_type>>::type;

	//! Creates an empty set associated with registry \p R.
	explicit RegisteredSet(const Registry& R) : m_it(&R) { }

	//! Creates an empty set associated with no registry.
	explicit RegisteredSet() : m_it() { }

	//! Reinitializes the set. Associates the set with no registry.
	void init() {
		m_it.init();
		m_elements.clear();
	}

	//! Reinitializes the set. Associates the set with registry \p R.
	void init(const Registry& R) {
		m_it.init(&R);
		m_elements.clear();
	}

	//! Inserts element \p v into this set.
	/**
	 * This operation has constant runtime.
	 * If the element is already contained in this set, nothing happens.
	 *
	 * \pre \p v is an element in the associated registry.
	 */
	void insert(element_type v) {
		ListIterator<element_type>& itV = m_it[v];
		if (!itV.valid()) {
			itV = m_elements.pushBack(v);
		}
	}

	//! Removes element \p v from this set.
	/**
	 * This operation has constant runtime.
	 * If the element is not contained in this set, nothing happens.
	 *
	 * \pre \p v is an element in the associated registry.
	 */
	void remove(element_type v) {
		ListIterator<element_type>& itV = m_it[v];
		if (itV.valid()) {
			m_elements.del(itV);
			itV = ListIterator<element_type>();
		}
	}

	//! Removes all elements from this set.
	/**
	 * After this operation, this set is empty and still associated with the same registry.
	 * The runtime of this operation is linear in the #size(). Implementation Detail:
	 * If less than 10% of all elements are part of this set, they will be individually cleared.
	 * Otherwise, \c std::vector::assign(...) will be used to clear all values.
	 */
	void clear() {
		if (!registeredAt()) {
			return;
		}
		if (size() * 10 < registeredAt()->getArraySize()) {
			while (!m_elements.empty()) {
				remove(m_elements.front());
			}
		} else {
			m_it.fill(ListIterator<element_type>());
			m_elements.clear();
		}
	}

	//! Returns \c true iff element \p v is contained in this set.
	/**
	 * This operation has constant runtime.
	 *
	 * \pre \p v is an element in the associated registry.
	 */
	bool isMember(element_type v) const { return m_it[v].valid(); }

	//! Returns the same as isMember() to use an RegisteredSet instance as filter function.
	bool operator()(element_type v) const { return isMember(v); }

	//! Returns a reference to the list of elements contained in this set.
	const list_type& elements() const { return m_elements; }

	//! Returns the associated registry.
	const Registry* registeredAt() const { return m_it.registeredAt(); }

	//! Returns the number of elements in this set.
	/**
	 * This operation has either linear or constant runtime, depending on \a SupportFastSizeQuery.
	 */
	int size() const { return m_elements.size(); }

	typename list_type::const_iterator begin() const { return m_elements.begin(); }

	typename list_type::const_iterator end() const { return m_elements.end(); }

	//! Copy constructor.
	template<bool OtherSupportsFastSizeQuery>
	explicit RegisteredSet(const RegisteredSet<Registry, OtherSupportsFastSizeQuery>& other) {
		this = other;
	}

	//! Assignment operator.
	template<bool OtherSupportsFastSizeQuery>
	RegisteredSet& operator=(const RegisteredSet<Registry, OtherSupportsFastSizeQuery>& other) {
		m_elements.clear();
		m_it.init(other.registeredAt());
		for (element_type v : other.elements()) {
			insert(v);
		}
	}

private:
	//! #m_it[\a v] contains the list iterator pointing to \a v if \a v is contained in this set,
	//! or an invalid list iterator otherwise.
	RegisteredArrayWithoutDefault<Registry, ListIterator<element_type>> m_it;

	//! The list of elements contained in this set.
	list_type m_elements;
};

}
