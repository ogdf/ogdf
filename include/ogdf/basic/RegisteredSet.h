/** \file
 * \brief Declaration and implementation of ogdf::RegisteredSet.
 *
 * \author Niko Fink, Matthias Pfretzschner
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

template<class ElementType, class RegistryType, bool SupportFastSizeQuery = true>
class RegisteredSet {
public:
	using ListType =
			typename std::conditional<SupportFastSizeQuery, List<ElementType>, ListPure<ElementType>>::type;

	//! Creates an empty set associated with registry \p R.
	explicit RegisteredSet(const RegistryType& R) : m_it(&R) { }

	//! Creates an empty set associated with no registry.
	explicit RegisteredSet() : m_it() { }

	//! Reinitializes the set. Associates the set with no registry.
	void init() {
		m_it.init();
		m_elements.clear();
	}

	//! Reinitializes the set. Associates the set with graph \p R.
	void init(const RegistryType& R) {
		m_it.init(&R);
		m_elements.clear();
	}

	//! Inserts element \p v into this set.
	/**
	 * This operation has constant runtime.
	 * If the node is already contained in this set, nothing happens.
	 *
	 * \pre \p v is an element in the associated registry.
	 */
	void insert(ElementType v) {
		ListIterator<ElementType>& itV = m_it[v];
		if (!itV.valid()) {
			itV = m_elements.pushBack(v);
		}
	}

	//! Removes element \p v from this set.
	/**
	 * This operation has constant runtime.
	 *
	 * \pre \p v is an element in the associated registry.
	 * If the element is not contained in this set, nothing happens.
	 */
	void remove(ElementType v) {
		ListIterator<ElementType>& itV = m_it[v];
		if (itV.valid()) {
			m_elements.del(itV);
			itV = ListIterator<ElementType>();
		}
	}

	//! Removes all elements from this set.
	/**
	 * After this operation, this set is empty and still associated with the same registry.
	 * The runtime of this operations is linear in the #size().
	 */
	void clear() {
		m_it.fill(ListIterator<ElementType>());
		m_elements.clear();
	}

	//! Returns \c true iff element \p v is contained in this set.
	/**
	 * This operation has constant runtime.
	 *
	 * \pre \p v is an element in the associated graph.
	 */
	bool isMember(ElementType v) const { return m_it[v].valid(); }

	//! Returns a reference to the list of elements contained in this set.
	const ListType& elements() const { return m_elements; }

	//! Returns the associated registry.
	const RegistryType* registeredAt() const { return m_it.registeredAt(); }

	//! Returns the number of elements in this set.
	/**
	 * This operation has either linear or constant runtime, depending on \a SupportFastSizeQuery.
	 */
	int size() const { return m_elements.size(); }

	//! Copy constructor.
	template<bool OtherSupportsFastSizeQuery>
	explicit RegisteredSet(
			const RegisteredSet<ElementType, RegistryType, OtherSupportsFastSizeQuery>& other) {
		this = other;
	}

	//! Assignment operator.
	template<bool OtherSupportsFastSizeQuery>
	RegisteredSet& operator=(
			const RegisteredSet<ElementType, RegistryType, OtherSupportsFastSizeQuery>& other) {
		m_elements.clear();
		m_it.init(other.registeredAt());
		for (ElementType v : other.elements()) {
			insert(v);
		}
	}

private:
	//! #m_it[\a v] contains the list iterator pointing to \a v if \a v is contained in this set,
	//! or an invalid list iterator otherwise.
	RegisteredArrayWithoutDefault<RegistryType, ListIterator<ElementType>> m_it;

	//! The list of elements contained in this set.
	ListType m_elements;
};

}
