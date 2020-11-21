/** \file
 * \brief Declaration and implementation of bounded queue class.
 *
 * \author Niko Fink
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

#include <ogdf/basic/Array.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/Math.h>

#ifndef OGDF_MEMORY_POOL_NTS

#	include <mutex>

#endif

namespace ogdf {

template<class Key, class Registry>
class RegisteredArrayBase {
	ListIterator<RegisteredArrayBase<Key, Registry>*> m_registration;
	const Registry* m_pRegistry;

public:
	//! Initializes an adjacency entry array not associated with a graph.
	RegisteredArrayBase() : m_pRegistry(nullptr) { }

	//! Initializes an adjacency entry array associated with \p pG.
	explicit RegisteredArrayBase(const Registry& registry) : m_pRegistry(&registry) {
		m_registration = registry.registerArray(this);
	}

	//! Initializes an adjacency entry array associated with \p pG.
	explicit RegisteredArrayBase(const Registry* registry) : m_pRegistry(registry) {
		if (registry) {
			m_registration = registry->registerArray(this);
		}
	}

	//! Moves adjacency entry array \p move_from to this adjacency entry array.
	RegisteredArrayBase(RegisteredArrayBase<Key, Registry>&& move_from) noexcept
		: m_registration(move_from.m_registration), m_pRegistry(move_from.m_pRegistry) {
		m_pRegistry->moveRegisterArray(m_registration, this);
		move_from.m_pRegistry = nullptr;
		move_from.m_registration = ListIterator<RegisteredArrayBase<Key, Registry>*>();
	}

	//! Destructor, unregisters the array
	virtual ~RegisteredArrayBase() {
		if (m_pRegistry) {
			m_pRegistry->unregisterArray(m_registration);
			disconnect();
		}
	}

	// event interface used by Graph
	//! Virtual function called when table size has to be enlarged.
	virtual void enlargeTable(int newTableSize) = 0;

	//! Virtual function called when table has to be reinitialized.
	virtual void reinit(int initTableSize) = 0;

	//! Virtual function called when array is disconnected from the graph.
	virtual void disconnect() { m_pRegistry = nullptr; }

	//! Virtual function called when the index of an adjacency entry is changed.
	virtual void resetIndex(int newIndex, int oldIndex) = 0;

	//! Associates the array with a new graph.
	void reregister(const Registry* registry) {
		if (m_pRegistry) {
			m_pRegistry->unregisterArray(m_registration);
		}
		m_pRegistry = registry;
		if (m_pRegistry != nullptr) {
			m_registration = m_pRegistry->registerArray(this);
		} else {
			m_registration = nullptr;
		}
	}

	//! Moves array registration from \p move_from to this array.
	void moveRegister(RegisteredArrayBase<Key, Registry>& move_from) {
		if (m_pRegistry) {
			m_pRegistry->unregisterArray(m_registration);
		}
		m_pRegistry = move_from.m_pRegistry;
		m_registration = move_from.m_registration;
		move_from.m_pRegistry = nullptr;
		move_from.m_registration = ListIterator<RegisteredArrayBase<Key, Registry>*>();
		if (m_pRegistry != nullptr) {
			m_pRegistry->moveRegisterArray(m_registration, this);
		}
	}

	//! Returns a pointer to the associated graph.
	const Registry* registeredAt() const { return m_pRegistry; }
};

template<class ArrayType, bool isConst = false>
class RegisteredArrayIterator {
public:
	//! Index type of the associated array.
	using registry_type = typename ArrayType::registry_type;
	using iterator_type = typename ArrayType::registry_type::iterator_type;
	using key_type = typename ArrayType::key_type;
	//! Value type of the associated array.
	using value_type = typename std::conditional<isConst, const typename ArrayType::value_type,
			typename ArrayType::value_type>::type;
	//! Type of the array.
	using array_pointer_type = typename std::conditional<isConst, const ArrayType*, ArrayType*>::type;

private:
	iterator_type m_it;
	array_pointer_type m_array; //!< The array.

public:
	RegisteredArrayIterator(iterator_type mIt, array_pointer_type mArray)
		: m_it(mIt), m_array(mArray) { }

	//! Index in #m_array.
	key_type key() const { return *m_it; }

	//! Value of #m_array at index #m_key.
	value_type& value() const { return (*m_array)[*m_it]; }

	//! Value of #m_array at index #m_key.
	value_type& operator*() const { return (*m_array)[*m_it]; }

	//! Equality operator.
	bool operator==(const RegisteredArrayIterator<ArrayType, isConst>& iter) const {
		return m_it == iter.m_it && m_array == iter.m_array;
	}

	//! Inequality operator.
	bool operator!=(const RegisteredArrayIterator<ArrayType, isConst>& iter) const {
		return !operator==(iter);
	}

	//! Increment operator (prefix).
	RegisteredArrayIterator<ArrayType, isConst>& operator++() {
		++m_it;
		return *this;
	}

	//! Increment operator (postfix).
	RegisteredArrayIterator<ArrayType, isConst> operator++(int) {
		RegisteredArrayIterator<ArrayType, isConst> iter = *this;
		++m_it;
		return iter;
	}

	//! Decrement operator (prefix).
	RegisteredArrayIterator<ArrayType, isConst>& operator--() {
		--m_it;
		return *this;
	}

	//! Decrement operator (postfix).
	RegisteredArrayIterator<ArrayType, isConst> operator--(int) {
		RegisteredArrayIterator<ArrayType, isConst> iter = *this;
		--m_it;
		return iter;
	}
};

template<class Registry, class Key, class Value>
class RegisteredArray : protected Array<Value>, protected RegisteredArrayBase<Key, Registry> {
	Value m_default_value; //!< The default value for array elements.

public:
	//! The type of the registry.
	using registry_type = Registry;
	//! The type for array keys.
	using key_type = Key;
	//! The type for array entries.
	using value_type = Value;

	using array = Array<Value>;
	using registered_array = RegisteredArray<Registry, Key, Value>;
	using registered_array_base = RegisteredArrayBase<Key, Registry>;

	//! The type for edge array iterators.
	using iterator = RegisteredArrayIterator<registered_array, false>;
	//! The type for edge array const iterators.
	using const_iterator = RegisteredArrayIterator<registered_array, true>;

	//! Constructs an empty adjacency entry array associated with no graph.
	RegisteredArray() : array(), registered_array_base() { }

	//! Constructs an adjacency entry array associated with \p G.
	explicit RegisteredArray(const Registry& pRegistry)
		: array(pRegistry.keyArrayTableSize()), registered_array_base(pRegistry) { }

	//! Constructs an adjacency entry array associated with \p G.
	/**
		 * @param R is the associated graph.
		 * @param default_value is the default value for all array elements.
		 */
	RegisteredArray(const Registry& pRegistry, const Value& default_value)
		: array(0, pRegistry.keyArrayTableSize() - 1, default_value)
		, registered_array_base(pRegistry)
		, m_default_value(default_value) { }

	//! Constructs an adjacency entry array that is a copy of \p A.
	/**
		 * Associates the array with the same graph as \p A and copies all elements.
		 */
	RegisteredArray(const registered_array& A)
		: array(A), registered_array_base(A.registeredAt()), m_default_value(A.m_default_value) { }

	//! Constructs an adjacency entry array containing the elements of \p A (move semantics).
	/**
		 * Adjacency entry array \p A is empty afterwards and not associated with any graph.
		 */
	RegisteredArray(registered_array&& A) noexcept
		: array(std::move(A))
		, registered_array_base(A.registeredAt())
		, m_default_value(A.m_default_value) { }

	/**
			  * @name Initialization and assignment
			  * These methods can be used to reinitialize the array, or to initialize all elements with a given value.
			  */
	//@{

	//! Reinitializes the array. Associates the array with no graph.
	void init() {
		array::init();
		registered_array::reregister(nullptr);
	}

	//! Reinitializes the array. Associates the array with \p G.
	void init(const Registry& pRegistry) {
		array::init(pRegistry.keyArrayTableSize());
		registered_array::reregister(&pRegistry);
	}

	//! Reinitializes the array. Associates the array with \p G.
	/**
		 * @param R is the associated graph.
		 * @param default_value is the default value.
		 */
	void init(const Registry& pRegistry, const Value& default_value) {
		array::init(0, pRegistry.keyArrayTableSize() - 1, m_default_value = default_value);
		registered_array::reregister(&pRegistry);
	}

	//! Sets all array elements to \p x.
	void fill(const Value& x) {
		if (!registeredAt()) {
			return;
		}
		int high = registeredAt()->maxKeyIndex();
		if (high >= 0) {
			array::fill(0, high, x);
		}
	}

	//! Assignment operator.
	registered_array& operator=(const registered_array& A) {
		array::operator=(A);
		m_default_value = A.m_default_value;
		registered_array_base::reregister(A.registeredAt());
		return *this;
	}

	//! Assignment operator (move semantics).
	/**
		 * Adjacency entry array \p a is empty afterwards and not associated with any graph.
		 */
	registered_array& operator=(registered_array&& a) {
		array::operator=(std::move(a));
		m_default_value = a.m_default_value;
		//			moveRegister(a.registeredAt());//TODO
		registered_array_base::reregister(a.registeredAt());
		a.init();
		return *this;
	}

	//@}
	/**
		 * @name Access methods
		 * These methods provide access to elements, size, and corresponding graph.
		 */
	//@{

	//! Returns true iff the array is associated with a graph.
	bool valid() const { return array::low() <= array::high(); }

	//! Returns a reference to the element with index \p adj.
	const Value& operator[](Key key) const {
		OGDF_ASSERT(registeredAt());
		OGDF_ASSERT(registeredAt()->isKeyAssociated(key));
		return array::operator[](registeredAt()->keyToIndex(key));
	}

	//! Returns a reference to the element with index \p adj.
	Value& operator[](Key key) {
		OGDF_ASSERT(registeredAt());
		OGDF_ASSERT(registeredAt()->isKeyAssociated(key));
		return array::operator[](registeredAt()->keyToIndex(key));
	}

	//! Returns a reference to the element with index \p adj.
	const Value& operator()(Key key) const {
		OGDF_ASSERT(registeredAt());
		OGDF_ASSERT(registeredAt()->isKeyAssociated(key));
		return array::operator[](registeredAt()->keyToIndex(key));
	}

	//! Returns a reference to the element with index \p adj.
	Value& operator()(Key key) {
		OGDF_ASSERT(registeredAt());
		OGDF_ASSERT(registeredAt()->isKeyAssociated(key));
		return array::operator[](registeredAt()->keyToIndex(key));
	}

	//@}
	/**
		 * @name Iterators
		 * These methods return bidirectional iterators to elements in the array.
		 */
	//@{

	//! Returns an iterator to the first entry in the edge array.
	/**
		 * If the edge array is empty, a null pointer iterator is returned.
		 */
	iterator begin() { return iterator(registeredAt()->begin(), this); }

	//! Returns a const iterator to the first entry in the edge array.
	/**
		 * If the edge array is empty, a null pointer iterator is returned.
		 */
	const_iterator begin() const { return const_iterator(registeredAt()->begin(), this); }

	//! Returns a const iterator to the first entry in the edge array.
	/**
		 * If the edge array is empty, a null pointer iterator is returned.
		 */
	const_iterator cbegin() const { return const_iterator(registeredAt()->begin(), this); }

	//! Returns an iterator to one-past-last entry in the edge array.
	/**
		 * This is always a null pointer iterator.
		 */
	iterator end() { return iterator(registeredAt()->end(), this); }

	//! Returns a const iterator to one-past-last entry in the edge array.
	/**
		 * This is always a null pointer iterator.
		 */
	const_iterator end() const { return const_iterator(registeredAt()->end(), this); }

	//! Returns a const iterator to one-past-last entry in the edge array.
	/**
		 * This is always a null pointer iterator.
		 */
	const_iterator cend() const { return const_iterator(registeredAt()->end(), this); }

	//@}

	using registered_array_base::registeredAt;

private:
	void enlargeTable(int newTableSize) override {
		array::grow(newTableSize - array::size(), m_default_value);
	}

	void reinit(int initTableSize) override { array::init(0, initTableSize - 1, m_default_value); }

	void resetIndex(int newIndex, int oldIndex) override {
		array::operator[](newIndex) = array::operator[](oldIndex);
	}

	void disconnect() override {
		array::init();
		registered_array_base::disconnect();
	}

	OGDF_NEW_DELETE
};

static constexpr int MIN_TABLE_SIZE = (1 << 4);

int calculateTableSize(int actualCount);

template<class Key, class Registry, typename Iterator = void>
class RegistryBase {
public:
	using registered_array_type = RegisteredArrayBase<Key, Registry>;
	using key_type = Key;
	using registry_type = Registry;
	using iterator_type = Iterator;

private:
	mutable ListPure<registered_array_type*> m_registeredArrays; //!< The registered arrays.

#ifndef OGDF_MEMORY_POOL_NTS
	mutable std::mutex m_mutexRegArrays; //!< The critical section for protecting shared acces to register/unregister methods.
#endif
public:
	virtual ~RegistryBase() {
		while (!m_registeredArrays.empty()) {
			m_registeredArrays.popFrontRet()->disconnect();
		}
	}

	ListIterator<registered_array_type*> registerArray(registered_array_type* pArray) const {
#ifndef OGDF_MEMORY_POOL_NTS
		std::lock_guard<std::mutex> guard(m_mutexRegArrays);
#endif
		return m_registeredArrays.pushBack(pArray);
	}

	void unregisterArray(ListIterator<registered_array_type*> it) const {
#ifndef OGDF_MEMORY_POOL_NTS
		std::lock_guard<std::mutex> guard(m_mutexRegArrays);
#endif
		m_registeredArrays.del(it);
	}

	void moveRegisterArray(ListIterator<registered_array_type*> it,
			registered_array_type* pArray) const {
#ifndef OGDF_MEMORY_POOL_NTS
		std::lock_guard<std::mutex> guard(m_mutexRegArrays);
#endif
		*it = pArray;
	}

	/**
		 * Returns whether a key is associated with this registry.
		 * @param key The key.
		 * @return Whether \p key is associated with this registry.
		 */
	[[nodiscard]] virtual bool isKeyAssociated(Key key) const = 0;

	/**
		 * Returns the associated index of a key.
		 * @param key The key.
		 * @return The index of the key.
		 */
	[[nodiscard]] virtual int keyToIndex(Key key) const = 0;

	[[nodiscard]] virtual int keyArrayTableSize() const = 0;

	[[nodiscard]] virtual int maxKeyIndex() const = 0;

	[[nodiscard]] virtual Iterator begin() const { }

	[[nodiscard]] virtual Iterator end() const { }

	void reinitArrays() {
		for (registered_array_type* ab : m_registeredArrays) {
			ab->reinit(this->keyArrayTableSize());
		}
	}

	void enlargeArrayTables() {
		for (registered_array_type* ab : m_registeredArrays) {
			ab->enlargeTable(this->keyArrayTableSize());
		}
	}

	void resetArrayEntryIndex(int newIndex, int oldIndex) {
		for (registered_array_type* ab : m_registeredArrays) {
			ab->resetIndex(newIndex, oldIndex);
		}
	}
};
}
