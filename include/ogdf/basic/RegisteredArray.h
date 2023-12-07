/** \file
 * \brief Declaration and implementation of RegisteredArray class.
 *
 * \author Simon D. Fink
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

#include <ogdf/basic/Math.h>

#include <list>
#include <memory>

#ifndef OGDF_MEMORY_POOL_NTS

#	include <mutex>

#endif

namespace ogdf {
template<typename Registry>
class RegisteredArrayBase;

//! The default minimum table size for registered arrays.
static constexpr int MIN_TABLE_SIZE = (1 << 4);

//! The default growth function for registered arrays.
/**
 * @return The smallest power of 2 that is no less than \p actualCount and #MIN_TABLE_SIZE.
 */
inline int calculateTableSize(int actualCount) {
	return Math::nextPower2(MIN_TABLE_SIZE, actualCount);
}

//! Abstract base class for registries.
/**
 * Defines the interface for event handling regarding the indexed keys. A registry manages one key type and stores all
 * registered arrays associated with that key. It determines the new size of all registered arrays when
 * keys are added or removed.
 *
 * The following methods must be implemented by all subclasses as they are used via the
 * <a href="https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern">CRTP</a>:
 * \code{.cpp}
 * //! Returns the index of \p key.
 * static inline int keyToIndex(Key key);
 *
 * //! Returns whether \p key is associated with this registry.
 * bool isKeyAssociated(Key key) const;
 *
 * //! Returns the maximum index of all keys managed by this registry.
 * int maxKeyIndex() const;
 *
 * //! Returns the array size currently requested by this registry.
 * int calculateArraySize(int add) const;
 * \endcode
 *
 * \remark To avoid frequent costly resize operations, the array size returned by calculateArraySize
 *         should grow in larger steps (e.g. powers of 2)
 *
 * @tparam Key The key type the registry manages.
 * @tparam Registry The class that implements the interface defined in RegistryBase.
 * @tparam Iterator An iterator for all managed keys. Can be \c void if iterating the keys through the registered array
 * is not required. To allow iterating over all keys, define \c begin() and \c end() methods.
 *
 * \sa RegisteredArray
 */
template<typename Key, typename Registry, typename Iterator = void>
class RegistryBase {
public:
	using registered_array_type = RegisteredArrayBase<Registry>;
	using key_type = Key;
	using registry_type = Registry;
	using iterator_type = Iterator;
	using registration_list_type =
			std::list<registered_array_type*, OGDFAllocator<registered_array_type*>>;
	using registration_iterator_type = typename registration_list_type::iterator;

private:
	mutable registration_list_type m_registeredArrays;
	bool m_autoShrink = false;
	int m_size = 0;

#ifndef OGDF_MEMORY_POOL_NTS
	mutable std::mutex m_mutexRegArrays;
#endif

protected:
	RegistryBase() = default;

public:
	//! Destructor. Unregisters all associated arrays.
	virtual ~RegistryBase() noexcept { unregisterArrays(); }

	RegistryBase(const RegistryBase& copy) = delete;

	RegistryBase(RegistryBase&& move) noexcept = delete;

	RegistryBase& operator=(const RegistryBase& other) = delete;

	RegistryBase& operator=(RegistryBase&& other) noexcept = delete;

	//! Registers a new array with this registry.
	/**
	 * @param pArray A pointer to the registered array.
	 * @return An iterator pointing to the entry for the registered array in the list of registered arrays.
	 *         This iterator is required for unregistering the array again.
	 */
	OGDF_NODISCARD registration_iterator_type registerArray(registered_array_type* pArray) const {
#ifndef OGDF_MEMORY_POOL_NTS
		std::lock_guard<std::mutex> guard(m_mutexRegArrays);
#endif
		return m_registeredArrays.emplace(m_registeredArrays.end(), pArray);
	}

	//! Unregisters an array associated with this registry.
	/**
	 * @param it An iterator pointing to the entry of the array in the list of all registered arrays.
	 */
	void unregisterArray(registration_iterator_type it) const noexcept {
#ifndef OGDF_MEMORY_POOL_NTS
		std::lock_guard<std::mutex> guard(m_mutexRegArrays);
#endif
		m_registeredArrays.erase(it);
	}

	//! Stores array \p pArray at position \p it in the list of registered arrays.
	void moveRegisterArray(registration_iterator_type it, registered_array_type* pArray) const {
#ifndef OGDF_MEMORY_POOL_NTS
		std::lock_guard<std::mutex> guard(m_mutexRegArrays);
#endif
		*it = pArray;
	}

	//! Records the addition of a new key and resizes all registered arrays if necessary.
	void keyAdded(Key key) {
		if (static_cast<Registry*>(this)->keyToIndex(key) >= m_size) {
			resizeArrays();
		}
	}

	//! Records the deletion of a key and resizes all registered arrays if auto shrink is enabled.
	void keyRemoved(Key key) {
		if (m_autoShrink) {
			resizeArrays();
		}
	}

	//! Records that all keys have been cleared. If auto shrink is enabled, all arrays are cleared and resized to 0.
	void keysCleared() { resizeArrays(0); }

	//! Resizes all arrays to the size requested by calculateArraySize(). Only shrinks the arrays if auto shrink is
	//! enabled
	void resizeArrays() { resizeArrays(static_cast<Registry*>(this)->calculateArraySize(0)); }

	//! Resizes all arrays to \p size. Only shrinks the arrays if auto shrink is enabled
	void resizeArrays(int size) { resizeArrays(size, m_autoShrink); }

	//! Resizes all arrays to \p size. If \p shrink is \c true, the arrays may also shrink.
	void resizeArrays(int size, bool shrink) {
		if (size == m_size) {
			return;
		}
		m_size = size = max(size, 0);
		for (registered_array_type* ab : m_registeredArrays) {
			ab->resize(size, shrink);
		}
	}

	//! Resizes all arrays to make space of \p new_keys new keys.
	void reserveSpace(int new_keys) {
		resizeArrays(static_cast<Registry*>(this)->calculateArraySize(new_keys));
	}

	//! Swaps the entries at \p index1 and \p index2 in all registered arrays.
	void swapArrayEntries(int index1, int index2) {
		for (registered_array_type* ab : m_registeredArrays) {
			ab->swapEntries(index1, index2);
		}
	}

	//! Copies the entry from \p fromIndex to \p toIndex in all registered arrays.
	void copyArrayEntries(int toIndex, int fromIndex) {
		for (registered_array_type* ab : m_registeredArrays) {
			ab->copyEntry(toIndex, fromIndex);
		}
	}

	//! Unregister all associated arrays.
	void unregisterArrays() noexcept {
		while (!m_registeredArrays.empty()) {
#ifdef OGDF_DEBUG
			auto size = m_registeredArrays.size();
#endif
			m_registeredArrays.front()->unregister();
			OGDF_ASSERT(m_registeredArrays.size() < size);
		}
	}

	//! Returns a reference to the list of all registered arrays.
	const registration_list_type& getRegisteredArrays() const { return m_registeredArrays; }

	//! Returns whether the registry allows arrays to shrink when keys are removed.
	bool isAutoShrink() const { return m_autoShrink; }

	//! Specifies whether the registry allows arrays to shrink when keys are removed.
	void setAutoShrink(bool mAutoShrink) { m_autoShrink = mAutoShrink; }

	//! Returns the current size of all registered arrays.
	int getArraySize() const { return m_size; }
};

//! Abstract base class for registered arrays.
/**
 * Defines the interface for event handling used by the registry.
 *
 * @tparam Registry The class which manages the registered keys. Must provide the functions defined in
 * class RegistryBase.
 */
template<class Registry>
class RegisteredArrayBase {
	using registry_type = Registry;
	using registration_iterator_type = typename Registry::registration_iterator_type;

	registration_iterator_type m_registration;
	const Registry* m_pRegistry = nullptr;

public:
	//! Creates a registered array associated with no registry.
	RegisteredArrayBase() = default;

	//! Creates a registered array associated with the same registry as \p copy.
	RegisteredArrayBase(const RegisteredArrayBase<Registry>& copy) { reregister(copy.m_pRegistry); }

	//! Moves the registration of \p move_from to this registered array.
	RegisteredArrayBase(RegisteredArrayBase<Registry>&& move_from) noexcept {
		moveRegister(move_from);
	}

	//! Assignment operator.
	RegisteredArrayBase& operator=(const RegisteredArrayBase<Registry>& copy) {
		reregister(copy.m_pRegistry);
		return *this;
	}

	//! Assignment operator (move semantics).
	RegisteredArrayBase& operator=(RegisteredArrayBase<Registry>&& move_from) noexcept {
		moveRegister(move_from);
		return *this;
	}

	//! Destructor.
	virtual ~RegisteredArrayBase() noexcept {
		if (m_pRegistry) {
			m_pRegistry->unregisterArray(m_registration);
		}
	}

	//! Resizes the registered array to \p size. The array will only shrink if \p shrink is \c true.
	virtual void resize(int size, bool shrink) = 0;

	//! Swaps the entries stored at \p index1 and \p index2.
	virtual void swapEntries(int index1, int index2) = 0;

	//! Copies the entry stored at \p oldIndex to \p newIndex.
	virtual void copyEntry(int newIndex, int oldIndex) = 0;

	//! Clears the array and associates it with no registry.
	void unregister() noexcept {
		resize(0, true);
		reregister(nullptr);
	}

protected:
	//! Associates the array with a new registry.
	void reregister(const Registry* registry) {
		if (m_pRegistry) {
			m_pRegistry->unregisterArray(m_registration);
		}
		m_pRegistry = registry;
		if (m_pRegistry != nullptr) {
			m_registration = m_pRegistry->registerArray(this);
		} else {
			m_registration = registration_iterator_type();
		}
	}

	//! Moves array registration from \p move_from to this array.
	void moveRegister(RegisteredArrayBase<Registry>& move_from) {
		if (m_pRegistry) {
			m_pRegistry->unregisterArray(m_registration);
		}
		m_pRegistry = move_from.m_pRegistry;
		m_registration = move_from.m_registration;
		move_from.m_pRegistry = nullptr;
		move_from.m_registration = registration_iterator_type();
		if (m_pRegistry != nullptr) {
			m_pRegistry->moveRegisterArray(m_registration, this);
		}
	}

public:
	//! Returns a pointer to the associated registry.
	const Registry* registeredAt() const { return m_pRegistry; }
};

//! Iterator for registered arrays.
/**
 * Provides an iterator for the key-value pairs stored in registered arrays.
 *
 * @tparam ArrayType The type of registered array.
 * @tparam KeyIterator An iterator for the keys in the registry. Determines the order of the key-value pairs.
 * @tparam isConst Whether the iterator allows modifying the data or not.
 */
template<class ArrayType, class KeyIterator, bool isConst = false>
class RegisteredArrayIterator {
public:
	using registry_type = typename ArrayType::registry_type;
	using key_type = typename ArrayType::key_type;
	using value_type = typename std::conditional<isConst, const typename ArrayType::value_type,
			typename ArrayType::value_type>::type;
	using array_pointer_type = typename std::conditional<isConst, const ArrayType*, ArrayType*>::type;

private:
	KeyIterator m_it;
	array_pointer_type m_array;

public:
	//! Creates a new iterator associated with no array.
	RegisteredArrayIterator() : m_it(), m_array(nullptr) { }

	//! Creates a new iterator.
	/**
	 * @param mIt An iterator pointing to the current key.
	 * @param mArray A pointer to the registered array containing the desired values.
	 */
	RegisteredArrayIterator(KeyIterator mIt, array_pointer_type mArray)
		: m_it(mIt), m_array(mArray) { }

	//! Returns the current key.
	key_type key() const { return *m_it; }

	//! Returns the value of key() in the registered array.
	value_type& value() const { return (*m_array)[*m_it]; }

	//! Returns the value of key() in the registered array.
	value_type& operator*() const { return (*m_array)[*m_it]; }

	//! Equality operator.
	bool operator==(const RegisteredArrayIterator<ArrayType, KeyIterator, isConst>& iter) const {
		return m_it == iter.m_it && m_array == iter.m_array;
	}

	//! Inequality operator.
	bool operator!=(const RegisteredArrayIterator<ArrayType, KeyIterator, isConst>& iter) const {
		return !operator==(iter);
	}

	//! Increment operator (prefix).
	RegisteredArrayIterator<ArrayType, KeyIterator, isConst>& operator++() {
		++m_it;
		return *this;
	}

	//! Increment operator (postfix).
	RegisteredArrayIterator<ArrayType, KeyIterator, isConst> operator++(int) {
		RegisteredArrayIterator<ArrayType, KeyIterator, isConst> iter = *this;
		++m_it;
		return iter;
	}

	//! Decrement operator (prefix).
	RegisteredArrayIterator<ArrayType, KeyIterator, isConst>& operator--() {
		--m_it;
		return *this;
	}

	//! Decrement operator (postfix).
	RegisteredArrayIterator<ArrayType, KeyIterator, isConst> operator--(int) {
		RegisteredArrayIterator<ArrayType, KeyIterator, isConst> iter = *this;
		--m_it;
		return iter;
	}
};

//! Registered arrays without default values.
/**
 * Registered arrays provide an efficient, constant-time mapping from indexed keys of a \a Registry to elements of
 * type \a Value. The storage automatically grows and shrinks when keys are added to or removed from the registry.
 * New values are initialized using the default constructor of \a Value.
 *
 * @tparam Registry The class which manages the registered keys. Must provide the functions defined in
 * class RegistryBase.
 * @tparam Value The type of the stored data.
 *
 * \sa RegistryBase, RegisteredArrayWithoutDefault, RegisteredArray
 */
template<class Registry, class Value>
class RegisteredArrayWithoutDefault : protected RegisteredArrayBase<Registry> {
protected:
	using key_iterator = typename Registry::iterator_type;
	using registered_array = RegisteredArrayWithoutDefault<Registry, Value>;
	using registered_array_base = RegisteredArrayBase<Registry>;

public:
	using registry_type = Registry;
	using key_type = typename Registry::key_type;
	using vector_type = std::vector<Value, OGDFAllocator<Value>>;
	using value_type = typename vector_type::value_type;
	using value_ref_type = typename vector_type::reference;
	using value_const_ref_type = typename vector_type::const_reference;

	using iterator = RegisteredArrayIterator<registered_array, key_iterator, false>;
	using const_iterator = RegisteredArrayIterator<registered_array, key_iterator, true>;

protected:
	vector_type m_data;

public:
	//! Creates a new registered array associated with no registry.
	RegisteredArrayWithoutDefault() = default;

	//! Creates a new registered array associated with \p registry.
	explicit RegisteredArrayWithoutDefault(const Registry* registry) {
		// during base class initialization, no part of the derived class exists, so this will always call our base init
		// so base classes should call their own init themselves
		registered_array::init(registry);
	}

	//! Associates the array with \p registry. All entries are initialized using the default constructor of \a Value.
	void init(const Registry* registry = nullptr) {
		if (registry == nullptr) {
			resize(0, true);
		} else {
			OGDF_ASSERT(registry->maxKeyIndex() < registry->getArraySize());
			resize(0, false);
			resize(registry->getArraySize(), true);
		}
		registered_array::reregister(registry);
	}

	//! Fills all entries with value \p x.
	void fill(value_const_ref_type x) { m_data.assign(m_data.size(), x); }

	//! Returns a const reference to the element associated with \p key.
	value_const_ref_type operator[](key_type key) const {
		OGDF_ASSERT(getRegistry().isKeyAssociated(key));
#ifdef OGDF_DEBUG
		return m_data.at(registeredAt()->keyToIndex(key));
#else
		return m_data[registeredAt()->keyToIndex(key)];
#endif
	}

	//! Returns a reference to the element associated with \p key.
	value_ref_type operator[](key_type key) {
		OGDF_ASSERT(getRegistry().isKeyAssociated(key));
#ifdef OGDF_DEBUG
		return m_data.at(registeredAt()->keyToIndex(key));
#else
		return m_data[registeredAt()->keyToIndex(key)];
#endif
	}

	//! Returns a const reference to the element associated with \p key.
	value_const_ref_type operator()(key_type key) const {
		OGDF_ASSERT(getRegistry().isKeyAssociated(key));
#ifdef OGDF_DEBUG
		return m_data.at(registeredAt()->keyToIndex(key));
#else
		return m_data[registeredAt()->keyToIndex(key)];
#endif
	}

	//! Returns a reference to the element associated with \p key.
	value_ref_type operator()(key_type key) {
		OGDF_ASSERT(getRegistry().isKeyAssociated(key));
#ifdef OGDF_DEBUG
		return m_data.at(registeredAt()->keyToIndex(key));
#else
		return m_data[registeredAt()->keyToIndex(key)];
#endif
	}

	//! Returns a const reference to the element with index \p idx.
	value_const_ref_type operator[](int idx) const {
#ifdef OGDF_DEBUG
		return m_data.at(idx);
#else
		return m_data[idx];
#endif
	}

	//! Returns a reference to the element with index \p idx.
	value_ref_type operator[](int idx) {
#ifdef OGDF_DEBUG
		return m_data.at(idx);
#else
		return m_data[idx];
#endif
	}

	//! Returns an iterator to the first key-value pair in the array.
	iterator begin() {
		using std::begin;
		return iterator(begin(getRegistry()), this);
	}

	//! Returns a const iterator to the first key-value pair in the array.
	const_iterator begin() const {
		using std::begin;
		return const_iterator(begin(getRegistry()), this);
	}

	//! Returns a const iterator to the first key-value pair in the array.
	const_iterator cbegin() const {
		using std::begin;
		return const_iterator(begin(getRegistry()), this);
	}

	//! Returns the past-the-end iterator of the array.
	iterator end() {
		using std::end;
		return iterator(end(getRegistry()), this);
	}

	//! Returns the const past-the-end iterator of the array.
	const_iterator end() const {
		using std::end;
		return const_iterator(end(getRegistry()), this);
	}

	//! Returns the const past-the-end iterator of the array.
	const_iterator cend() const {
		using std::end;
		return const_iterator(end(getRegistry()), this);
	}

	using registered_array_base::registeredAt;

	//! Returns true iff the array is associated with a registry.
	bool valid() const {
		OGDF_ASSERT(registeredAt() == nullptr || registeredAt()->maxKeyIndex() < 0
				|| ((size_t)registeredAt()->maxKeyIndex()) < m_data.size());
		return registeredAt();
	}

protected:
	//! Returns a reference to the associated registry.
	inline const Registry& getRegistry() const {
		OGDF_ASSERT(registeredAt());
		OGDF_ASSERT(valid());
		return *registeredAt();
	}

	void resize(int size, bool shrink) override {
		m_data.resize(size);
		if (shrink) {
			m_data.shrink_to_fit();
		}
	}

	void swapEntries(int index1, int index2) override {
		std::swap(m_data.at(index1), m_data.at(index2));
	}

	//! This operation is not supported for registered arrays without default.
	void copyEntry(int toIndex, int fromIndex) override {
		// silently ignored
	}
};

//! Registered arrays with default values.
/**
 * Extends the functionality of RegisteredArrayWithoutDefault by adding the possibility to set a specific default value
 * for new keys added to the registry.
 *
 * \pre Type \a Value must be copy-constructible.
 *
 * @tparam Registry The class which manages the registered keys. Must provide the functions defined in
 * class RegistryBase.
 * @tparam Value The type of the stored data.
 *
 * \sa RegistryBase, RegisteredArrayWithoutDefault, RegisteredArray
 */
template<class Registry, class Value>
class RegisteredArrayWithDefault : public RegisteredArrayWithoutDefault<Registry, Value> {
	using RA = RegisteredArrayWithoutDefault<Registry, Value>;
	Value m_default;

public:
	//! Creates a new registered array associated with no registry and a default-constructed default value.
	explicit RegisteredArrayWithDefault() : RA(), m_default() {};

	//! Creates a new registered array associated with \p registry and a default-constructed default value.
	explicit RegisteredArrayWithDefault(const Registry* registry) : RA(), m_default() {
		// call init from here, as our virtual override of init is not available during initialization of the base class
		RA::init(registry);
	};

	//! Creates a new registered array associated with no registry and default value \p def.
	explicit RegisteredArrayWithDefault(const Value& def) : RA(), m_default(def) {};

	//! Creates a new registered array associated with \p registry and default value \p def.
	explicit RegisteredArrayWithDefault(const Registry* registry, const Value& def)
		: RA(), m_default(def) {
		// call init from here, as our virtual override of init is not available during initialization of the base class
		RA::init(registry);
	};

	//! Creates a new registered array associated with no registry and default value \p def.
	explicit RegisteredArrayWithDefault(Value&& def) : RA(), m_default(std::forward<Value>(def)) {};

	//! Creates a new registered array associated with \p registry and default value \p def.
	explicit RegisteredArrayWithDefault(const Registry* registry, Value&& def)
		: RA(), m_default(std::forward<Value>(def)) {
		// call init from here, as our virtual override of init is not available during initialization of the base class
		RA::init(registry);
	};

	//! Sets a new default value for new keys.
	void setDefault(Value&& def) { m_default = std::forward<Value>(def); }

	//! Sets a new default value for new keys.
	void setDefault(const Value& def) { m_default = def; }

	//! Returns the current default value for new keys.
	const Value& getDefault() const { return m_default; }

	//! Returns the current default value for new keys.
	Value& getDefault() { return m_default; }

	//! Overwrites all values with the current default value.
	void fillWithDefault() { RA::m_data.assign(RA::getRegistry().getArraySize(), m_default); }

protected:
	//! Copies the entry stored at \p oldIndex to \p newIndex.
	void copyEntry(int toIndex, int fromIndex) override {
		RA::m_data.at(toIndex) = RA::m_data.at(fromIndex);
	}

	void resize(int size, bool shrink) override {
		RA::m_data.resize(size, m_default);
		if (shrink) {
			RA::m_data.shrink_to_fit();
		}
	}
};

//! Dynamic arrays indexed with arbitrary keys.
/**
 * Registered arrays provide an efficient, constant-time mapping from indexed keys of a \a Registry to elements of
 * type \a Value. The storage automatically grows and shrinks when keys are added to or removed from the registry.
 *
 * \warning When the array grows or shrinks, all pointers to its entries become invalid.
 *
 * @tparam Registry The class which manages the registered keys. Must provide the functions defined in
 * class RegistryBase.
 * @tparam Value The type of the stored data.
 * @tparam WithDefault Determines whether the registered array inherits from RegisteredArrayWithDefault
 * or RegisteredArrayWithoutDefault. With \a WithDefault \= \c true, the array can be initialized with specific default
 * values, but this requires \a Value to be a copy-constructible type. With \a WithDefault \= \c false, the array uses
 * the default constructor of \a Value to initialize new storage.
 * @tparam Base The class that manages multiple related registries. \a Base must be convertible
 * to \a Registry. If only one such registry exists, \a Base and \a Registry can be the
 * same class (i.e. \a Base directly inherits from class RegistryBase)
 *
 * ### Class Interaction
 * - **Key**
 *
 * 	Used to index the registered array. Every key must have a unique non-negative index. Must either provide a public
 * 	method called \c index() or the template function keyToIndex() must offer a specialization for \a Key to give access
 * 	to its index.
 *
 * - **Value**
 *
 * 	The type of the elements stored in the array.
 *
 * - **RegistryBase**
 *
 * 	Defines the interface for the \a Registry.
 *
 * - **Registry**
 *
 * 	Implements the abstract functionality defined in RegistryBase. Manages the objects of type \a Key and stores a list
 * 	of associated registered arrays for \a Key. Determines the growth rate of the arrays. When keys are added or
 * 	removed, the functions RegistryBase::keyAdded(), RegistryBase::keyRemoved(), and RegistryBase::keysCleared()
 * 	should be called so that the size of all arrays will be adjusted accordingly.
 *
 * - **RegisteredArrayBase**
 *
 * 	Abstract base class for all registered arrays. The \a Registry communicates with its registered arrays using
 * 	this interface.
 *
 * - **RegisteredArrayWithoutDefault**
 *
 * 	Provides the core functionality for accessing the values stored in the array. New entries are initialized using the
 * 	default constructor of type \a Value.
 *
 * - **RegisteredArrayWithDefault**
 *
 * 	Extends the functionality of RegisteredArrayWithoutDefault by adding the possibility to set a specific default value
 * 	for new keys added to the registry. This requires type \a Value to be copy-constructible.
 *
 * - **RegisteredArray**
 *
 * 	Used in user code. Inherits from RegisteredArrayWithoutDefault or RegisteredArrayWithDefault, depending on the
 * 	template parameter \a WithDefault.
 *
 *
 * ### Example Setup
 *
 * A simple registry that only allows addition of keys:
 * \code
 * class ExampleKey {
 * 	int m_index;
 *
 * public:
 * 	explicit ExampleKey(int index) : m_index(index) {}
 * 	int index() const { return m_index; }
 * };
 *
 * class ExampleRegistry : public RegistryBase<ExampleKey *, ExampleRegistry> {
 * 	std::list<std::unique_ptr<ExampleKey>> m_keys;
 *
 * public:
 * 	ExampleKey *newKey() {
 * 		m_keys.push_back(std::unique_ptr<ExampleKey>(new ExampleKey(m_keys.size())));
 * 		keyAdded(m_keys.back().get());
 * 		return m_keys.back().get();
 * 	}
 *
 * 	bool isKeyAssociated(ExampleKey *key) const override { return true; }
 * 	int maxKeyIndex() const override { return m_keys.size() - 1; }
 * 	int calculateArraySize() const override { return calculateTableSize(m_keys.size()); }
 * 	void begin() const override {}
 * 	void end() const override {}
 * };
 * \endcode
 * With this setup, registering an array and modifying its values works as follows:
 * \code
 * ExampleRegistry G;
 * RegisteredArray<ExampleRegistry, int> R(G);
 * ExampleKey *key = G.newKey();
 * R[key] = 42;
 * \endcode
 *
 * \extends RegisteredArrayWithDefault
 *
 * \sa RegisteredArrayWithoutDefault, RegisteredArrayWithDefault, RegistryBase, NodeArray
 */
template<class Registry, class Value, bool WithDefault = true, class Base = Registry>
class RegisteredArray
	: public std::conditional<WithDefault, RegisteredArrayWithDefault<Registry, Value>,
			  RegisteredArrayWithoutDefault<Registry, Value>>::type {
	using RA = typename std::conditional<WithDefault, RegisteredArrayWithDefault<Registry, Value>,
			RegisteredArrayWithoutDefault<Registry, Value>>::type;

	static inline const Registry* cast(const Base* base) {
		if (base != nullptr) {
			// unpack the pointer to invoke the conversion operator
			return &((const Registry&)*base);
		} else {
			return nullptr;
		}
	}

public:
	//! Creates a new registered array associated with no registry.
	RegisteredArray() : RA() {};

	//! Creates a new registered array associated with the matching registry of \p base.
	explicit RegisteredArray(const Base& base) : RA(cast(&base)) {};

	/**
	 * Creates a new registered array associated with the matching registry of \p base and initializes all values with \p def.
	 *
	 * \remarks This constructor is only available with \a WithDefault \= \c true.
	 */
	RegisteredArray(const Base& base, const Value& def) : RA(cast(&base), def) {};

	//! Creates a new registered array associated with the matching registry of \p base.
	explicit RegisteredArray(const Base* base) : RA(cast(base)) {};

	/**
	 * Creates a new registered array associated with the matching registry of \p base and initializes all values with \p def.
	 *
	 * \remarks This constructor is only available with \a WithDefault \= \c true.
	 */
	RegisteredArray(const Base* base, const Value& def) : RA(cast(base), def) {};

	//! Reinitializes the array. Associates the array with the matching registry of \p base.
	void init(const Base* base = nullptr) { RA::init(cast(base)); }

	//! Reinitializes the array. Associates the array with the matching registry of \p base.
	void init(const Base& base) { RA::init(cast(&base)); }

	/**
	 * Reinitializes the array with default value \p new_default. Associates the array with the matching registry of \p base.
	 *
	 * \remarks This method is only available with \a WithDefault \= \c true.
	 */
	void init(const Base& base, const Value& new_default) {
		RA::setDefault(new_default);
		RA::init(cast(&base));
	}

	void init(const Base* base, const Value& new_default) {
		RA::setDefault(new_default);
		RA::init(cast(base));
	}
};

//! Specialization to work around vector<bool>.
template<class Registry, bool WithDefault, class Base>
class RegisteredArray<Registry, bool, WithDefault, Base>
	: public RegisteredArray<Registry, unsigned char, WithDefault, Base> {
	using RA = RegisteredArray<Registry, unsigned char, WithDefault, Base>;
	using BRA = RegisteredArray<Registry, bool, WithDefault, Base>;

public:
	using RA::RA;

	using key_type = typename RA::key_type;
	using value_type = bool;
	using value_const_ref_type = const bool&;
	using value_ref_type = bool&;
	using iterator = RegisteredArrayIterator<BRA, typename RA::key_iterator, false>;
	using const_iterator = RegisteredArrayIterator<BRA, typename RA::key_iterator, true>;

	value_const_ref_type operator[](key_type key) const {
		return reinterpret_cast<value_const_ref_type>(RA::operator[](key));
	}

	value_ref_type operator[](key_type key) {
		return reinterpret_cast<value_ref_type>(RA::operator[](key));
	}

	value_const_ref_type operator()(key_type key) const {
		return reinterpret_cast<value_const_ref_type>(RA::operator()(key));
	}

	value_ref_type operator()(key_type key) {
		return reinterpret_cast<value_ref_type>(RA::operator()(key));
	}

	value_const_ref_type operator[](int idx) const {
		return reinterpret_cast<value_const_ref_type>(RA::operator[](idx));
	}

	value_ref_type operator[](int idx) {
		return reinterpret_cast<value_ref_type>(RA::operator[](idx));
	}

	value_ref_type getDefault() { return reinterpret_cast<value_ref_type>(RA::getDefault()); }

	value_const_ref_type getDefault() const {
		return reinterpret_cast<value_const_ref_type>(RA::getDefault());
	}

	iterator begin() {
		using std::begin;
		return iterator(begin(RA::getRegistry()), this);
	}

	const_iterator begin() const {
		using std::begin;
		return const_iterator(begin(RA::getRegistry()), this);
	}

	const_iterator cbegin() const {
		using std::begin;
		return const_iterator(begin(RA::getRegistry()), this);
	}

	iterator end() {
		using std::end;
		return iterator(end(RA::getRegistry()), this);
	}

	const_iterator end() const {
		using std::end;
		return const_iterator(end(RA::getRegistry()), this);
	}

	const_iterator cend() const {
		using std::end;
		return const_iterator(end(RA::getRegistry()), this);
	}
};
}

#define OGDF_DECL_REG_ARRAY(NAME)                              \
	template<typename Value, bool WithDefault = true>          \
	using NAME = OGDF_DECL_REG_ARRAY_TYPE(Value, WithDefault); \
	template<typename Value>                                   \
	using NAME##P = NAME<std::unique_ptr<Value>, false>;
