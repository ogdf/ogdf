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

#include <ogdf/basic/Math.h>

#include <list>

#ifndef OGDF_MEMORY_POOL_NTS

#	include <mutex>

#endif

#define OGDF_NODISCARD

// TODO define OGDF_NODISCARD properly and in a sensible place

namespace ogdf {
template<typename Registry>
class RegisteredArrayBase;

static constexpr int MIN_TABLE_SIZE = (1 << 4);

int calculateTableSize(int actualCount);

template<class Registry>
class RegistryObserver {
public:
	using key_type = typename Registry::key_type;

	virtual void registered(Registry* registry) = 0;

	virtual void unregistered(Registry* registry) noexcept = 0;

	virtual void keyAdded(Registry* registry, key_type key) = 0;

	virtual void keyRemoved(Registry* registry, key_type key) = 0;

	virtual void keysCleared(Registry* registry) = 0;
};

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
	virtual ~RegistryBase() noexcept { unregisterArrays(); }

	RegistryBase(const RegistryBase& copy) = delete;

	RegistryBase(RegistryBase&& move) noexcept = delete;

	RegistryBase& operator=(const RegistryBase& other) = delete;

	RegistryBase& operator=(RegistryBase&& other) noexcept = delete;

	OGDF_NODISCARD registration_iterator_type registerArray(registered_array_type* pArray) const {
#ifndef OGDF_MEMORY_POOL_NTS
		std::lock_guard<std::mutex> guard(m_mutexRegArrays);
#endif
		return m_registeredArrays.emplace(m_registeredArrays.end(), pArray);
	}

	void unregisterArray(registration_iterator_type it) const noexcept {
#ifndef OGDF_MEMORY_POOL_NTS
		std::lock_guard<std::mutex> guard(m_mutexRegArrays);
#endif
		m_registeredArrays.erase(it);
	}

	void moveRegisterArray(registration_iterator_type it, registered_array_type* pArray) const {
#ifndef OGDF_MEMORY_POOL_NTS
		std::lock_guard<std::mutex> guard(m_mutexRegArrays);
#endif
		*it = pArray;
	}

	OGDF_NODISCARD virtual bool isKeyAssociated(Key key) const = 0;

	OGDF_NODISCARD virtual int keyToIndex(Key key) const = 0;

	OGDF_NODISCARD virtual int maxKeyIndex() const = 0;

	OGDF_NODISCARD virtual int calculateArraySize() const = 0;

	OGDF_NODISCARD virtual Iterator begin() const = 0;

	OGDF_NODISCARD virtual Iterator end() const = 0;

	void keyAdded(Key key) { resizeArrays(); }

	void keyRemoved(Key key) { resizeArrays(); }

	void keysCleared() { resizeArrays(0); }

	void resizeArrays() { resizeArrays(calculateArraySize()); }

	void resizeArrays(int size) { resizeArrays(size, m_autoShrink); }

	void resizeArrays(int size, bool shrink) {
		if (size == m_size) {
			return;
		}
		m_size = size = max(size, 0);
		for (registered_array_type* ab : m_registeredArrays) {
			ab->resize(size, shrink);
		}
	}

	void shrinkArrays(int size) {
		m_size = size = max(size, 0);
		for (registered_array_type* ab : m_registeredArrays) {
			ab->resize(size, true);
		}
	}

	void swapArrayEntries(int newIndex, int oldIndex) {
		for (registered_array_type* ab : m_registeredArrays) {
			ab->swapEntries(newIndex, oldIndex);
		}
	}

	void copyArrayEntries(int toIndex, int fromIndex) {
		for (registered_array_type* ab : m_registeredArrays) {
			ab->copyEntry(toIndex, fromIndex);
		}
	}

	void unregisterArrays() noexcept {
		while (!m_registeredArrays.empty()) {
#ifdef OGDF_DEBUG
			auto size = m_registeredArrays.size();
#endif
			m_registeredArrays.front()->unregister();
			OGDF_ASSERT(m_registeredArrays.size() < size);
		}
	}

	const registration_list_type& getRegisteredArrays() const { return m_registeredArrays; }

	bool isAutoShrink() const { return m_autoShrink; }

	void setAutoShrink(bool mAutoShrink) { m_autoShrink = mAutoShrink; }

	int getArraySize() const { return m_size; }
};

template<class Registry>
class RegisteredArrayBase {
	using registry_type = Registry;
	using registration_iterator_type = typename Registry::registration_iterator_type;

	registration_iterator_type m_registration;
	const Registry* m_pRegistry = nullptr;

public:
	RegisteredArrayBase() = default;

	RegisteredArrayBase(const RegisteredArrayBase<Registry>& copy) { reregister(copy.m_pRegistry); }

	RegisteredArrayBase(RegisteredArrayBase<Registry>&& move_from) noexcept {
		moveRegister(move_from);
	}

	RegisteredArrayBase& operator=(const RegisteredArrayBase<Registry>& copy) {
		reregister(copy.m_pRegistry);
		return *this;
	}

	RegisteredArrayBase& operator=(RegisteredArrayBase<Registry>&& move_from) noexcept {
		moveRegister(move_from);
		return *this;
	}

	virtual ~RegisteredArrayBase() noexcept {
		if (m_pRegistry) {
			m_pRegistry->unregisterArray(m_registration);
		}
	}

	virtual void resize(int size, bool shrink) = 0;

	virtual void swapEntries(int newIndex, int oldIndex) = 0;

	virtual void copyEntry(int newIndex, int oldIndex) = 0;

	void unregister() noexcept {
		resize(0, true);
		reregister(nullptr);
	}

protected:
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
	const Registry* registeredAt() const { return m_pRegistry; }
};

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
	RegisteredArrayIterator(KeyIterator mIt, array_pointer_type mArray)
		: m_it(mIt), m_array(mArray) { }

	key_type key() const { return *m_it; }

	value_type& value() const { return (*m_array)[*m_it]; }

	value_type& operator*() const { return (*m_array)[*m_it]; }

	bool operator==(const RegisteredArrayIterator<ArrayType, KeyIterator, isConst>& iter) const {
		return m_it == iter.m_it && m_array == iter.m_array;
	}

	bool operator!=(const RegisteredArrayIterator<ArrayType, KeyIterator, isConst>& iter) const {
		return !operator==(iter);
	}

	RegisteredArrayIterator<ArrayType, KeyIterator, isConst>& operator++() {
		++m_it;
		return *this;
	}

	RegisteredArrayIterator<ArrayType, KeyIterator, isConst> operator++(int) {
		RegisteredArrayIterator<ArrayType, KeyIterator, isConst> iter = *this;
		++m_it;
		return iter;
	}

	RegisteredArrayIterator<ArrayType, KeyIterator, isConst>& operator--() {
		--m_it;
		return *this;
	}

	RegisteredArrayIterator<ArrayType, KeyIterator, isConst> operator--(int) {
		RegisteredArrayIterator<ArrayType, KeyIterator, isConst> iter = *this;
		--m_it;
		return iter;
	}
};

template<class Registry, class Value>
class RegisteredArray : protected RegisteredArrayBase<Registry> {
protected:
	using key_iterator = typename Registry::iterator_type;
	using registered_array = RegisteredArray<Registry, Value>;
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
	RegisteredArray() = default;

	explicit RegisteredArray(const Registry* registry) {
		// during base class initialization, no part of the derived class exists, so this will always call our base init
		// so base classes should call their own init themselves
		registered_array::init(registry);
	}

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

	void fill(value_const_ref_type x) { m_data.assign(m_data.size(), x); }

	value_const_ref_type operator[](key_type key) const {
		OGDF_ASSERT(getRegistry().isKeyAssociated(key));
#ifdef OGDF_DEBUG
		return m_data.at(getRegistry().keyToIndex(key));
#else
		return m_data[getRegistry().keyToIndex(key)];
#endif
	}

	value_ref_type operator[](key_type key) {
		OGDF_ASSERT(getRegistry().isKeyAssociated(key));
#ifdef OGDF_DEBUG
		return m_data.at(getRegistry().keyToIndex(key));
#else
		return m_data[getRegistry().keyToIndex(key)];
#endif
	}

	value_const_ref_type operator()(key_type key) const {
		OGDF_ASSERT(getRegistry().isKeyAssociated(key));
#ifdef OGDF_DEBUG
		return m_data.at(getRegistry().keyToIndex(key));
#else
		return m_data[getRegistry().keyToIndex(key)];
#endif
	}

	value_ref_type operator()(key_type key) {
		OGDF_ASSERT(getRegistry().isKeyAssociated(key));
#ifdef OGDF_DEBUG
		return m_data.at(getRegistry().keyToIndex(key));
#else
		return m_data[getRegistry().keyToIndex(key)];
#endif
	}

	value_const_ref_type operator[](int idx) const {
#ifdef OGDF_DEBUG
		return m_data.at(idx);
#else
		return m_data[idx];
#endif
	}

	value_ref_type operator[](int idx) {
#ifdef OGDF_DEBUG
		return m_data.at(idx);
#else
		return m_data[idx];
#endif
	}

	iterator begin() { return iterator(getRegistry().begin(), this); }

	const_iterator begin() const { return const_iterator(getRegistry().begin(), this); }

	const_iterator cbegin() const { return const_iterator(getRegistry().begin(), this); }

	iterator end() { return iterator(getRegistry().end(), this); }

	const_iterator end() const { return const_iterator(getRegistry().end(), this); }

	const_iterator cend() const { return const_iterator(getRegistry().end(), this); }

	bool valid() const {
		OGDF_ASSERT(registeredAt() == nullptr || registeredAt()->maxKeyIndex() < 0
				|| ((size_t)registeredAt()->maxKeyIndex()) < m_data.size());
		return registered_array_base::registeredAt();
	}

	using registered_array_base::registeredAt;

protected:
	inline const Registry& getRegistry() const {
		OGDF_ASSERT(registeredAt());
		OGDF_ASSERT(registeredAt()->maxKeyIndex() < ((int)m_data.size()));
		return *registeredAt();
	}

	void resize(int size, bool shrink) override {
		m_data.resize(size);
		if (shrink) {
			m_data.shrink_to_fit();
		}
	}

	void swapEntries(int newIndex, int oldIndex) override {
		std::swap(m_data.at(newIndex), m_data.at(oldIndex));
	}

	void copyEntry(int toIndex, int fromIndex) override {
		m_data.at(toIndex) = m_data.at(fromIndex);
	}
};

template<class Registry, class Value>
class RegisteredArrayWithDefault : public RegisteredArray<Registry, Value> {
	using RA = RegisteredArray<Registry, Value>;
	Value m_default;

public:
	explicit RegisteredArrayWithDefault(Value&& def) : RA(), m_default(std::forward<Value>(def)) {};

	explicit RegisteredArrayWithDefault(const Registry* registry, Value&& def)
		: RA(), m_default(std::forward<Value>(def)) {
		// call init from here, as our virtual override of init is not available during initialization of the base class
		RA::init(registry);
	};

	explicit RegisteredArrayWithDefault(const Value& def) : RA(), m_default(def) {};

	explicit RegisteredArrayWithDefault(const Registry* registry, const Value& def)
		: RA(), m_default(def) {
		// call init from here, as our virtual override of init is not available during initialization of the base class
		RA::init(registry);
	};

	void setDefault(Value&& def) { m_default = std::forward<Value>(def); }

	void setDefault(const Value& def) { m_default = def; }

	const Value& getDefault() const { return m_default; }

	Value& getDefault() { return m_default; }

	void fillWithDefault() { RA::m_data.assign(RA::getRegistry().getArraySize(), m_default); }

protected:
	void resize(int size, bool shrink) override {
		RA::m_data.resize(size, m_default);
		if (shrink) {
			RA::m_data.shrink_to_fit();
		}
	}
};
}
