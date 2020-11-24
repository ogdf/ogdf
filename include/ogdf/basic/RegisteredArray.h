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

namespace ogdf {
template<typename Registry>
class RegisteredArrayBase;

static constexpr int MIN_TABLE_SIZE = (1 << 4);

int calculateTableSize(int actualCount);

template<typename Key, typename Registry, typename Iterator = void>
class RegistryBase {
public:
	using registered_array_type = RegisteredArrayBase<Registry>;
	using key_type = Key;
	using registry_type = Registry;
	using iterator_type = Iterator;
	using registration_list_type = std::list<registered_array_type*>;
	using registration_iterator_type = typename registration_list_type::iterator;

private:
	mutable registration_list_type m_registeredArrays;
	bool m_lazy = false;
	bool m_auto_shrink = false;
	int m_size = 0;

#ifndef OGDF_MEMORY_POOL_NTS
	mutable std::mutex m_mutexRegArrays;
#endif

protected:
	RegistryBase() = default;

public:
	virtual ~RegistryBase() { unregisterArrays(); }

	[[nodiscard]] virtual registration_iterator_type registerArray(
			registered_array_type* pArray) const {
#ifndef OGDF_MEMORY_POOL_NTS
		std::lock_guard<std::mutex> guard(m_mutexRegArrays);
#endif
		return m_registeredArrays.emplace(m_registeredArrays.end(), pArray);
	}

	virtual void unregisterArray(registration_iterator_type it) const {
#ifndef OGDF_MEMORY_POOL_NTS
		std::lock_guard<std::mutex> guard(m_mutexRegArrays);
#endif
		m_registeredArrays.erase(it);
	}

	virtual void moveRegisterArray(registration_iterator_type it,
			registered_array_type* pArray) const {
#ifndef OGDF_MEMORY_POOL_NTS
		std::lock_guard<std::mutex> guard(m_mutexRegArrays);
#endif
		*it = pArray;
	}

	[[nodiscard]] virtual bool isKeyAssociated(Key key) const = 0;

	[[nodiscard]] virtual int keyToIndex(Key key) const = 0;

	[[nodiscard]] virtual int maxKeyIndex() const = 0;

	[[nodiscard]] virtual int arraySize() const = 0;

	[[nodiscard]] virtual Iterator begin() const = 0;

	[[nodiscard]] virtual Iterator end() const = 0;

	virtual void keyAdded(Key key) {
		// TODO use m_lazy, m_auto_shrink and maxKeyIndex() / arraySize() to call resizeArrays(...)
		// TODO notify observers
	}

	virtual void keyRemoved(Key key) {
		// TODO use m_lazy, m_auto_shrink and maxKeyIndex() / arraySize() to call resizeArrays(...)
		// TODO notify observers
	}

	virtual void keysCleared() {
		// TODO use m_lazy, m_auto_shrink and maxKeyIndex() / arraySize() to call resizeArrays(...)
		// TODO notify observers
	}

	void resizeArrays(int size) {
		m_size = size;
		for (registered_array_type* ab : m_registeredArrays) {
			ab->resize(size, m_auto_shrink);
		}
	}

	void resizeArrays(int size, bool shrink) {
		m_size = size;
		for (registered_array_type* ab : m_registeredArrays) {
			ab->resize(size, shrink);
		}
	}

	void swapArrayEntries(int newIndex, int oldIndex) {
		for (registered_array_type* ab : m_registeredArrays) {
			ab->swapEntry(newIndex, oldIndex);
		}
	}

	void unregisterArrays() {
		while (!m_registeredArrays.empty()) {
#ifdef OGDF_DEBUG
			auto size = m_registeredArrays.size();
#endif
			m_registeredArrays.front()->unregister();
			OGDF_ASSERT(m_registeredArrays.size() < size);
		}
	}
};

// TODO Add Observers
//	template<class Registry>
//	class RegistryObserver {
//		using key_type = typename Registry::key_type;
//
//		virtual void registered(Registry &registry) = 0;
//
//		virtual void unregistered(Registry &registry) = 0;
//
//		virtual void keyAdded(Registry &registry, key_type key) = 0;
//
//		virtual void keyRemoved(Registry &registry, key_type key) = 0;
//
//		virtual void keysCleared(Registry &registry) = 0;
//	};

template<class Registry>
class RegisteredArrayBase {
	using registry_type = Registry;
	using registration_iterator_type = typename Registry::registration_iterator_type;

	registration_iterator_type m_registration;
	const Registry* m_pRegistry = nullptr;

protected:
	RegisteredArrayBase() = default;

	// TODO copy / move assignment and constructors
	//		RegisteredArrayBase(const RegisteredArrayBase<Registry> &copy) : m_pRegistry(copy.m_pRegistry) {
	//			if (m_pRegistry)
	//				m_registration = m_pRegistry->registerArray(this);
	//		}
	//
	//		RegisteredArrayBase(RegisteredArrayBase<Registry> &&move_from) noexcept
	//				: m_registration(move_from.m_registration), m_pRegistry(move_from.m_pRegistry) {
	//			m_pRegistry->moveRegisterArray(m_registration, this);
	//			move_from.m_pRegistry = nullptr;
	//			move_from.m_registration = ListIterator<RegisteredArrayBase<Registry> *>();
	//		}

public:
	virtual ~RegisteredArrayBase() {
		if (m_pRegistry) {
			m_pRegistry->unregisterArray(m_registration);
		}
	}

	virtual void resize(int size, bool shrink) = 0;

	virtual void swapEntries(int newIndex, int oldIndex) = 0;

	virtual void unregister() {
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
		move_from.m_registration = nullptr;
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
	//		using vector_type = typename std::conditional<std::is_same<Value, bool>::value,
	//				std::vector<unsigned char>, std::vector<Value>>::type;
	using vector_type = std::vector<Value>;
	using value_type = typename vector_type::value_type;
	using value_ref_type = typename vector_type::reference;
	using value_const_ref_type = typename vector_type::const_reference;

	using iterator = RegisteredArrayIterator<registered_array, key_iterator, false>;
	using const_iterator = RegisteredArrayIterator<registered_array, key_iterator, true>;

protected:
	vector_type m_data;

public:
	explicit RegisteredArray(const Registry* registry = nullptr) { init(registry); }

	// TODO copy / move assignment and constructors
	//		registered_array &operator=(const registered_array &A) {
	//			array::operator=(A);
	//			m_default_value = A.m_default_value;
	//			registered_array_base::reregister(A.registeredAt());
	//			return *this;
	//		}
	//
	//		registered_array &operator=(registered_array &&a) {
	//			array::operator=(std::move(a));
	//			m_default_value = a.m_default_value;
	//			moveRegister(a.registeredAt());
	//			registered_array_base::reregister(a.registeredAt());
	//			a.init();
	//			return *this;
	//		}

	~RegisteredArray() override = default;

	void init(const Registry* registry = nullptr) {
		if (registry == nullptr) {
			resize(0, true);
		} else {
			resize(registry->arraySize(), true);
		}
		registered_array::reregister(registry);
	}

	void fill(value_const_ref_type x) { m_data.assign(getRegistry().arraySize(), x); }

	value_const_ref_type operator[](key_type key) const {
		OGDF_ASSERT(getRegistry().isKeyAssociated(key));
		return m_data[getRegistry().keyToIndex(key)];
	}

	value_ref_type operator[](key_type key) {
		OGDF_ASSERT(getRegistry().isKeyAssociated(key));
		return m_data[getRegistry().keyToIndex(key)];
	}

	value_const_ref_type operator()(key_type key) const {
		OGDF_ASSERT(getRegistry().isKeyAssociated(key));
		return m_data[getRegistry().keyToIndex(key)];
	}

	value_ref_type operator()(key_type key) {
		OGDF_ASSERT(getRegistry().isKeyAssociated(key));
		return m_data[getRegistry().keyToIndex(key)];
	}

	value_const_ref_type operator[](int idx) const { return m_data[idx]; }

	value_ref_type operator[](int idx) { return m_data[idx]; }

	iterator begin() { return iterator(getRegistry().begin(), this); }

	const_iterator begin() const { return const_iterator(getRegistry().begin(), this); }

	const_iterator cbegin() const { return const_iterator(getRegistry().begin(), this); }

	iterator end() { return iterator(getRegistry().end(), this); }

	const_iterator end() const { return const_iterator(getRegistry().end(), this); }

	const_iterator cend() const { return const_iterator(getRegistry().end(), this); }

	bool valid() const { return registered_array_base::registeredAt(); }

	using registered_array_base::registeredAt;

protected:
	inline const Registry& getRegistry() const {
		OGDF_ASSERT(registeredAt());
		return *registeredAt();
	}

	void resize(int size, bool shrink) override {
		m_data.resize(size);
		if (shrink) {
			m_data.shrink_to_fit();
		}
	}

	void swapEntries(int newIndex, int oldIndex) override {
		std::swap(m_data[newIndex], m_data[oldIndex]);
	}
};

template<class Registry, class Value>
class RegisteredArrayWithDefault : public RegisteredArray<Registry, Value> {
	using RA = RegisteredArray<Registry, Value>;
	Value m_default;

public:
	// TODO use varargs or move default?
	//		explicit RegisteredArrayWithDefault(const Registry *registry, Value &&def)
	//				: RA(registry), m_default(std::forward<Value>(def)) {};

	explicit RegisteredArrayWithDefault(const Registry* registry, const Value& def)
		: RA(registry), m_default(def) {};

	//		explicit RegisteredArrayWithDefault(const Registry *registry, typename RA::value_const_ref_type def)
	//				: RA(registry), m_default(def) {};

	~RegisteredArrayWithDefault() override = default;

	//		void setDefault(Value &&def) {
	//			m_default = std::move<Value>(def);
	//		}

	void setDefault(const Value& def) { m_default = def; }

	//		void setDefault(typename RA::value_const_ref_type def) {
	//			m_default = def;
	//		}

	void fillWithDefault() { RA::m_data.assign(RA::getRegistry().arraySize(), m_default); }

protected:
	void resize(int size, bool shrink) override {
		RA::m_data.resize(size, m_default);
		if (shrink) {
			RA::m_data.shrink_to_fit();
		}
	}
};
}
