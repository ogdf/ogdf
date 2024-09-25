/** \file
 * \brief Data structure for two-dimensional mappings that are sparse in the second dimension. TODO should be moved to a central location.
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
#include <ogdf/basic/internal/copy_move.h>

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace ogdf {

//! A wrapper around std::map that uses a constant-size array (or only a single value) plus linear search until the map grows too big.
/**
 * Automatically grows from storing a single value, over using an std::array of at most \p array_max values, to using an unbounded std::map.
 * This is used for managing the second dimension of RegisteredMultiArray for a single value in the first dimension.
 */
template<typename Key2, typename Value, int array_max>
struct UsuallySmallMap {
	using ValuePairType = std::pair<Key2, Value>;
	using ValueArrayType = std::array<ValuePairType, array_max>;
	using ValueMapType = std::unordered_map<Key2, Value>;

	uint8_t m_size = 0;
	void* m_value = nullptr;


	UsuallySmallMap() = default;

	~UsuallySmallMap() {
		if (m_value != nullptr) {
			if (m_size <= 0) {
				delete (ValueMapType*)m_value;
			} else if (m_size == 1) {
				delete (ValuePairType*)m_value;
			} else {
				delete (ValueArrayType*)m_value;
			}
		}
	}

	UsuallySmallMap(const UsuallySmallMap& copy) : m_size(copy.m_size) {
		if (copy.m_value == nullptr) {
			m_value = nullptr;
		} else if (m_size <= 0) {
			m_value = new ValueMapType(copy.getValueMap());
		} else if (m_size == 1) {
			m_value = new ValuePairType(copy.getValueScalar());
		} else {
			m_value = new ValueArrayType(copy.getValueArray());
		}
	}

	OGDF_COPY_MOVE_BY_SWAP(UsuallySmallMap)

	OGDF_SWAP_OP(UsuallySmallMap) {
		using std::swap;
		swap(first.m_size, second.m_size);
		swap(first.m_value, second.m_value);
	}

	void unset(const Key2& key) {
		if (!contains(key)) {
			return;
		}
		if (m_size <= 0) {
			ValueMapType& map = getValueMap();
			int cnt = map.erase(key);
			OGDF_ASSERT(cnt == 1);
		} else if (m_size == 1) {
			delete (ValuePairType*)m_value;
			m_value = nullptr;
			m_size = 0;
		} else {
			ValueArrayType& array = getValueArray();
			for (int i = 0; i < m_size; ++i) {
				if (array[i].first == key) {
					if (i != m_size - 1) {
						using std::swap;
						swap(array[i], array[m_size - 1]);
					}
					m_size--;
					break;
				}
			}
		}
	}

	Value& get_or_create(const Key2& key, const Value& def = Value()) {
		if (m_value == nullptr) {
			// create scalar value
			m_size = 1;
			ValuePairType* pair = new ValuePairType(key, def);
			m_value = pair;
			return pair->second;
		} else if (m_size <= 0) {
			ValueMapType& map = getValueMap();
			auto it = map.find(key);
			if (it != map.end()) {
				return it->second;
			}

			// insert into map
			auto ins = map.insert({key, def});
			OGDF_ASSERT(ins.second);
			return ins.first->second;
		} else if (m_size == 1) {
			ValuePairType& pair = getValueScalar();
			if (pair.first == key) {
				return pair.second;
			}

			// convert to array
			m_size = 2;
			ValueArrayType* array = new ValueArrayType {std::move(pair), ValuePairType(key, def)};
			delete (ValuePairType*)m_value;
			m_value = array;
			return (*array)[1].second;
		} else {
			ValueArrayType& array = getValueArray();
			for (int i = 0; i < m_size; ++i) {
				if (array[i].first == key) {
					return array[i].second;
				}
			}

			if (m_size < array_max) {
				// insert into array
				int i = m_size;
				m_size++;
				array[i] = ValuePairType(key, def);
				return array[i].second;
			} else {
				// convert to map
				ValueMapType* map = new ValueMapType();
				for (int i = 0; i < array_max; ++i) {
					map->insert(std::move(array[i]));
				}
				delete (ValueArrayType*)m_value;
				m_value = map;
				m_size = 0;
				return get_or_create(key);
			}
		}
	}

	Value& get_or_raise(const Key2& key) const {
		if (m_value == nullptr) {
			throw std::out_of_range("no keys stored");
		} else if (m_size <= 0) {
			ValueMapType& map = getValueMap();
			auto it = map.find(key);
			if (it == map.end()) {
				throw std::out_of_range("key not in map");
			}
			return it->second;
		} else if (m_size == 1) {
			ValuePairType& pair = getValueScalar();
			if (pair.first == key) {
				return pair.second;
			} else {
				throw std::out_of_range("key not in scalar");
			}
		} else {
			ValueArrayType& array = getValueArray();
			for (int i = 0; i < m_size; ++i) {
				if (array[i].first == key) {
					return array[i].second;
				}
			}
			throw std::out_of_range("key not in array");
		}
	}

	bool contains(const Key2& key) const {
		if (m_value == nullptr) {
			return false;
		} else if (m_size <= 0) {
			return getValueMap().count(key);
		} else if (m_size == 1) {
			return getValueScalar().first == key;
		} else {
			ValueArrayType& array = getValueArray();
			for (int i = 0; i < m_size; ++i) {
				if (array[i].first == key) {
					return true;
				}
			}
			return false;
		}
	}

	bool empty() const { return m_value == nullptr; }

	size_t size() const {
		if (m_value == nullptr) {
			return 0;
		} else if (m_size > 0) {
			return m_size;
		} else {
			return getValueMap().size();
		}
	}

	bool usesMap() const { return m_size == 0 && m_value != nullptr; }

	ValuePairType& getValueScalar() const {
		OGDF_ASSERT(m_size == 1);
		OGDF_ASSERT(m_value != nullptr);
		return *((ValuePairType*)m_value);
	}

	ValueArrayType& getValueArray() const {
		OGDF_ASSERT(m_size > 1);
		OGDF_ASSERT(m_value != nullptr);
		return *((ValueArrayType*)m_value);
	}

	ValueMapType& getValueMap() const {
		OGDF_ASSERT(m_size == 0);
		OGDF_ASSERT(m_value != nullptr);
		return *((ValueMapType*)m_value);
	}
};

template<typename Key1, typename Key2, typename Value, template<typename...> class BaseArray,
		int array_max = 64>
//! Data structure for two-dimensional mappings that are sparse in the second dimension.
/**
 * This is effectively a RegisteredArray that has two indexing dimensions, the first one uses a regular RegisteredArray,
 * while the second uses UsuallySmallMap as compact representation.
 * We assume that, for most values in the first dimension, there are only very few values in the second dimension.
 * Thus, for most values we only keep the single value in the second dimension or a short std::array through which we can easily seek.
 * If more than \p array_max values are stored for a first dimension value, we use a std::map to manage the second dimension for this value.
 */
class RegisteredMultiArray {
public:
	using EntryType = UsuallySmallMap<Key2, Value, array_max>;

	RegisteredMultiArray() = default;

	OGDF_DEFAULT_COPY(RegisteredMultiArray)
	OGDF_DEFAULT_MOVE(RegisteredMultiArray)

	template<class... T>
	explicit RegisteredMultiArray(T&&... t) : m_array(std::forward<T>(t)...) {};

	Value& operator()(const Key1& k1, const Key2& k2) { return m_array[k1].get_or_create(k2); }

	Value& get_or_create(const Key1& k1, const Key2& k2) { return m_array[k1].get_or_create(k2); }

	Value& get_or_raise(const Key1& k1, const Key2& k2) { return m_array[k1].get_or_raise(k2); }

	const Value& get_or_raise(const Key1& k1, const Key2& k2) const {
		return m_array[k1].get_or_raise(k2);
	}

	EntryType& get_all(const Key1& k1) { return m_array[k1]; }

	const EntryType& get_all(const Key1& k1) const { return m_array[k1]; }

	void remove(const Key1& k1, const Key2& k2) { return m_array[k1].unset(k2); }

	bool contains(const Key1& k1, const Key2& k2) const { return m_array[k1].contains(k2); }

	size_t count(const Key1& k1) const { return m_array[k1].size(k1); }

	bool has(const Key1& k1) const { return !m_array[k1].empty(k1); }

private:
	BaseArray<EntryType> m_array;
};

}
