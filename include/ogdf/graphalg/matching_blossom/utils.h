/** \file
 * \brief Utility functions and classes regarding map access and iteration.
 *
 * \author Joshua Sangmeister
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

#include <ogdf/basic/GraphAttributes.h>

#include <unordered_map>

namespace ogdf {
namespace matching_blossom {

//! Helper function to get the maximum value for a given weight type.
template<class TWeight>
TWeight infinity() {
	return std::numeric_limits<TWeight>::has_infinity ? std::numeric_limits<TWeight>::infinity()
													  : std::numeric_limits<TWeight>::max();
}

//! Return the pointer belonging to key \p key int the given map \p map, or nullptr if \p key does
//! not exist.
template<class K, class V>
V* tryGetPointerFromMap(const std::unordered_map<K, V*>& map, const K& key) {
	auto it = map.find(key);
	if (it != map.end()) {
		return it->second;
	} else {
		return nullptr;
	}
}

//! Helper function to get the edge weight of \p e from the EdgeArray \p weights.
template<class TWeight>
TWeight getWeight(edge e, const EdgeArray<TWeight>& weights) {
	return weights[e];
}

//! Helper function to get the edge weight of \p e from the GraphAttributes \p weights.
template<class TWeight>
TWeight getWeight(edge e, const GraphAttributes& GA) {
	if (std::numeric_limits<TWeight>::is_integer) {
		if (GA.has(GraphAttributes::edgeIntWeight)) {
			return GA.intWeight(e);
		}
	} else if (GA.has(GraphAttributes::edgeIntWeight)) {
		return GA.doubleWeight(e);
	}
	return 1;
}

//! Iterator to access the keys of a std::unordered_map.
template<typename Key, typename Value>
class MapKeyIterator : public std::unordered_map<Key, Value>::iterator {
	using MapIterator = typename std::unordered_map<Key, Value>::iterator;

public:
	MapKeyIterator() : MapIterator() {};
	MapKeyIterator(MapIterator it_) : MapIterator(it_) {};

	Key* operator->() { return (Key* const)&(MapIterator::operator->()->first); }

	Key operator*() { return MapIterator::operator*().first; }
};

//! Iterator to access the values of a std::unordered_map.
template<typename Key, typename Value>
class MapValueIterator : public std::unordered_map<Key, Value>::iterator {
	using MapIterator = typename std::unordered_map<Key, Value>::iterator;

public:
	MapValueIterator() : MapIterator() {};
	MapValueIterator(MapIterator it_) : MapIterator(it_) {};

	Value* operator->() { return (Value* const)&(MapIterator::operator->()->second); }

	Value operator*() { return MapIterator::operator*().second; }
};

//! Dummy class for scoped iteration of a std::unordered_map.
template<template<typename, typename> class Iterator, typename Key, typename Value>
class BaseIteratorContainer {
	using iterator = Iterator<Key, Value>;

	std::unordered_map<Key, Value>& m_map;

public:
	BaseIteratorContainer(std::unordered_map<Key, Value>& map) : m_map(map) { }

	iterator begin() { return iterator(m_map.begin()); }

	iterator end() { return iterator(m_map.end()); }

	size_t size() { return m_map.size(); }
};

template<typename Key, typename Value>
using KeyIteratorContainer = BaseIteratorContainer<MapKeyIterator, Key, Value>;

template<typename Key, typename Value>
using ValueIteratorContainer = BaseIteratorContainer<MapValueIterator, Key, Value>;

}
}
