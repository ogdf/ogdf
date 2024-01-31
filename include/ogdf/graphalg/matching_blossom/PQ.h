/** \file
 * \brief An extension of the standard priority queue with additional
 * features.
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

#include <ogdf/basic/PriorityQueue.h>
#include <ogdf/graphalg/matching_blossom/utils.h>

#include <limits>
#include <unordered_map>

namespace ogdf {
namespace matching_blossom {

//! A custom priority queue for the blossom algorithm, based on the PrioritizedMapQueue. It uses a
//! std::unordered_map to store the handles, handles merge operations correctly, can iterate over its elements
//! and offers a remove() function.
template<typename E, typename TWeight, typename C = std::less<TWeight>,
		template<typename, class> class Impl = PairingHeap>
class BlossomPQ
	: public PrioritizedQueue<E, TWeight, C, Impl>,
	  public KeyIteratorContainer<E, typename PrioritizedQueue<E, TWeight, C, Impl>::Handle> {
protected:
	using ThisQueue = BlossomPQ<E, TWeight, C, Impl>;
	using SuperQueue = PrioritizedQueue<E, TWeight, C, Impl>;
	using Handle = typename SuperQueue::Handle;

	std::unordered_map<E, Handle> m_handles;

public:
	BlossomPQ() : KeyIteratorContainer<E, Handle>(m_handles) { }

	//! Returns whether this queue contains that key
	bool contains(const E& element) const { return m_handles.find(element) != m_handles.end(); }

	/*
	 * Returns the priority of the key.
	 * Note that the behaviour is undefined if the key is not present.
	 */
	TWeight priority(const E& element) const {
		auto it = m_handles.find(element);
		OGDF_ASSERT(it != m_handles.end());
		return this->value(it->second).priority();
	}

	/**
	 * Adds a new element to the queue.
	 * Note that the behaviour is undefined if the key is already present.
	 */
	void push(const E& element, const TWeight priority) {
		OGDF_ASSERT(!contains(element));
		m_handles[element] = SuperQueue::push(element, priority);
	}

	//! Removes the topmost element from the queue.
	void pop() {
		m_handles.erase(SuperQueue::topElement());
		SuperQueue::pop();
	}

	/**
	 * Decreases the priority of the given element.
	 * Note that the behaviour is undefined if the key is not present
	 * or the given priority is greater than the former one.
	 */
	void decrease(const E& element, const TWeight priority) {
		Handle pos = m_handles[element];
		SuperQueue::decrease(pos, priority);
	}

	//! Merge this Priority Queue with another one. The other queue is cleared afterwards.
	void merge(ThisQueue& other) {
		SuperQueue::merge(other);
		for (auto entry : other.m_handles) {
			m_handles[entry.first] = entry.second;
		}
		other.m_handles.clear();
	}

	//! Removes all elements from this queue.
	void clear() {
		SuperQueue::clear();
		m_handles.clear();
	}

	//! Remove \p e from the queue by decreasing its priority to the minimum. \p e must be present
	//! in the queue; if that's not guaranteed, use tryRemove() instead.
	void remove(const E& e) {
		OGDF_ASSERT(this->contains(e));
		this->decrease(e, -infinity<TWeight>());
		OGDF_ASSERT(this->topElement() == e);
		this->pop();
	}
};

}
}
