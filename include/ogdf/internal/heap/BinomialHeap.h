/** \file
 * \brief Implementation of binomial heap data structure.
 *
 * \author Łukasz Hanuszczak
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.txt in the root directory of the OGDF installation for details.
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
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/

#pragma once

#include <functional>
#include <utility>

#include <ogdf/internal/heap/HeapBase.h>


namespace ogdf {


//! Binomial heap node.
template<typename T>
struct BinomialHeapNode {
	template<typename, typename> friend class BinomialHeap;
protected:
	T value; //!< Value contained in the node.

	size_t rank; //!< Determines rank of a node.

	BinomialHeapNode<T> *parent; //!< Parent of the node.
	BinomialHeapNode<T> *next; //!< Next sibling of the node.
	BinomialHeapNode<T> *child; //!< First child of the node.

	//! Creates heap node with a given \a value.
	BinomialHeapNode(const T &value)
	: value(value),
	  rank(0), parent(nullptr), next(nullptr), child(nullptr)
	{
	}
};


//! Binomial heap implementation.
/**
 * @ingroup containers
 *
 * Code is mainly based on samples and ideas provided in "Introduction to
 * Algorithms" book (aka "Cormen").
 *
 * @tparam T Denotes value type of inserted elements.
 * @tparam C Denotes comparison functor determining value ordering.
 */
template<typename T, typename C = std::less<T>>
class BinomialHeap : public HeapBase<BinomialHeap, BinomialHeapNode<T>, T, C> {
public:

	/**
	 * Creates empty binomial heap.
	 *
	 * @param cmp Comparison functor determining value ordering.
	 * @param initialSize ignored by this implementation.
	 */
	BinomialHeap(const C &cmp = C(), int initialSize = -1);

	/**
	 * Destructs the heap.
	 *
	 * If the heap is not empty, destructors of contained elements are called
	 * and used storage is deallocated.
	 */
	virtual ~BinomialHeap();

	//! Returns reference to the top element in the heap.
	const T &top() const;

	/**
	 * Inserts a new node with given \a value into a heap.
	 *
	 * @param value A value to be inserted.
	 * @return Handle to the inserted node.
	 */
	BinomialHeapNode<T> *push(const T &value);

	/**
	 * Removes the top element from the heap.
	 *
	 * Behaviour of this function is undefined if the heap is empty.
	 */
	void pop();

	/**
	 * Decreases value of the given \a node to \a value.
	 *
	 * Behaviour of this function is undefined if node does not belong to a the
	 * heap or new value is greater than old one.
	 *
	 * @param node A node for which the value is to be decreased.
	 * @param value A new value for the node.
	 */
	void decrease(BinomialHeapNode<T> *node, const T &value);

	/**
	 * Merges in values of \a other heap.
	 *
	 * After merge \a other heap becomes empty and is valid for further usage.
	 *
	 * @param other A heap to be merged in.
	 */
	void merge(BinomialHeap<T, C> &other);

	/*
	* Retuns the value of the node
	*
	* @param node The nodes handle
	* @return the value of the node
	*/
	const T &value(BinomialHeapNode<T> *node) const {
		return node->value;
	}

private:
	BinomialHeapNode<T> *m_root; //!< Root node of the heap.

	//! Joins heap lists \a a and \a b into single list sorted by the ranks.
	BinomialHeapNode<T> *join(BinomialHeapNode<T> *a, BinomialHeapNode<T> *b);
	//! Merges in \a other heap list into the heap.
	void merge(BinomialHeapNode<T> *other);

	//! Makes \a child node a child of \a parent node.
	static void link(BinomialHeapNode<T> *parent, BinomialHeapNode<T> *child);

	//! Releases memory occupied by list of heaps given as \a node.
	static void release(BinomialHeapNode<T> *node);
};


template<typename T, typename C>
BinomialHeap<T, C>::BinomialHeap(const C &cmp, int initialSize) : m_root(nullptr)
{
	this->m_comp = cmp;
}


template<typename T, typename C>
BinomialHeap<T, C>::~BinomialHeap()
{
	release(m_root);
	m_root = nullptr;
}


template<typename T, typename C>
void BinomialHeap<T, C>::release(BinomialHeapNode<T> *node)
{
	while(node != nullptr) {
		release(node->child);

		BinomialHeapNode<T> *next = node->next;
		delete node;
		node = next;
	}
}


template<typename T, typename C>
inline const T &BinomialHeap<T, C>::top() const
{
	BinomialHeapNode<T> *min = m_root;
	for(BinomialHeapNode<T> *it = m_root->next; it != nullptr; it = it->next) {
		if(this->m_comp(it->value, min->value)) {
			min = it;
		}
	}

	return min->value;
}


template<typename T, typename C>
BinomialHeapNode<T> *BinomialHeap<T, C>::push(const T &value)
{
	BinomialHeapNode<T> *node = new BinomialHeapNode<T>(value);

	merge(node);
	return node;
}


template<typename T, typename C>
void BinomialHeap<T, C>::pop()
{
	BinomialHeapNode<T> *curr = m_root, *min = m_root, *minPrev = nullptr;

	while(curr->next != nullptr) {
		if(this->m_comp(curr->next->value, min->value)) {
			min = curr->next;
			minPrev = curr;
		}
		curr = curr->next;
	}

	if(min == m_root) {
		m_root = min->next;
	} else {
		minPrev->next = min->next;
	}

	// Children list has to be reversed before it can be merged.
	BinomialHeapNode<T> *reversed = nullptr, *child = min->child;
	while(child != nullptr) {
		BinomialHeapNode<T> *next = child->next;
		child->next = reversed;
		reversed = child;
		child = next;
	}
	merge(reversed);
	delete min;
}


template<typename T, typename C>
void BinomialHeap<T, C>::decrease(BinomialHeapNode<T> *node, const T &value)
{
	// BinomialHeap::decrease is not supported
	OGDF_ASSERT(false);

	node->value = value;

	while(node->parent != nullptr &&
	      this->m_comp(node->value, node->parent->value))
	{
		std::swap(node->value, node->parent->value);
		node = node->parent;
	}
}


template<typename T, typename C>
void BinomialHeap<T, C>::merge(BinomialHeap<T, C> &other)
{
	merge(other.m_root);
	other.m_root = nullptr;
}


template<typename T, typename C>
inline BinomialHeapNode<T> *BinomialHeap<T, C>::join(
	BinomialHeapNode<T> *a, BinomialHeapNode<T> *b)
{
	if(a == nullptr) {
		return b;
	}
	if(b == nullptr) {
		return a;
	}

	if(b->rank < a->rank) {
		std::swap(a, b);
	}

	BinomialHeapNode<T> *head = a;
	while(b != nullptr) {
		if(a->next == nullptr) {
			a->next = b;
			break;
		}

		if(b->rank < a->next->rank) {
			BinomialHeapNode<T> *nextB = b->next;
			b->next = a->next;
			a->next = b;

			a = b;
			b = nextB;
		} else {
			a = a->next;
		}
	}

	return head;
}


template<typename T, typename C>
inline void BinomialHeap<T, C>::merge(BinomialHeapNode<T> *other)
{
	m_root = join(m_root, other);
	if(m_root == nullptr) {
		return;
	}

	BinomialHeapNode<T> *prev = nullptr, *curr = m_root, *next = m_root->next;
	while(next != nullptr) {
		if(curr->rank != next->rank || (next->next != nullptr &&
		                                next->next->rank == curr->rank))
		{
			prev = curr;
			curr = next;
			next = curr->next;
			continue;
		}

		if(this->m_comp(curr->value, next->value)) {
			curr->next = next->next;
			link(curr, next);
		} else if(prev == nullptr) {
			m_root = next;
			link(next, curr);
			curr = next;
		} else {
			prev->next = next;
			link(next, curr);
			curr = next;
		}
		next = curr->next;
	}
}


template<typename T, typename C>
inline void BinomialHeap<T, C>::link(
	BinomialHeapNode<T> *parent, BinomialHeapNode<T> *child)
{
	child->next = parent->child;
	child->parent = parent;
	parent->child = child;
	parent->rank++;
}

} // end namespace ogdf
