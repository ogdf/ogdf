/** \file
 * \brief Implementation of randomized meldable heap data structure.
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

#ifdef _MSC_VER
#pragma once
#endif


#ifndef OGDF_IMPL_RM_HEAP_H
#define OGDF_IMPL_RM_HEAP_H


#include <functional>
#include <random>

#include <ogdf/internal/heap/HeapBase.h>


namespace ogdf {


//! Randomized meldable heap node.
template<typename T>
struct RMHeapNode {
	template<typename, typename> friend class RMHeap;
protected:
	T value; //!< Value contained in the node.

	RMHeapNode<T> *parent; //!< Parent of the node.
	RMHeapNode<T> *left; //!< Left child of the node.
	RMHeapNode<T> *right; //!< Right child of the node.

	//! Creates heap node with a given \a value.
	RMHeapNode(const T &value)
	: value(value),
	  parent(nullptr), left(nullptr), right(nullptr)
	{
	}
};


//! Randomized meldable heap implementation.
/**
 * @ingroup containers
 *
 * Code of meld (also known as merge) operation is solely based on Wikipedia
 * article. Other methods are based on my intuitions and make use of melding.
 *
 * For random number generation it uses default generator provided by the C++11
 * standard. In the future, it should be possible to provide custom random
 * device, generator and seed.
 *
 * @tparam T Denotes value type of inserted elements.
 * @tparam C Denotes comparison functor determining value ordering.
 */
template<typename T, typename C=std::less<T>>
class RMHeap : public HeapBase<RMHeap, RMHeapNode<T>, T, C> {
public:

	/**
	 * Creates empty randomized meldable heap.
	 *
	 * @param cmp Comparison functor determining value ordering.
	 * @param initialSize ignored by this implementation.
	 */
	RMHeap(const C &cmp = C(), int initialSize = -1);

	/**
	 * Destructs the heap.
	 *
	 * If the heap is not empty, destructors of contained elements are called
	 * and used storage is deallocated.
	 */
	virtual ~RMHeap();

	//! Returns reference to the top element in the heap.
	const T &top() const;

	/**
	 * Inserts a new node with given \a value into a heap.
	 *
	 * @param value A value to be inserted.
	 * @return Handle to the inserted node.
	 */
	RMHeapNode<T> *push(const T &value);

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
	void decrease(RMHeapNode<T> *node, const T &value);

	/**
	 * Merges in values of \a other heap.
	 *
	 * After merge \a other heap becomes empty and is valid for further usage.
	 *
	 * @param other A heap to be merged in.
	 */
	void merge(RMHeap<T, C> &other);

	/*
	* Retuns the value of the node
	*
	* @param node The nodes handle
	* @return the value of the node
	*/
	const T &value(RMHeapNode<T> *node) const {
		return node->value;
	}

private:
	std::default_random_engine m_rand; //!< Random values generator.
	RMHeapNode<T> *m_root; //!< Root node of the heap.

	//! Recursively merges heaps \a a and \a b. Returns resulting heap.
	RMHeapNode<T> *merge(RMHeapNode<T> *a, RMHeapNode<T> *b);
	//! Removes given \a node from the main tree (does not free memory!).
	void remove(RMHeapNode<T> *node);

	//! Recursively releases memory occupied by heap pointed by \a node.
	static void release(RMHeapNode<T> *node);
};


template<typename T, typename C>
RMHeap<T, C>::RMHeap(const C &cmp, int initialSize)
: m_rand((std::random_device())()), m_root(nullptr)
{
	this->m_comp = cmp;
}


template<typename T, typename C>
RMHeap<T, C>::~RMHeap()
{
	release(m_root);
}


template<typename T, typename C>
const T &RMHeap<T, C>::top() const
{
	return m_root->value;
}


template<typename T, typename C>
RMHeapNode<T> *RMHeap<T, C>::push(const T &value)
{
	RMHeapNode<T> *node = new RMHeapNode<T>(value);
	m_root = merge(m_root, node);

	return node;
}


template<typename T, typename C>
void RMHeap<T, C>::pop()
{
	RMHeapNode<T> *root = m_root;
	m_root = merge(m_root->left, m_root->right);
	if(m_root != nullptr) {
		m_root->parent = nullptr;
	}
	delete root;
}


template<typename T, typename C>
void RMHeap<T, C>::decrease(RMHeapNode<T> *node, const T &value)
{
	node->value = value;
	if(node == m_root) {
		return;
	}

	remove(node);
	node->left = nullptr;
	node->right = nullptr;
	node->parent = nullptr;

	m_root = merge(m_root, node);
}


template<typename T, typename C>
void RMHeap<T, C>::merge(RMHeap<T, C> &other)
{
	m_root = merge(m_root, other.m_root);
	other.m_root = nullptr;
}


template<typename T, typename C>
RMHeapNode<T> *RMHeap<T, C>::merge(RMHeapNode<T> *a, RMHeapNode<T> *b)
{
	if(a == nullptr) {
		return b;
	}
	if(b == nullptr) {
		return a;
	}

	if(this->m_comp(a->value, b->value)) {
		if(m_rand() % 2 == 0) {
			a->left = merge(a->left, b);
			if(a->left != nullptr) {
				a->left->parent = a;
			}
		} else {
			a->right = merge(a->right, b);
			if(a->right != nullptr) {
				a->right->parent = a;
			}
		}
		return a;
	} else {
		if(m_rand() % 2 == 0) {
			b->left = merge(b->left, a);
			if(b->left != nullptr) {
				b->left->parent = b;
			}
		} else {
			b->right = merge(b->right, a);
			if(b->right != nullptr) {
				b->right->parent = b;
			}
		}
		return b;
	}
}


template<typename T, typename C>
void RMHeap<T, C>::remove(RMHeapNode<T> *node)
{
	RMHeapNode<T> *merged = merge(node->left, node->right);
	if(node == node->parent->left) {
		node->parent->left = merged;
	} else {
		node->parent->right = merged;
	}
	if(merged != nullptr) {
		merged->parent = node->parent;
	}
}


template<typename T, typename C>
void RMHeap<T, C>::release(RMHeapNode<T> *node)
{
	if(node == nullptr) {
		return;
	}

	release(node->left);
	release(node->right);
	delete node;
}

} // end namespace ogdf


#endif
