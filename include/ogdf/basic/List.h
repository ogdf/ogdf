/** \file
 * \brief Declaration of doubly linked lists and iterators
 *
 * \author Carsten Gutwenger and Sebastian Leipert
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

#ifndef OGDF_LIST_H
#define OGDF_LIST_H


#include <ogdf/internal/basic/list_templates.h>
#include <random>


namespace ogdf {


template<class E> class List;
template<class E> class ListPure;
template<class E> class ListIterator;
template<class E> class ListConstIterator;


//! The parameterized class \a ListElement<E> represents the structure for elements of doubly linked lists.
template<class E>
class ListElement {
	friend class ListPure<E>;
	friend class List<E>;
	friend class ListIterator<E>;
	friend class ListConstIterator<E>;

	ListElement<E> *m_next; //!< Pointer to successor element.
	ListElement<E> *m_prev; //!< Pointer to predecessor element.
	E m_x; //!< Stores the content.

	//! Constructs a ListElement.
	ListElement() : m_next(nullptr), m_prev(nullptr) { }
	//! Constructs a ListElement.
	ListElement(const E &x) : m_next(nullptr), m_prev(nullptr), m_x(x) { }
	//! Constructs a ListElement.
	ListElement(const E &x, ListElement<E> *next, ListElement<E> *prev) :
		m_next(next), m_prev(prev), m_x(x) { }
	//! Constructs a ListElement with given arguments \a args for constructor of element type.
	template<class ... Args>
	ListElement(ListElement<E> *next, ListElement<E> *prev, Args && ... args)
		: m_next(next), m_prev(prev), m_x(std::forward<Args>(args)...) { }

	OGDF_NEW_DELETE
}; // class ListElement



//! The parameterized class \a ListIterator<E> encapsulates a pointer to a dlist element.
/**
 * It is used in order to iterate over doubly linked lists,
 * and to specify a position in a doubly linked list. It is possible that
 * an iterator encapsulates a null pointer.
 */

template<class E> class ListIterator {
	ListElement<E> *m_pX; // pointer to associated list element

	friend class ListConstIterator<E>;
	friend class ListPure<E>;

	//! Conversion to pointer to list element.
	operator ListElement<E> *() { return m_pX; }
	//! Conversion to pointer to list element.
	operator const ListElement<E> *() const { return m_pX; }

public:
	//! Constructs an iterator pointing to no element.
	ListIterator() : m_pX(nullptr) { }
	//! Constructs an iterator pointing to \a pX.
	ListIterator(ListElement<E> *pX) : m_pX(pX) { }
	//! Constructs an iterator that is a copy of \a it.
	ListIterator(const ListIterator<E> &it) : m_pX(it.m_pX) { }

	//! Returns true iff the iterator points to an element.
	bool valid() const { return m_pX != nullptr; }

	//! Equality operator.
	bool operator==(const ListIterator<E> &it) const {
		return m_pX == it.m_pX;
	}

	//! Inequality operator.
	bool operator!=(const ListIterator<E> &it) const {
		return m_pX != it.m_pX;
	}

	//! Returns successor iterator.
	ListIterator<E> succ() const { return m_pX->m_next; }

	//! Returns predecessor iterator.
	ListIterator<E> pred() const { return m_pX->m_prev; }

	//! Returns a reference to the element content.
	E &operator*() const { return m_pX->m_x; }

	//! Assignment operator.
	ListIterator<E> &operator=(const ListIterator<E> &it) {
		m_pX = it.m_pX;
		return *this;
	}

	//! Increment operator (prefix).
	ListIterator<E> &operator++() {
		m_pX = m_pX->m_next;
		return *this;
	}

	//! Increment operator (postfix).
	ListIterator<E> operator++(int) {
		ListIterator<E> it = *this;
		m_pX = m_pX->m_next;
		return it;
	}

	//! Decrement operator (prefix).
	ListIterator<E> &operator--() {
		m_pX = m_pX->m_prev;
		return *this;
	}

	//! Decrement operator (postfix).
	ListIterator<E> operator--(int) {
		ListIterator<E> it = *this;
		m_pX = m_pX->m_prev;
		return it;
	}

	OGDF_NEW_DELETE
}; // class ListIterator



//---------------------------------------------------------
// ListConstIterator<E>
// const iterator for doubly linked lists
//---------------------------------------------------------
//! The parameterized class \a ListIterator<E> encapsulates a constant pointer to a list element.
/**
 * It is used in order to iterate over doubly linked lists,
 * and to specify a position in a doubly linked list. It is possible that
 * an iterator encapsulates a null pointer. In contrast to ListIterator,
 * it is not possible to change the list element pointed to.
 */

template<class E> class ListConstIterator {
	const ListElement<E> *m_pX; // pointer to list element

	friend class ListPure<E>;

	//! Conversion to pointer to list element.
	operator const ListElement<E> *() { return m_pX; }

public:
	//! Constructs an iterator pointing to no element.
	ListConstIterator() : m_pX(nullptr) { }

	//! Constructs an iterator pointing to \a pX.
	ListConstIterator(const ListElement<E> *pX) : m_pX(pX) { }

	//! Constructs an iterator that is a copy of \a it.
	ListConstIterator(const ListIterator<E> &it) : m_pX((const ListElement<E> *)it) { }
	//! Constructs an iterator that is a copy of \a it.
	ListConstIterator(const ListConstIterator &it) : m_pX(it.m_pX) { }

	//! Returns true iff the iterator points to an element.
	bool valid() const { return m_pX != nullptr; }

	//! Equality operator.
	bool operator==(const ListConstIterator<E> &it) const {
		return m_pX == it.m_pX;
	}

	//! Inequality operator.
	bool operator!=(const ListConstIterator<E> &it) const {
		return m_pX != it.m_pX;
	}

	//! Returns successor iterator.
	ListConstIterator<E> succ() const { return m_pX->m_next; }

	//! Returns predecessor iterator.
	ListConstIterator<E> pred() const { return m_pX->m_prev; }

	//! Returns a reference to the element content.
	const E &operator*() const { return m_pX->m_x; }

	//! Assignment operator.
	ListConstIterator<E> &operator=(const ListConstIterator<E> &it) {
		m_pX = it.m_pX;
		return *this;
	}

	//! Increment operator (prefix).
	ListConstIterator<E> &operator++() {
		m_pX = m_pX->m_next;
		return *this;
	}

	//! Increment operator (postfix).
	ListConstIterator<E> operator++(int) {
		ListConstIterator<E> it = *this;
		m_pX = m_pX->m_next;
		return it;
	}

	//! Decrement operator (prefix).
	ListConstIterator<E> &operator--() {
		m_pX = m_pX->m_prev;
		return *this;
	}

	//! Decrement operator (postfix).
	ListConstIterator<E> operator--(int) {
		ListConstIterator<E> it = *this;
		m_pX = m_pX->m_prev;
		return it;
	}

	OGDF_NEW_DELETE
}; // class ListConstIterator



//! Doubly linked lists.
/**
 * @ingroup containers
 *
 * Elements of the list are instances of type ListElement.
 * Use ListConstIterator or ListIterator in order to iterate over the list.
 *
 * In contrast to List, instances of ListPure do not store the length of the list.
 *
 * @tparam E is the data type stored in list elements.
 */

template<class E> class ListPure {
protected:

	ListElement<E> *m_head; //!< Pointer to first element.
	ListElement<E> *m_tail; //!< Pointer to last element.

public:
	//! Represents the data type stored in a list element.
	typedef E value_type;
	//! Provides a reference to an element stored in a list.
	typedef E &reference;
	//! Provides a reference to a const element stored in a list for reading and performing const operations.
	typedef const E &const_reference;
	//! Provides a bidirectional iterator that can read a const element in a list.
	typedef ListConstIterator<E> const_iterator;
	//! Provides a bidirectional iterator that can read or modify any element in a list.
	typedef ListIterator<E> iterator;

	//! Constructs an empty doubly linked list.
	ListPure() : m_head(nullptr), m_tail(nullptr) { }

	//! Constructs a doubly linked list containing the elements in \a init.
	ListPure(std::initializer_list<E> init) : m_head(nullptr), m_tail(nullptr) {
		for (const E &x : init)
			pushBack(x);
	}

	//! Constructs a doubly linked list that is a copy of \a L.
	ListPure(const ListPure<E> &L) : m_head(nullptr), m_tail(nullptr) {
		copy(L);
	}

	//! Constructs a doubly linked list containing the elements of \a L (move semantics).
	/**
	 * The list \a L is empty afterwards.
	 */
	ListPure(ListPure<E> &&L) : m_head(L.m_head), m_tail(L.m_tail) {
		L.m_head = L.m_tail = nullptr;
	}

	//! Destructor.
	~ListPure() { clear(); }


	/**
	 * @name Access methods
	 * These methods provide simple access without changing the list.
	 */
	//@{

	//! Returns true iff the list is empty.
	bool empty() const { return m_head == nullptr; }

	//! Returns the number of elements in the list.
	/**
	 * Notice that this method requires to iterate over the whole list and takes linear running time!
	 * If you require frequent access to the size of the list, consider using List instead of ListPure.
	 */
	int size() const {
		int count = 0;
		for (ListElement<E> *pX = m_head; pX; pX = pX->m_next)
			++count;
		return count;
	}

	//! Returns a const reference to the first element.
	/**
	 * \pre The list is not empty!
	 */
	const_reference front() const {
		OGDF_ASSERT(m_head != nullptr)
		return m_head->m_x;
	}

	//! Returns a reference to the first element.
	/**
	 * \pre The list is not empty!
	 */
	reference front() {
		OGDF_ASSERT(m_head != nullptr)
		return m_head->m_x;
	}

	//! Returns a const reference to the last element.
	/**
	 * \pre The list is not empty!
	 */
	const_reference back() const {
		OGDF_ASSERT(m_tail != nullptr)
		return m_tail->m_x;
	}

	//! Returns a reference to the last element.
	/**
	 * \pre The list is not empty!
	 */
	reference back() {
		OGDF_ASSERT(m_tail != nullptr)
		return m_tail->m_x;
	}

	//! Returns a const iterator pointing to the element at position \a pos.
	/**
	 * The running time of this method is linear in \a pos.
	 */
	const_iterator get(int pos) const {
		ListElement<E> *pX;
		for(pX = m_head; pX != nullptr; pX = pX->m_next)
			if (pos-- == 0) break;
		return pX;
	}

	//! Returns an iterator pointing to the element at position \a pos.
	/**
	 * The running time of this method is linear in \a pos.
	 */
	iterator get(int pos) {
		ListElement<E> *pX;
		for(pX = m_head; pX != nullptr; pX = pX->m_next)
			if (pos-- == 0) break;
		return pX;
	}

	//! Returns the position (starting with 0) of iterator \a it in the list.
	/**
	 * \pre \a it is a valid iterator pointing to an element in this list!
	 */
	int pos(const_iterator it) const {
		OGDF_ASSERT(it.valid())
		int p = 0;
		for(ListElement<E> *pX = m_head; pX != nullptr; pX = pX->m_next, ++p)
			if (pX == it) break;
		return p;
	}


	//@}
	/**
	 * @name Iterators
	 * These methods return bidirectional iterators to elements in the list and allow to iterate over the elements in linear or cyclic order.
	 */
	//@{

	//! Returns an iterator to the first element of the list.
	/**
	 * If the list is empty, a null pointer iterator is returned.
	 */
	iterator begin() { return m_head; }

	//! Returns a const iterator to the first element of the list.
	/**
	 * If the list is empty, a null pointer iterator is returned.
	 */
	const_iterator begin() const { return m_head; }

	//! Returns a const iterator to the first element of the list.
	/**
	 * If the list is empty, a null pointer iterator is returned.
	 */
	const_iterator cbegin() const { return m_head; }

	//! Returns an iterator to one-past-last element of the list.
	/**
	 * This is always a null pointer iterator.
	 */
	iterator end() { return iterator(); }

	//! Returns a const iterator to one-past-last element of the list.
	/**
	 * This is always a null pointer iterator.
	 */
	const_iterator end() const { return const_iterator(); }

	//! Returns a const iterator to one-past-last element of the list.
	/**
	 * This is always a null pointer iterator.
	 */
	const_iterator cend() const { return const_iterator(); }

	//! Returns an iterator to the last element of the list.
	/**
	 * If the list is empty, a null pointer iterator is returned.
	 */
	iterator rbegin() { return m_tail; }

	//! Returns a const iterator to the last element of the list.
	/**
	 * If the list is empty, a null pointer iterator is returned.
	 */
	const_iterator rbegin() const { return m_tail; }

	//! Returns a const iterator to the last element of the list.
	/**
	 * If the list is empty, a null pointer iterator is returned.
	 */
	const_iterator crbegin() const { return m_tail; }

	//! Returns an iterator to one-before-first element of the list.
	/**
	 * This is always a null pointer iterator.
	 */
	iterator rend() { return iterator(); }

	//! Returns a const iterator to one-before-first element of the list.
	/**
	 * This is always a null pointer iterator.
	 */
	const_iterator rend() const { return const_iterator(); }

	//! Returns a const iterator to one-before-first element of the list.
	/**
	 * This is always a null pointer iterator.
	 */
	const_iterator crend() const { return const_iterator(); }

	//! Returns a const iterator to the cyclic successor of \a it.
	/**
	 * \pre \a it points to an element in this list or to nullptr!
	 */
	const_iterator cyclicSucc(const_iterator it) const {
		const ListElement<E> *pX = it;
		return (pX && pX->m_next) ? pX->m_next : m_head;
	}

	//! Returns an iterator to the cyclic successor of \a it.
	/**
	 * \pre \a it points to an element in this list or to nullptr!
	 */
	iterator cyclicSucc(iterator it) {
		ListElement<E> *pX = it;
		return (pX && pX->m_next) ? pX->m_next : m_head;
	}

	//! Returns a const iterator to the cyclic predecessor of \a it.
	/**
	 * \pre \a it points to an element in this list or to nullptr!
	 */
	const_iterator cyclicPred(const_iterator it) const {
		const ListElement<E> *pX = it;
		return (pX && pX->m_prev) ? pX->m_prev : m_tail;
	}

	//! Returns an iterator to the cyclic predecessor of \a it.
	/**
	 * \pre \a it points to an element in this list or to nullptr!
	 */
	iterator cyclicPred(iterator it) {
		ListElement<E> *pX = it;
		return (pX && pX->m_prev) ? pX->m_prev : m_tail;
	}


	//@}
	/**
	 * @name Operators
	 * The following operators are provided by lists.
	 */
	//@{

	//! Assignment operator.
	ListPure<E> &operator=(const ListPure<E> &L) {
		clear(); copy(L);
		return *this;
	}

	//! Assignment operator (move semantics).
	/**
	 * The list \a L is empty afterwards.
	 */
	ListPure<E> &operator=(ListPure<E> &&L) {
		clear();
		m_head = L.m_head;
		m_tail = L.m_tail;
		L.m_head = L.m_tail = nullptr;
		return *this;
	}

	//! Equality operator.
	bool operator==(const ListPure<E> &L) const {
		ListElement<E> *pX = m_head, *pY = L.m_head;
		while(pX != nullptr && pY != nullptr) {
			if(pX->m_x != pY->m_x)
				return false;
			pX = pX->m_next;
			pY = pY->m_next;
		}
		return (pX == nullptr && pY == nullptr);
	}

	//! Inequality operator.
	bool operator!=(const ListPure<E> &L) const {
		return !operator==(L);
	}


	//@}
	/**
	 * @name Adding elements
	 * These method add elements to the list.
	 */
	//@{

	//! Adds element \a x at the beginning of the list.
	iterator pushFront(const E &x) {
		ListElement<E> *pX = OGDF_NEW ListElement<E>(x,m_head,nullptr);
		if (m_head)
			m_head = m_head->m_prev = pX;
		else
			m_head = m_tail = pX;
		return m_head;
	}

	//! Adds a new element at the beginning of the list.
	/**
	 * The element is constructed in-place with exactly the same arguments \a args as supplied to the function.
	 */
	template<class ... Args>
	iterator emplaceFront(Args && ... args) {
		ListElement<E> *pX = OGDF_NEW ListElement<E>(m_head, nullptr, std::forward<Args>(args)...);
		if (m_head)
			m_head = m_head->m_prev = pX;
		else
			m_head = m_tail = pX;
		return m_head;
	}

	//! Adds element \a x at the end of the list.
	iterator pushBack(const E &x) {
		ListElement<E> *pX = OGDF_NEW ListElement<E>(x,nullptr,m_tail);
		if (m_head)
			m_tail = m_tail->m_next = pX;
		else
			m_tail = m_head = pX;
		return m_tail;
	}

	//! Adds a new element at the end of the list.
	/**
	* The element is constructed in-place with exactly the same arguments \a args as supplied to the function.
	*/
	template<class ... Args>
	iterator emplaceBack(Args && ... args) {
		ListElement<E> *pX = OGDF_NEW ListElement<E>(nullptr, m_tail, std::forward<Args>(args)...);
		if (m_head)
			m_tail = m_tail->m_next = pX;
		else
			m_tail = m_head = pX;
		return m_tail;
	}

	//! Inserts element \a x before or after \a it.
	/**
	 * @param x is the element to be inserted.
	 * @param it is a list iterator in this list.
	 * @param dir determines if \a x is inserted before or after \a it.
	 *   Possible values are \c ogdf::before and \c ogdf::after.
	 * \pre \a it points to an element in this list.
	 */
	iterator insert(const E &x, iterator it, Direction dir = after) {
		OGDF_ASSERT(it.valid())
		OGDF_ASSERT(dir == after || dir == before)
		ListElement<E> *pY = it, *pX;
		if (dir == after) {
			ListElement<E> *pYnext = pY->m_next;
			pY->m_next = pX = OGDF_NEW ListElement<E>(x,pYnext,pY);
			if (pYnext) pYnext->m_prev = pX;
			else m_tail = pX;
		} else {
			ListElement<E> *pYprev = pY->m_prev;
			pY->m_prev = pX = OGDF_NEW ListElement<E>(x,pY,pYprev);
			if (pYprev) pYprev->m_next = pX;
			else m_head = pX;
		}
		return pX;
	}

	//! Inserts element \a x before \a it.
	/**
	 * \pre \a it points to an element in this list.
	 */
	iterator insertBefore(const E &x, iterator it) {
		OGDF_ASSERT(it.valid())
		ListElement<E> *pY = it, *pX;
		ListElement<E> *pYprev = pY->m_prev;
		pY->m_prev = pX = OGDF_NEW ListElement<E>(x,pY,pYprev);
		if (pYprev) pYprev->m_next = pX;
		else m_head = pX;
		return pX;
	}

	//! Inserts element \a x after \a it.
	/**
	 * \pre \a it points to an element in this list.
	 */
	iterator insertAfter(const E &x, iterator it) {
		OGDF_ASSERT(it.valid())
		ListElement<E> *pY = it, *pX;
		ListElement<E> *pYnext = pY->m_next;
		pY->m_next = pX = OGDF_NEW ListElement<E>(x,pYnext,pY);
		if (pYnext) pYnext->m_prev = pX;
		else m_tail = pX;
		return pX;
	}


	//@}
	/**
	 * @name Removing elements
	 * These method remove elements from the list.
	 */
	//@{

	//! Removes the first element from the list.
	/**
	 * \pre The list is not empty!
	 */
	void popFront() {
		OGDF_ASSERT(m_head != nullptr)
		ListElement<E> *pX = m_head;
		m_head = m_head->m_next;
		delete pX;
		if (m_head) m_head->m_prev = nullptr;
		else m_tail = nullptr;
	}

	//! Removes the first element from the list and returns it.
	/**
	 * \pre The list is not empty!
	 */
	E popFrontRet() {
		E el = front();
		popFront();
		return el;
	}

	//! Removes the last element from the list.
	/**
	 * \pre The list is not empty!
	 */
	void popBack() {
		OGDF_ASSERT(m_tail != nullptr)
		ListElement<E> *pX = m_tail;
		m_tail = m_tail->m_prev;
		delete pX;
		if (m_tail) m_tail->m_next = nullptr;
		else m_head = nullptr;
	}

	//! Removes the last element from the list and returns it.
	/**
	 * \pre The list is not empty!
	 */
	E popBackRet() {
		E el = back();
		popBack();
		return el;
	}

	//! Removes \a it from the list.
	/**
	 * \pre \a it points to an element in this list.
	 */
	void del(iterator it) {
		OGDF_ASSERT(it.valid())
		ListElement<E> *pX = it, *pPrev = pX->m_prev, *pNext = pX->m_next;
		delete pX;
		if (pPrev) pPrev->m_next = pNext;
		else m_head = pNext;
		if (pNext) pNext->m_prev = pPrev;
		else m_tail = pPrev;
	}

	//! Removes the first occurrence of \a x (if any) from the list.
	/**
	 * If the list contains \a x several times, only the first element
	 * containing \a x is removed.
	 *
	 * \return true if one element has been removed, false otherwise.
	 */
	bool removeFirst(const E &x) {
		for(ListElement<E> *pX = m_head; pX != nullptr; pX = pX->m_next)
			if(pX->m_x == x) {
				del(pX); return true;
			}
		return false;
	}

	//! Removes all elements from the list.
	void clear() {
		if (m_head == nullptr) return;

		if (doDestruction((E*)nullptr)) {
			for(ListElement<E> *pX = m_head; pX != nullptr; pX = pX->m_next)
				pX->m_x.~E();
		}
		OGDF_ALLOCATOR::deallocateList(sizeof(ListElement<E>),m_head,m_tail);

		m_head = m_tail = nullptr;
	}


	//@}
	/**
	 * @name Moving elements
	 * The method allow to change the order of elements within the list, or to move elements to another list.
	 */
	//@{

	//! Exchanges the positions of \a it1 and \a it2 in the list.
	/**
	 * \pre \a it1 and \a it2 point to elements in this list.
	 */
	void exchange(iterator it1, iterator it2) {
		OGDF_ASSERT(it1.valid());
		OGDF_ASSERT(it2.valid());
		OGDF_ASSERT(it1 != it2);
		ListElement<E> *pX = it1, *pY = it2;

		std::swap(pX->m_next,pY->m_next);
		std::swap(pX->m_prev,pY->m_prev);

		if(pX->m_next == pX) {
			pX->m_next = pY; pY->m_prev = pX;
		}
		if(pX->m_prev == pX) {
			pX->m_prev = pY; pY->m_next = pX;
		}

		if(pX->m_prev) pX->m_prev->m_next = pX;
		else m_head = pX;

		if(pY->m_prev) pY->m_prev->m_next = pY;
		else m_head = pY;

		if(pX->m_next) pX->m_next->m_prev = pX;
		else m_tail = pX;

		if(pY->m_next) pY->m_next->m_prev = pY;
		else m_tail = pY;
	}

	//! Moves \a it to the begin of the list.
	/**
	 * \pre \a it points to an element in this list.
	 */
	void moveToFront(iterator it) {
		OGDF_ASSERT(it.valid())
		// remove it
		ListElement<E> *pX = it, *pPrev = pX->m_prev, *pNext = pX->m_next;
		//already at front
		if (!pPrev) return;

		//update old position
		if (pPrev) pPrev->m_next = pNext;
		if (pNext) pNext->m_prev = pPrev;
		else m_tail = pPrev;
		// insert it at front
		pX->m_prev = nullptr;
		pX->m_next = m_head;
		m_head = m_head->m_prev = pX;
	}

	//! Moves \a it to the end of the list.
	/**
	 * \pre \a it points to an element in this list.
	 */
	void moveToBack(iterator it) {
		OGDF_ASSERT(it.valid())
		// remove it
		ListElement<E> *pX = it, *pPrev = pX->m_prev, *pNext = pX->m_next;
		//already at back
		if (!pNext) return;

		//update old position
		if (pPrev) pPrev->m_next = pNext;
		else m_head = pNext;
		if (pNext) pNext->m_prev = pPrev;
		// insert it at back
		pX->m_prev = m_tail;
		pX->m_next = nullptr;
		m_tail = m_tail->m_next = pX;
	}

	//! Moves \a it after \a itBefore.
	/**
	 * \pre \a it and \a itBefore point to elements in this list.
	 */
	void moveToSucc(iterator it, iterator itBefore) {
		OGDF_ASSERT(it.valid());
		OGDF_ASSERT(itBefore.valid());
		// move it
		ListElement<E> *pX = it, *pPrev = pX->m_prev, *pNext = pX->m_next;
		//the same of already in place
		ListElement<E> *pY = itBefore;
		if(pX == pY || pPrev == pY) return;

		// update old position
		if (pPrev) pPrev->m_next = pNext;
		else m_head = pNext;
		if (pNext) pNext->m_prev = pPrev;
		else m_tail = pPrev;
		// move it after itBefore
		ListElement<E> *pYnext = pX->m_next = pY->m_next;
		(pX->m_prev = pY)->m_next = pX;
		if (pYnext) pYnext->m_prev = pX;
		else m_tail = pX;
	}

	//! Moves \a it before \a itAfter.
	/**
	 * \pre \a it and \a itAfter point to elements in this list.
	 */
	void moveToPrec(iterator it, iterator itAfter) {
		OGDF_ASSERT(it.valid());
		OGDF_ASSERT(itAfter.valid());
		// move it
		ListElement<E> *pX = it, *pPrev = pX->m_prev, *pNext = pX->m_next;
		//the same of already in place
		ListElement<E> *pY = itAfter;
		if(pX == pY || pNext == pY) return;

		// update old position
		if (pPrev) pPrev->m_next = pNext;
		else m_head = pNext;
		if (pNext) pNext->m_prev = pPrev;
		else m_tail = pPrev;
		// move it before itAfter
		ListElement<E> *pYprev = pX->m_prev = pY->m_prev;
		(pX->m_next = pY)->m_prev = pX;
		if (pYprev) pYprev->m_next = pX;
		else m_head = pX;
	}

	//! Moves \a it to the begin of \a L2.
	/**
	 * \pre \a it points to an element in this list.
	 */
	void moveToFront(iterator it, ListPure<E> &L2) {
		OGDF_ASSERT(it.valid())
		OGDF_ASSERT(this != &L2)
		// remove it
		ListElement<E> *pX = it, *pPrev = pX->m_prev, *pNext = pX->m_next;
		if (pPrev) pPrev->m_next = pNext;
		else m_head = pNext;
		if (pNext) pNext->m_prev = pPrev;
		else m_tail = pPrev;
		// insert it at front of L2
		pX->m_prev = nullptr;
		if ((pX->m_next = L2.m_head) != nullptr)
			L2.m_head = L2.m_head->m_prev = pX;
		else
			L2.m_head = L2.m_tail = pX;
	}

	//! Moves \a it to the end of \a L2.
	/**
	 * \pre \a it points to an element in this list.
	 */
	void moveToBack(iterator it, ListPure<E> &L2) {
		OGDF_ASSERT(it.valid())
		OGDF_ASSERT(this != &L2)
		// remove it
		ListElement<E> *pX = it, *pPrev = pX->m_prev, *pNext = pX->m_next;
		if (pPrev) pPrev->m_next = pNext;
		else m_head = pNext;
		if (pNext) pNext->m_prev = pPrev;
		else m_tail = pPrev;
		// insert it at back of L2
		pX->m_next = nullptr;
		if ((pX->m_prev = L2.m_tail) != nullptr)
			L2.m_tail = L2.m_tail->m_next = pX;
		else
			L2.m_head = L2.m_tail = pX;
	}

	//! Moves \a it to list \a L2 and inserts it after \a itBefore.
	/**
	 * \pre \a it points to an element in this list, and \a itBefore
	 *      points to an element in \a L2.
	 */
	void moveToSucc(iterator it, ListPure<E> &L2, iterator itBefore) {
		OGDF_ASSERT(it.valid());
		OGDF_ASSERT(itBefore.valid());
		OGDF_ASSERT(this != &L2)
		// remove it
		ListElement<E> *pX = it, *pPrev = pX->m_prev, *pNext = pX->m_next;
		if (pPrev) pPrev->m_next = pNext;
		else m_head = pNext;
		if (pNext) pNext->m_prev = pPrev;
		else m_tail = pPrev;
		// insert it in list L2 after itBefore
		ListElement<E> *pY = itBefore;
		ListElement<E> *pYnext = pX->m_next = pY->m_next;
		(pX->m_prev = pY)->m_next = pX;
		if (pYnext) pYnext->m_prev = pX;
		else L2.m_tail = pX;
	}

	//! Moves \a it to list \a L2 and inserts it before \a itAfter.
	/**
	 * \pre \a it points to an element in this list, and \a itAfter
	 *      points to an element in \a L2.
	 */
	void moveToPrec(iterator it, ListPure<E> &L2, iterator itAfter) {
		OGDF_ASSERT(it.valid());
		OGDF_ASSERT(itAfter.valid());
		OGDF_ASSERT(this != &L2)
		// remove it
		ListElement<E> *pX = it, *pPrev = pX->m_prev, *pNext = pX->m_next;
		if (pPrev) pPrev->m_next = pNext;
		else m_head = pNext;
		if (pNext) pNext->m_prev = pPrev;
		else m_tail = pPrev;
		// insert it in list L2 after itBefore
		ListElement<E> *pY = itAfter;
		ListElement<E> *pYprev = pX->m_prev = pY->m_prev;
		(pX->m_next = pY)->m_prev = pX;
		if (pYprev) pYprev->m_next = pX;
		else L2.m_head = pX;
	}

	//! Appends \a L2 to this list and makes \a L2 empty.
	void conc(ListPure<E> &L2) {
		OGDF_ASSERT(this != &L2)
		if (m_head)
			m_tail->m_next = L2.m_head;
		else
			m_head = L2.m_head;
		if (L2.m_head) {
			L2.m_head->m_prev = m_tail;
			m_tail = L2.m_tail;
		}
		L2.m_head = L2.m_tail = nullptr;
	}

	//! Prepends \a L2 to this list and makes \a L2 empty.
	void concFront(ListPure<E> &L2) {
		OGDF_ASSERT(this != &L2)
		if (m_head)
			m_head->m_prev = L2.m_tail;
		else
			m_tail = L2.m_tail;
		if (L2.m_head) {
			L2.m_tail->m_next = m_head;
			m_head = L2.m_head;
		}
		L2.m_head = L2.m_tail = nullptr;
	}

	//! Exchanges the contents of this list and \a other in constant time.
	void swap(ListPure<E> &other) {
		std::swap(m_head, other.m_head);
		std::swap(m_tail, other.m_tail);
	}

	//! Deprecated, use swap instead.
	/**
	 * \deprecated This method has been marked as deprecated and might be removed in a future version of the library.
	 *             Use the equivalent swap() method instead.
	 */
	OGDF_DEPRECATED_BEGIN
	void exchange(List<E>& L2)
	OGDF_DEPRECATED_END
		{ swap(L2); }

	//! Splits the list at element \a it into lists \a L1 and \a L2.
	/**
	 * If \a it is not a null pointer and \a L = x1,...,x{k-1}, \a it,x_{k+1},xn, then
	 * \a L1 = x1,...,x{k-1} and \a L2 = \a it,x{k+1},...,xn if \a dir = \c before.
	 * If \a it is a null pointer, then \a L1 is made empty and \a L2 = \a L. Finally
	 * \a L is made empty if it is not identical to \a L1 or \a L2.
	 *
	 * \pre \a it points to an element in this list.
	 */

	void split(iterator it,ListPure<E> &L1,ListPure<E> &L2,Direction dir = before) {
		if (&L1 != this) L1.clear();
		if (&L2 != this) L2.clear();

		if (it.valid()){
			L1.m_head = m_head;
			L2.m_tail = m_tail;
			if (dir == before){
				L2.m_head = it;
				L1.m_tail = L2.m_head->m_prev;
			}
			else {
				L1.m_tail = it;
				L2.m_head = L1.m_tail->m_next;
			}
			L2.m_head->m_prev = L1.m_tail->m_next = nullptr;

		} else {
			L1.m_head = L1.m_tail = nullptr;
			L2.m_head = m_head;
			L2.m_tail = m_tail;
		}

		if (this != &L1 && this != &L2) {
			m_head = m_tail = nullptr;
		}
	}

	//! Splits the list after \a it.
	void splitAfter(iterator it, ListPure<E> &L2) {
		OGDF_ASSERT(it.valid())
		OGDF_ASSERT(this != &L2)
		L2.clear();
		ListElement<E> *pX = it;
		if (pX != m_tail) {
			(L2.m_head = pX->m_next)->m_prev = nullptr;
			pX->m_next = nullptr;
			L2.m_tail = m_tail;
			m_tail = pX;
		}
	}

	//! Splits the list before \a it.
	void splitBefore(iterator it, ListPure<E> &L2) {
		OGDF_ASSERT(it.valid())
		OGDF_ASSERT(this != &L2)
		L2.clear();
		ListElement<E> *pX = it;
		L2.m_head = pX; L2.m_tail = m_tail;
		if ((m_tail = pX->m_prev) == nullptr)
			m_head = nullptr;
		else
			m_tail->m_next = nullptr;
		pX->m_prev = nullptr;
	}

	//! Reverses the order of the list elements.
	void reverse() {
		ListElement<E> *pX = m_head;
		m_head = m_tail;
		m_tail = pX;
		while(pX) {
			ListElement<E> *pY = pX->m_next;
			pX->m_next = pX->m_prev;
			pX = pX->m_prev = pY;
		}
	}


	//@}
	/**
	 * @name Searching and sorting
	 * These methods provide searching for values and sorting the list.
	 */
	//@{

	//! Scans the list for the specified element and returns an iterator to the first occurrence in the list, or an invalid iterator if not found.
	ListConstIterator<E> search(const E& e) const {
		ListConstIterator<E> i;
		for (i = begin(); i.valid(); ++i)
			if (*i == e) return i;
		return i;
	}

	//! Scans the list for the specified element and returns an iterator to the first occurrence in the list, or an invalid iterator if not found.
	ListIterator<E> search(const E& e) {
		ListIterator<E> i;
		for (i = begin(); i.valid(); ++i)
			if (*i == e) return i;
		return i;
	}

	//! Scans the list for the specified element (using the user-defined comparer) and returns an iterator to the first occurrence in the list, or an invalid iterator if not found.
	template<class COMPARER>
	ListConstIterator<E> search(const E &e, const COMPARER &comp) const {
		ListConstIterator<E> i;
		for (i = begin(); i.valid(); ++i)
			if (comp.equal(*i, e)) return i;
		return i;
	}

	//! Scans the list for the specified element (using the user-defined comparer) and returns an iterator to the first occurrence in the list, or an invalid iterator if not found.
	template<class COMPARER>
	ListIterator<E> search(const E &e, const COMPARER &comp) {
		ListIterator<E> i;
		for (i = begin(); i.valid(); ++i)
			if (comp.equal(*i, e)) return i;
		return i;
	}

	//! Sorts the list using Quicksort.
	void quicksort() {
		ogdf::quicksortTemplate(*this);
	}

	//! Sorts the list using Quicksort and comparer \a comp.
	template<class COMPARER>
	void quicksort(const COMPARER &comp) {
		ogdf::quicksortTemplate(*this,comp);
	}

	//! Sorts the list using bucket sort.
	/**
	 * @param l is the lowest bucket that will occur.
	 * @param h is the highest bucket that will occur.
	 * @param f returns the bucket for each element.
	 * \pre The bucket function \a f will only return bucket values between \a l
	 * and \a h for this list.
	 */
	void bucketSort(int l, int h, BucketFunc<E> &f);


	//@}
	/**
	 * @name Random elements and permutations
	 * These methods allow to select a random element in the list, or to randomly permute the list.
	 */
	//@{

	//! Returns a const iterator to a random element in the list (or an invalid iterator if the list is empty).
	/**
	 * This method takes linear time.
	 */
	const_iterator chooseIterator() const {
		return empty() ? const_iterator() : get(randomNumber(0,size()-1));
	}

	//! Returns an iterator to a random element in the list (or an invalid iterator if the list is empty).
	/**
	 * This method takes linear time.
	 */
	iterator chooseIterator() {
		return empty() ? iterator() : get(randomNumber(0,size()-1));
	}

	//! Returns a random element from the list.
	/**
	 * \pre The list is not empty!
	 *
	 * This method takes linear time.
	 */
	const_reference chooseElement() const {
		OGDF_ASSERT(m_head != nullptr)
		return *chooseIterator();
	}

	//! Returns a random element from the list.
	/**
	 * \pre The list is not empty!
	 *
	 * This method takes linear time.
	 */
	reference chooseElement() {
		return *chooseIterator();
	}

	//! Randomly permutes the elements in the list.
	void permute() {
		std::minstd_rand rng(randomSeed());
		permute(size(), rng);
	}

	//! Randomly permutes the elements in the list using random number generator \a rng.
	template<class RNG>
	void permute(RNG &rng) {
		permute(size(), rng);
	}

	//@}

protected:
	void copy(const ListPure<E> &L) {
		for(ListElement<E> *pX = L.m_head; pX != nullptr; pX = pX->m_next)
			pushBack(pX->m_x);
	}

	template<class RNG>
	void permute(const int n, RNG &rng);

	OGDF_NEW_DELETE
}; // class ListPure



//! Doubly linked lists (maintaining the length of the list).
/**
 * @ingroup containers
 *
 * Elements of the list are instances of type ListElement.
 * Use ListConstIterator or ListIterator in order to iterate over the list.
 *
 * In contrast to ListPure, instances of List store the length of the list.
 *
 * @tparam E is the data type stored in list elements.
 */
template<class E>
class List : private ListPure<E> {

	int m_count; //!< The length of the list.

public:
	//! Represents the data type stored in a list element.
	typedef E value_type;
	//! Provides a reference to an element stored in a list.
	typedef E &reference;
	//! Provides a reference to a const element stored in a list for reading and performing const operations.
	typedef const E &const_reference;
	//! Provides a bidirectional iterator that can read a const element in a list.
	typedef ListConstIterator<E> const_iterator;
	//! Provides a bidirectional iterator that can read or modify any element in a list.
	typedef ListIterator<E> iterator;

	//! Constructs an empty doubly linked list.
	List() : m_count(0) { }

	//! Constructs a doubly linked list containing the elements in \a init.
	List(std::initializer_list<E> init) : ListPure<E>(init), m_count((int)init.size()) { }

	//! Constructs a doubly linked list that is a copy of \a L.
	List(const List<E> &L) : ListPure<E>(L), m_count(L.m_count) { }

	//! Constructs a doubly linked list containing the elements of \a L (move semantics).
	/**
	 * The list \a L is empty afterwards.
	 */
	List(List<E> &&L) : ListPure<E>(std::move(L)), m_count(L.m_count) {
		L.m_count = 0;
	}

	//! Destructor.
	~List() { }


	/**
	 * @name Access methods
	 * These methods provide simple access without changing the list.
	 */
	//@{

	//! Returns true iff the list is empty.
	bool empty() const { return ListPure<E>::empty(); }

	//! Returns the number of elements in the list.
	/**
	 * This method has constant runtime (in contrast to ListPure::size()), since the list maintains the current size.
	 */
	int size() const { return m_count; }

	//! Returns a const reference to the first element.
	/**
	 * \pre The list is not empty!
	 */
	const_reference front() const { return ListPure<E>::front(); }

	//! Returns a reference to the first element.
	/**
	 * \pre The list is not empty!
	 */
	reference front() { return ListPure<E>::front(); }

	//! Returns a const reference to the last element.
	/**
	 * \pre The list is not empty!
	 */
	const_reference back() const { return ListPure<E>::back(); }

	//! Returns a reference to the last element.
	/**
	 * \pre The list is not empty!
	 */
	reference back() { return ListPure<E>::back(); }

	//! Returns a const iterator pointing to the element at position \a pos.
	/**
	 * The running time of this method is linear in \a pos.
	 */
	const_iterator get(int pos) const {
		OGDF_ASSERT(0 <= pos);
		OGDF_ASSERT(pos < m_count)
		return ListPure<E>::get(pos);
	}

	//! Returns an iterator pointing to the element at position \a pos.
	/**
	 * The running time of this method is linear in \a pos.
	 */
	iterator get(int pos) {
		OGDF_ASSERT(0 <= pos);
		OGDF_ASSERT(pos < m_count)
		return ListPure<E>::get(pos);
	}

	//! Returns the position (starting with 0) of iterator \a it in the list.
	/**
	 * \pre \a it is a valid iterator pointing to an element in this list!
	 */
	int pos(const_iterator it) const {
		OGDF_ASSERT(it.valid())
		return ListPure<E>::pos(it);
	}

	//! Conversion to const ListPure.
	const ListPure<E> &getListPure() const { return *this; }


	//@}
	/**
	 * @name Iterators
	 * These methods return bidirectional iterators to elements in the list and allow to iterate over the elements in linear or cyclic order.
	 */
	//@{

	//! Returns an iterator to the first element of the list.
	/**
	 * If the list is empty, a null pointer iterator is returned.
	 */
	iterator begin() { return ListPure<E>::begin(); }

	//! Returns a const iterator to the first element of the list.
	/**
	 * If the list is empty, a null pointer iterator is returned.
	 */
	const_iterator begin() const { return ListPure<E>::begin(); }

	//! Returns a const iterator to the first element of the list.
	/**
	 * If the list is empty, a null pointer iterator is returned.
	 */
	const_iterator cbegin() const { return ListPure<E>::cbegin(); }

	//! Returns an iterator to one-past-last element of the list.
	/**
	 * This is always a null pointer iterator.
	 */
	iterator end() { return iterator(); }

	//! Returns a const iterator to one-past-last element of the list.
	/**
	 * This is always a null pointer iterator.
	 */
	const_iterator end() const { return const_iterator(); }

	//! Returns a const iterator to one-past-last element of the list.
	/**
	 * This is always a null pointer iterator.
	 */
	const_iterator cend() const { return const_iterator(); }

	//! Returns an iterator to the last element of the list.
	/**
	 * If the list is empty, a null pointer iterator is returned.
	 */
	iterator rbegin() { return ListPure<E>::rbegin(); }

	//! Returns a const iterator to the last element of the list.
	/**
	 * If the list is empty, a null pointer iterator is returned.
	 */
	const_iterator rbegin() const { return ListPure<E>::rbegin(); }

	//! Returns a const iterator to the last element of the list.
	/**
	 * If the list is empty, a null pointer iterator is returned.
	 */
	const_iterator crbegin() const { return ListPure<E>::crbegin(); }

	//! Returns an iterator to one-before-first element of the list.
	/**
	 * This is always a null pointer iterator.
	 */
	iterator rend() { return iterator(); }

	//! Returns a const iterator to one-before-first element of the list.
	/**
	 * This is always a null pointer iterator.
	 */
	const_iterator rend() const { return const_iterator(); }

	//! Returns a const iterator to one-before-first element of the list.
	/**
	 * This is always a null pointer iterator.
	 */
	const_iterator crend() const { return const_iterator(); }

	//! Returns a const iterator to the cyclic successor of \a it.
	/**
	 * \pre \a it points to an element in this list!
	 */
	const_iterator cyclicSucc(const_iterator it) const {
		return ListPure<E>::cyclicSucc(it);
	}

	//! Returns an iterator to the cyclic successor of \a it.
	/**
	 * \pre \a it points to an element in this list!
	 */
	iterator cyclicSucc(iterator it) {
		return ListPure<E>::cyclicSucc(it);
	}

	//! Returns a const iterator to the cyclic predecessor of \a it.
	/**
	 * \pre \a it points to an element in this list!
	 */
	const_iterator cyclicPred(const_iterator it) const {
		return ListPure<E>::cyclicPred(it);
	}

	//! Returns an iterator to the cyclic predecessor of \a it.
	/**
	 * \pre \a it points to an element in this list!
	 */
	iterator cyclicPred(iterator it) {
		return ListPure<E>::cyclicPred(it);
	}


	//@}
	/**
	 * @name Operators
	 * The following operators are provided by lists.
	 */
	//@{

	//! Assignment operator
	List<E> &operator=(const List<E> &L) {
		ListPure<E>::operator=(L);
		m_count = L.m_count;
		return *this;
	}

	//! Assignment operator (move semantics)
	/**
	 * The list \a L is empty afterwards.
	 */
	List<E> &operator=(List<E> &&L) {
		m_count = L.m_count;
		ListPure<E>::operator=(std::move(L));
		L.m_count = 0;
		return *this;
	}

	//! Equality operator.
	bool operator==(const List<E> &L) const {
		return (m_count == L.m_count) && ListPure<E>::operator==(L);
	}

	//! Inequality operator.
	bool operator!=(const List<E> &L) const {
		return !operator==(L);
	}


	//@}
	/**
	 * @name Adding elements
	 * These method add elements to the list.
	 */
	//@{

	//! Adds element \a x at the beginning of the list.
	iterator pushFront(const E &x) {
		++m_count;
		return ListPure<E>::pushFront(x);
	}

	//! Adds a new element at the beginning of the list.
	/**
	* The element is constructed in-place with exactly the same arguments \a args as supplied to the function.
	*/
	template<class ... Args>
	iterator emplaceFront(Args && ... args) {
		++m_count;
		return ListPure<E>::emplaceFront(std::forward<Args>(args)...);
	}

	//! Adds element \a x at the end of the list.
	iterator pushBack(const E &x) {
		++m_count;
		return ListPure<E>::pushBack(x);
	}

	//! Adds a new element at the end of the list.
	/**
	* The element is constructed in-place with exactly the same arguments \a args as supplied to the function.
	*/
	template<class ... Args>
	iterator emplaceBack(Args && ... args) {
		++m_count;
		return ListPure<E>::emplaceBack(std::forward<Args>(args)...);
	}

	//! Inserts element \a x before or after \a it.
	/**
	 * @param x is the element to be inserted.
	 * @param it is a list iterator in this list.
	 * @param dir determines if \a x is inserted before or after \a it.
	 *   Possible values are \c ogdf::before and \c ogdf::after.
	 * \pre \a it points to an element in this list.
	 */
	iterator insert(const E &x, iterator it, Direction dir = after) {
		++m_count;
		return ListPure<E>::insert(x,it,dir);
	}

	//! Inserts element \a x before \a it.
	/**
	 * \pre \a it points to an element in this list.
	 */
	iterator insertBefore(const E &x, iterator it) {
		++m_count;
		return ListPure<E>::insertBefore(x,it);
	}

	//! Inserts element \a x after \a it.
	/**
	 * \pre \a it points to an element in this list.
	 */
	iterator insertAfter(const E &x, iterator it) {
		++m_count;
		return ListPure<E>::insertAfter(x,it);
	}


	//@}
	/**
	 * @name Removing elements
	 * These method remove elements from the list.
	 */
	//@{

	//! Removes the first element from the list.
	/**
	 * \pre The list is not empty!
	 */
	void popFront() {
		--m_count;
		ListPure<E>::popFront();
	}

	//! Removes the first element from the list and returns it.
	/**
	 * \pre The list is not empty!
	 */
	E popFrontRet() {
		E el = front();
		popFront();
		return el;
	}

	//! Removes the last element from the list.
	/**
	 * \pre The list is not empty!
	 */
	void popBack() {
		--m_count;
		ListPure<E>::popBack();
	}

	//! Removes the last element from the list and returns it.
	/**
	 * \pre The list is not empty!
	 */
	E popBackRet() {
		E el = back();
		popBack();
		return el;
	}

	//! Removes \a it from the list.
	/**
	 * \pre \a it points to an element in this list.
	 */
	void del(iterator it) {
		--m_count;
		ListPure<E>::del(it);
	}

	//! Removes the first occurrence of \a x (if any) from the list.
	/**
	 * If the list contains \a x several times, only the first element
	 * containing \a x is removed.
	 */
	bool removeFirst(const E &x) {
		bool hasRemoved = ListPure<E>::removeFirst(x);
		if(hasRemoved)
			--m_count;
		return hasRemoved;
	}

	//! Removes all elements from the list.
	void clear() {
		m_count = 0;
		ListPure<E>::clear();
	}


	//@}
	/**
	 * @name Moving elements
	 * The method allow to change the order of elements within the list, or to move elements to another list.
	 */
	//@{

	//! Exchanges the positions of \a it1 and \a it2 in the list.
	/**
	 * \pre \a it1 and \a it2 point to elements in this list.
	 */
	void exchange(iterator it1, iterator it2) {
		ListPure<E>::exchange(it1,it2);
	}

	//! Moves \a it to the beginning of the list
	/**
	 * \pre \a it points to an element in the list.
	 */
	void moveToFront(iterator it) {
		ListPure<E>::moveToFront(it);
	}
	//! Moves \a it to the end of the list
	/**
	 * \pre \a it points to an element in the list.
	 */
	void moveToBack(iterator it) {
		ListPure<E>::moveToBack(it);
	}
	//! Moves \a it after \a itBefore.
	/**
	 * \pre \a it and \a itBefore point to elements in this list.
	 */
	void moveToSucc(iterator it, iterator itBefore) {
		ListPure<E>::moveToSucc(it,itBefore);
	}
	//! Moves \a it before \a itAfter.
	/**
	 * \pre \a it and \a itAfter point to elements in this list.
	 */
	void moveToPrec(iterator it, iterator itAfter) {
		ListPure<E>::moveToPrec(it,itAfter);
	}

	//! Moves \a it to the beginning of \a L2.
	/**
	 * \pre \a it points to an element in this list.
	 */
	void moveToFront(iterator it, List<E> &L2) {
		ListPure<E>::moveToFront(it,L2);
		--m_count; ++L2.m_count;
	}
	//! Moves \a it to the end of \a L2.
	/**
	 * \pre \a it points to an element in this list.
	 */
	void moveToBack(iterator it, List<E> &L2) {
		ListPure<E>::moveToBack(it,L2);
		--m_count; ++L2.m_count;
	}

	//! Moves \a it to list \a L2 and inserts it after \a itBefore.
	/**
	 * \pre \a it points to an element in this list, and \a itBefore
	 *      points to an element in \a L2.
	 */
	void moveToSucc(iterator it, List<E> &L2, iterator itBefore) {
		ListPure<E>::moveToSucc(it,L2,itBefore);
		--m_count; ++L2.m_count;
	}
	//! Moves \a it to list \a L2 and inserts it after \a itBefore.
	/**
	 * \pre \a it points to an element in this list, and \a itBefore
	 *      points to an element in \a L2.
	 */
	void moveToPrec(iterator it, List<E> &L2, iterator itAfter) {
		ListPure<E>::moveToPrec(it,L2,itAfter);
		--m_count; ++L2.m_count;
	}


	//! Appends \a L2 to this list and makes \a L2 empty.
	void conc(List<E> &L2) {
		ListPure<E>::conc(L2);
		m_count += L2.m_count;
		L2.m_count = 0;
	}

	//! Prepends \a L2 to this list and makes \a L2 empty.
	void concFront(List<E> &L2) {
		ListPure<E>::concFront(L2);
		m_count += L2.m_count;
		L2.m_count = 0;
	}

	//! Exchanges the contents of this list and \a other in constant time.
	void swap(List<E> &other) {
		ListPure<E>::swap(other);
		std::swap(m_count, other.m_count);
	}

	//! Deprecated, use swap instead.
	/**
	* \deprecated This method has been marked as deprecated and might be removed in a future version of the library.
	*             Use the equivalent swap() method instead.
	*/
	OGDF_DEPRECATED_BEGIN
	void exchange(List<E>& L2)
	OGDF_DEPRECATED_END
		{ swap(L2); }

	//! Splits the list at element \a it into lists \a L1 and \a L2.
	/**
	 * If \a it is not a null pointer and \a L = x1,...,x{k-1}, \a it,x_{k+1},xn, then
	 * \a L1 = x1,...,x{k-1} and \a L2 = \a it,x{k+1},...,xn if \a dir = \c before.
	 * If \a it is a null pointer, then \a L1 is made empty and \a L2 = \a L. Finally
	 * \a L is made empty if it is not identical to \a L1 or \a L2.
	 *
	 * \pre \a it points to an element in this list.
	 */
	void split(iterator it, List<E> &L1, List<E> &L2, Direction dir = before) {
		ListPure<E>::split(it,L1,L2,dir);
		int countL = m_count, countL1 = 0;
		for(ListElement<E> *pX = L1.m_head; pX != nullptr; pX = pX->m_next)
			++countL1;

		L1.m_count = countL1;
		L2.m_count = countL - countL1;
		if (this->m_head == nullptr) m_count = 0;
	}

	//! Reverses the order of the list elements.
	void reverse() { ListPure<E>::reverse(); }


	//@}
	/**
	 * @name Searching and sorting
	 * These methods provide searching for values and sorting the list.
	 */
	//@{

	//! Scans the list for the specified element and returns an iterator to the first occurrence in the list, or an invalid iterator if not found.
	ListConstIterator<E> search(const E &e) const {
		return ListPure<E>::search(e);
	}

	//! Scans the list for the specified element and returns an iterator to the first occurrence in the list, or an invalid iterator if not found.
	ListIterator<E> search(const E &e) {
		return ListPure<E>::search(e);
	}

	//! Scans the list for the specified element (using the user-defined comparer) and returns an iterator to the first occurrence in the list, or an invalid iterator if not found.
	template<class COMPARER>
	ListConstIterator<E> search(const E &e, const COMPARER &comp) const {
		return ListPure<E>::search(e, comp);
	}

	//! Scans the list for the specified element (using the user-defined comparer) and returns an iterator to the first occurrence in the list, or an invalid iterator if not found.
	template<class COMPARER>
	ListIterator<E> search(const E &e, const COMPARER &comp) {
		return ListPure<E>::search(e, comp);
	}

	//! Sorts the list using Quicksort.
	void quicksort() {
		ogdf::quicksortTemplate(*this);
	}

	//! Sorts the list using Quicksort and comparer \a comp.
	template<class COMPARER>
	void quicksort(const COMPARER &comp) {
		ogdf::quicksortTemplate(*this,comp);
	}

	//! Sorts the list using bucket sort.
	/**
	 * @param l is the lowest bucket that will occur.
	 * @param h is the highest bucket that will occur.
	 * @param f returns the bucket for each element.
	 * \pre The bucket function \a f will only return bucket values between \a l
	 * and \a h for this list.
	 */
	void bucketSort(int l, int h, BucketFunc<E> &f) {
		ListPure<E>::bucketSort(l,h,f);
	}


	//@}
	/**
	 * @name Random elements and permutations
	 * These methods allow to select a random element in the list, or to randomly permute the list.
	 */
	//@{

	//! Returns a const iterator to a random element in the list (or an invalid iterator if the list is empty).
	/**
	 * This method takes linear time.
	 */
	const_iterator chooseIterator() const {
		return (m_count > 0) ? get(randomNumber(0,m_count-1)) : const_iterator();
	}

	//! Returns an iterator to a random element in the list (or an invalid iterator if the list is empty).
	/**
	 * This method takes linear time.
	 */
	iterator chooseIterator() {
		return (m_count > 0) ? get(randomNumber(0,m_count-1)) : iterator();
	}

	//! Returns a random element from the list.
	/**
	 * \pre The list is not empty!
	 *
	 * This method takes linear time.
	 */
	const_reference chooseElement() const {
		OGDF_ASSERT(!empty());
		return *chooseIterator();
	}

	//! Returns a random element from the list.
	/**
	 * \pre The list is not empty!
	 *
	 * This method takes linear time.
	 */
	reference chooseElement() {
		OGDF_ASSERT(!empty());
		return *chooseIterator();
	}

	//! Randomly permutes the elements in the list.
	void permute() {
		std::minstd_rand rng(randomSeed());
		permute(rng);
	}

	//! Randomly permutes the elements in the list using random number generator \a rng.
	template<class RNG>
	void permute(RNG &rng) {
		ListPure<E>::permute(m_count,rng);
	}

	//@}

	OGDF_NEW_DELETE
}; // class List



template<class E>
void ListPure<E>::bucketSort(int l, int h, BucketFunc<E> &f)
{
	if (m_head == m_tail) return;

	Array<ListElement<E> *> head(l,h,0), tail(l,h);

	ListElement<E> *pX;
	for (pX = m_head; pX; pX = pX->m_next) {
		int i = f.getBucket(pX->m_x);
		if (head[i])
			tail[i] = ((pX->m_prev = tail[i])->m_next = pX);
		else
			head[i] = tail[i] = pX;
	}

	ListElement<E> *pY = nullptr;
	for (int i = l; i <= h; i++) {
		pX = head[i];
		if (pX) {
			if (pY) {
				(pY->m_next = pX)->m_prev = pY;
			} else
				(m_head = pX)->m_prev = nullptr;
			pY = tail[i];
		}
	}

	m_tail = pY;
	pY->m_next = nullptr;
}


// permutes elements in list randomly; n is the length of the list
template<class E>
template<class RNG>
void ListPure<E>::permute(const int n, RNG &rng)
{
	//if n==0 do nothing
	if (n == 0) { return; }

	Array<ListElement<E> *> A(n+2);
	A[0] = A[n+1] = nullptr;

	int i = 1;
	ListElement<E> *pX;
	for (pX = m_head; pX; pX = pX->m_next)
		A[i++] = pX;

	A.permute(1,n,rng);

	for (i = 1; i <= n; i++) {
		pX = A[i];
		pX->m_next = A[i+1];
		pX->m_prev = A[i-1];
	}

	m_head = A[1];
	m_tail = A[n];
}


// prints list L to output stream os using delimiter delim
template<class E>
void print(ostream &os, const ListPure<E> &L, char delim = ' ')
{
	typename ListPure<E>::const_iterator pX = L.begin();
	if (pX.valid()) {
		os << *pX;
		for(++pX; pX.valid(); ++pX)
			os << delim << *pX;
	}
}

// prints list L to output stream os using delimiter delim
template<class E>
void print(ostream &os, const List<E> &L, char delim = ' ')
{
	print(os, L.getListPure(), delim);
}

// prints list L to output stream os
template<class E>
ostream &operator<<(ostream &os, const ListPure<E> &L)
{
	print(os,L);
	return os;
}

// prints list L to output stream os
template<class E>
ostream &operator<<(ostream &os, const List<E> &L)
{
	return os << L.getListPure();
}


template<class E, class Master>
class ListContainer : public List<E> {

	friend Master;

public:
	//! Provides a bidirectional iterator to an object in the container.
	typedef typename List<E>::const_iterator iterator;

	//! Returns an iterator to the first element in the container.
	iterator begin() const { return List<E>::cbegin(); }

	//! Returns an iterator to the one-past-last element in the container.
	iterator end() const { return List<E>::cend(); }

	//! Returns an iterator to the last element in the container.
	iterator rbegin() const { return List<E>::crbegin(); }

	//! Returns an iterator to the one-before-first element in the container.
	iterator rend() const { return List<E>::crend(); }

	//! Returns the number of elements in the container.
	int size() const { return List<E>::size(); }
};



} // end namespace ogdf


#endif
