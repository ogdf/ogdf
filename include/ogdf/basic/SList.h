/** \file
 * \brief Declaration and implementation of singly linked lists
 * (SListPure<E> and SList<E>) and iterators (SListConstIterator<E>
 * and SListIterator<E>).
 *
 * \author Carsten Gutwenger
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

#include <ogdf/internal/basic/list_templates.h>


namespace ogdf {


template<class E> class SListPure;
template<class E> class StackPure;
template<class E> class SListIterator;
template<class E> class SListConstIterator;


//! The parameterized class \a SListElement<E> represents the structure for elements of singly linked lists.
template<class E>
class SListElement {
	friend class SListPure<E>;
	friend class StackPure<E>;
	friend class SListIterator<E>;
	friend class SListConstIterator<E>;

	SListElement<E> *m_next; //!< Pointer to successor element.
	E m_x; //!< Stores the content.

	//! Constructs an SListElement.
	SListElement() : m_next(nullptr) { }
	//! Constructs an SListElement.
	SListElement(const E &x) : m_next(nullptr), m_x(x) { }
	//! Constructs an SListElement.
	SListElement(const E &x, SListElement<E> *next) :
		m_next(next), m_x(x) { }
	//! Constructs an SListElement with given arguments \a args for constructor of element type.
	template<class ... Args>
	SListElement(SListElement<E> *next, Args && ... args)
		: m_next(next), m_x(std::forward<Args>(args)...) { }

	OGDF_NEW_DELETE
}; // class SListElement



//! The parameterized class \a SListIterator<E> encapsulates a pointer to an slist element.
/**
 * It is used in order to iterate over singly linked lists,
 * and to specify a position in a singly linked list. It is possible that
 * an iterator encapsulates a null pointer.
 */

template<class E> class SListIterator {
	SListElement<E> *m_pX; //!< Pointer to slist element.

	friend class SListConstIterator<E>;
	friend class SListPure<E>;

	//! Conversion to pointer to slist element.
	operator SListElement<E> *() { return m_pX; }
	//! Conversion to pointer to slist element.
	operator const SListElement<E> *() const { return m_pX; }

public:
	//! Constructs an iterator pointing to no element.
	SListIterator() : m_pX(nullptr) { }
	//! Constructs an iterator pointing to \a pX.
	SListIterator(SListElement<E> *pX) : m_pX(pX) { }
	//! Constructs an iterator that is a copy of \a it.
	SListIterator(const SListIterator<E> &it) : m_pX(it.m_pX) { }

	//! Returns true iff the iterator points to an element.
	bool valid() const { return m_pX != nullptr; }

	//! Equality operator.
	bool operator==(const SListIterator<E> &it) const {
		return m_pX == it.m_pX;
	}

	//! Inequality operator.
	bool operator!=(const SListIterator<E> &it) const {
		return m_pX != it.m_pX;
	}

	//! Returns successor iterator.
	SListIterator<E> succ() const { return m_pX->m_next; }

	//! Returns a reference to the element content.
	E &operator*() const { return m_pX->m_x; }

	//! Assignment operator.
	SListIterator<E> &operator=(const SListIterator<E> &it) {
		m_pX = it.m_pX;
		return *this;
	}

	//! Increment operator (prefix).
	SListIterator<E> &operator++() {
		m_pX = m_pX->m_next;
		return *this;
	}

	//! Increment operator (postfix).
	SListIterator<E> operator++(int) {
		SListIterator<E> it = *this;
		m_pX = m_pX->m_next;
		return it;
	}

	OGDF_NEW_DELETE
}; // class SListIterator



//! The parameterized class \a SListIterator<E> encapsulates a constant pointer to an slist element.
/**
 * It is used in order to iterate over singly linked lists,
 * and to specify a position in a singly linked list. It is possible that
 * an iterator encapsulates a null pointer. In contrast to SListIterator,
 * it is not possible to change the slist element pointed to.
 */

template<class E> class SListConstIterator {
	const SListElement<E> *m_pX; //!< Pointer to slist element.

	friend class SListPure<E>;

	//! Conversion to pointer to slist element.
	operator const SListElement<E> *() { return m_pX; }

public:
	//! Constructs an iterator pointing to no element.
	SListConstIterator() : m_pX(nullptr) { }

	//! Constructs an iterator pointing to \a pX.
	SListConstIterator(const SListElement<E> *pX) : m_pX(pX) { }

	//! Constructs an iterator that is a copy of \a it.
	SListConstIterator(const SListIterator<E> &it) : m_pX((const SListElement<E> *)it) { }
	//! Constructs an iterator that is a copy of \a it.
	SListConstIterator(const SListConstIterator &it) : m_pX(it.m_pX) { }

	//! Returns true iff the iterator points to an element.
	bool valid() const { return m_pX != nullptr; }

	//! Equality operator.
	bool operator==(const SListConstIterator<E> &it) const {
		return m_pX == it.m_pX;
	}

	//! Inequality operator.
	bool operator!=(const SListConstIterator<E> &it) const {
		return m_pX != it.m_pX;
	}

	//! Returns successor iterator.
	SListConstIterator<E> succ() const { return m_pX->m_next; }

	//! Returns a reference to the element content.
	const E &operator*() const { return m_pX->m_x; }

	//! Assignment operator.
	SListConstIterator<E> &operator=(const SListConstIterator<E> &it) {
		m_pX = it.m_pX;
		return *this;
	}


	//! Increment operator (prefix).
	SListConstIterator<E> &operator++() {
		m_pX = m_pX->m_next;
		return *this;
	}

	//! Increment operator (postfix).
	SListConstIterator<E> operator++(int) {
		SListConstIterator<E> it = *this;
		m_pX = m_pX->m_next;
		return it;
	}

	OGDF_NEW_DELETE
}; // class SListConstIterator


//! Singly linked lists.
/**
 * @ingroup containers
 *
 * Elements of the list are instances of type SListElement.
 * Use SListConstIterator or SListIterator in order to iterate over the list.
 *
 * In contrast to SList, instances of SListPure do not store the length of the list.
 *
 * @tparam E is the data type stored in list elements.
 */

template<class E> class SListPure {
	SListElement<E> *m_head; //!< Pointer to first element.
	SListElement<E> *m_tail; //!< Pointer to last element.

public:
	//! Represents the data type stored in a list element.
	typedef E value_type;
	//! Provides a reference to an element stored in a list.
	typedef E &reference;
	//! Provides a reference to a const element stored in a list for reading and performing const operations.
	typedef const E &const_reference;
	//! Provides a forward iterator that can read a const element in a list.
	typedef SListConstIterator<E> const_iterator;
	//! Provides a forward iterator that can read or modify any element in a list.
	typedef SListIterator<E> iterator;

	//! Constructs an empty singly linked list.
	SListPure() : m_head(nullptr), m_tail(nullptr) { }

	//! Constructs a singly linked list containing the elements in \a init.
	SListPure(std::initializer_list<E> init) : m_head(nullptr), m_tail(nullptr) {
		for (const E &x : init)
			pushBack(x);
	}
	//! Constructs a singly linked list that is a copy of \a L.
	SListPure(const SListPure<E> &L) : m_head(nullptr), m_tail(nullptr) {
		copy(L);
	}

	//! Constructs a singly linked list containing the elements of \a L (move semantics).
	/**
	 * The list \a L is empty afterwards.
	 */
	SListPure(SListPure<E> &&L) : m_head(L.m_head), m_tail(L.m_tail) {
		L.m_head = L.m_tail = nullptr;
	}

	//! Destructor.
	~SListPure() { clear(); }


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
	 * If you require frequent access to the size of the list, consider using SList instead of SListPure.
	 */
	int size() const {
		int count = 0;
		for (SListElement<E> *pX = m_head; pX != nullptr; pX = pX->m_next)
			++count;
		return count;
	}

	//! Returns a reference to the first element.
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

	//! Returns a reference to the last element.
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

	//! Returns an iterator pointing to the element at position \a pos.
	/**
	 * The running time of this method is linear in \a pos.
	 */
	const_iterator get(int pos) const {
		SListElement<E> *pX;
		for(pX = m_head; pX != nullptr; pX = pX->m_next)
			if (pos-- == 0) break;
		return pX;
	}

	//! Returns an iterator pointing to the element at position \a pos.
	/**
	 * The running time of this method is linear in \a pos.
	 */
	iterator get(int pos) {
		SListElement<E> *pX;
		for(pX = m_head; pX != nullptr; pX = pX->m_next)
			if (pos-- == 0) break;
		return pX;
	}

	//! Returns the position (starting with 0) of \a it in the list.
	/**
	 * Positions are numbered 0,1,...
	 * \pre \a it is an iterator pointing to an element in this list.
	 */
	int pos(const_iterator it) const {
		OGDF_ASSERT(it.valid())
		int p = 0;
		for(SListElement<E> *pX = m_head; pX != nullptr; pX = pX->m_next, ++p)
			if (pX == it) break;
		return p;
	}


	//@}
	/**
	 * @name Iterators
	 * These methods return forward iterators to elements in the list and allow to iterate over the elements in linear or cyclic order.
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
	iterator end() { return SListIterator<E>(); }

	//! Returns a const iterator to one-past-last element of the list.
	/**
	 * This is always a null pointer iterator.
	 */
	const_iterator end() const { return SListConstIterator<E>(); }

	//! Returns a const iterator to one-past-last element of the list.
	/**
	 * This is always a null pointer iterator.
	 */
	const_iterator cend() const { return SListConstIterator<E>(); }

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

	//! Returns an iterator to the cyclic successor of \a it.
	/**
	 * \pre \a it points to an element in this list!
	 */
	const_iterator cyclicSucc(const_iterator it) const {
		const SListElement<E> *pX = it;
		return (pX->m_next) ? pX->m_next : m_head;
	}

	//! Returns an iterator to the cyclic successor of \a it.
	/**
	 * \pre \a it points to an element in this list!
	 */
	iterator cyclicSucc(iterator it) {
		SListElement<E> *pX = it;
		return (pX->m_next) ? pX->m_next : m_head;
	}


	//@}
	/**
	 * @name Operators
	 * The following operators are provided by lists.
	 */
	//@{

	//! Assignment operator.
	SListPure<E> &operator=(const SListPure<E> &L) {
		clear(); copy(L);
		return *this;
	}

	//! Assignment operator (move semantics).
	/**
	 * The list \a L is empty afterwards.
	 */
	SListPure<E> &operator=(SListPure<E> &&L) {
		clear();
		m_head = L.m_head;
		m_tail = L.m_tail;
		L.m_head = L.m_tail = nullptr;
		return *this;
	}

	//! Equality operator.
	bool operator==(const SListPure<E> &L) const {
		SListElement<E> *pX = m_head, *pY = L.m_head;
		while(pX != nullptr && pY != nullptr) {
			if(pX->m_x != pY->m_x)
				return false;
			pX = pX->m_next;
			pY = pY->m_next;
		}
		return (pX == nullptr && pY == nullptr);
	}

	//! Inequality operator.
	bool operator!=(const SListPure<E> &L) const {
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
		m_head = OGDF_NEW SListElement<E>(x,m_head);
		if (m_tail == nullptr) m_tail = m_head;
		return m_head;
	}

	//! Adds a new element at the beginning of the list.
	/**
	* The element is constructed in-place with exactly the same arguments \a args as supplied to the function.
	*/
	template<class ... Args>
	iterator emplaceFront(Args && ... args) {
		m_head = OGDF_NEW SListElement<E>(m_head, std::forward<Args>(args)...);
		if (m_tail == nullptr) m_tail = m_head;
		return m_head;
	}

	//! Adds element \a x at the end of the list.
	iterator pushBack(const E &x) {
		SListElement<E> *pNew = OGDF_NEW SListElement<E>(x);
		if (m_head == nullptr)
			m_head = m_tail = pNew;
		else
			m_tail = m_tail->m_next = pNew;
		return m_tail;
	}

	//! Adds a new element at the end of the list.
	/**
	* The element is constructed in-place with exactly the same arguments \a args as supplied to the function.
	*/
	template<class ... Args>
	iterator emplaceBack(Args && ... args) {
		SListElement<E> *pNew = OGDF_NEW SListElement<E>(nullptr, std::forward<Args>(args)...);
		if (m_head == nullptr)
			m_head = m_tail = pNew;
		else
			m_tail = m_tail->m_next = pNew;
		return m_tail;
	}

	//! Inserts element \a x after \a pBefore.
	/**
	 * \pre \a pBefore points to an element in this list.
	 */
	iterator insertAfter(const E &x, iterator itBefore) {
		SListElement<E> *pBefore = itBefore;
		OGDF_ASSERT(pBefore != nullptr)
		SListElement<E> *pNew = OGDF_NEW SListElement<E>(x,pBefore->m_next);
		if (pBefore == m_tail) m_tail = pNew;
		return (pBefore->m_next = pNew);
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
		SListElement<E> *pX = m_head;
		if ((m_head = m_head->m_next) == nullptr) m_tail = nullptr;
		delete pX;
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

	//! Removes the succesor of \a pBefore.
	/**
	 * \pre \a pBefore points to an element in this list.
	 */
	void delSucc(iterator itBefore) {
		SListElement<E> *pBefore = itBefore;
		OGDF_ASSERT(pBefore != nullptr)
		SListElement<E> *pDel = pBefore->m_next;
		OGDF_ASSERT(pDel != nullptr)
		if ((pBefore->m_next = pDel->m_next) == nullptr) m_tail = pBefore;
		delete pDel;
	}

	//! Removes all elements from the list.
	void clear() {
		if (m_head == nullptr) return;

		if (doDestruction((E*)nullptr)) {
			for(SListElement<E> *pX = m_head; pX != nullptr; pX = pX->m_next)
				pX->m_x.~E();
		}
		OGDF_ALLOCATOR::deallocateList(sizeof(SListElement<E>),m_head,m_tail);

		m_head = m_tail = nullptr;
	}


	//@}
	/**
	 * @name Moving elements
	 * The method allow to change the order of elements within the list, or to move elements to another list.
	 */
	//@{

	//! Moves the first element of this list to the begin of list \a L2.
	void moveFrontToFront(SListPure<E> &L2) {
		OGDF_ASSERT(m_head != nullptr)
		OGDF_ASSERT(this != &L2)
		SListElement<E> *pX = m_head;
		if ((m_head = m_head->m_next) == nullptr) m_tail = nullptr;
		pX->m_next = L2.m_head;
		L2.m_head = pX;
		if (L2.m_tail == nullptr) L2.m_tail = L2.m_head;
	}

	//! Moves the first element of this list to the end of list \a L2.
	void moveFrontToBack(SListPure<E> &L2) {
		OGDF_ASSERT(m_head != nullptr)
		OGDF_ASSERT(this != &L2)
		SListElement<E> *pX = m_head;
		if ((m_head = m_head->m_next) == nullptr) m_tail = nullptr;
		pX->m_next = nullptr;
		if (L2.m_head == nullptr)
			L2.m_head = L2.m_tail = pX;
		else
			L2.m_tail = L2.m_tail->m_next = pX;
	}

	//! Moves the first element of this list to list \a L2 inserted after \a itBefore.
	/**
	 * \pre \a itBefore points to an element in \a L2.
	 */
	void moveFrontToSucc(SListPure<E> &L2, iterator itBefore) {
		OGDF_ASSERT(m_head != nullptr)
		OGDF_ASSERT(this != &L2)
		SListElement<E> *pBefore = itBefore;
		SListElement<E> *pX = m_head;
		if ((m_head = m_head->m_next) == nullptr) m_tail = nullptr;
		pX->m_next = pBefore->m_next;
		pBefore->m_next = pX;
		if (pBefore == L2.m_tail) L2.m_tail = pX;
	}

	//! Appends \a L2 to this list and makes \a L2 empty.
	void conc(SListPure<E> &L2) {
		if (m_head)
			m_tail->m_next = L2.m_head;
		else
			m_head = L2.m_head;
		if (L2.m_tail != nullptr) m_tail = L2.m_tail;
		L2.m_head = L2.m_tail = nullptr;
	}

	//! Reverses the order of the list elements.
	void reverse() {
		SListElement<E> *p, *pNext, *pPred = nullptr;
		for(p = m_head; p; p = pNext) {
			pNext = p->m_next;
			p->m_next = pPred;
			pPred = p;
		}
		swap(m_head,m_tail);
	}


	//@}
	/**
	 * @name Searching and sorting
	 * These methods provide searching for values and sorting the list.
	 */
	//@{

	//! Scans the list for the specified element and returns an iterator to the first occurrence in the list, or an invalid iterator if not found.
	SListConstIterator<E> search(const E &e) const {
		SListConstIterator<E> i;
		for (i = begin(); i.valid(); ++i)
			if (*i == e) return i;
		return i;
	}

	//! Scans the list for the specified element and returns an iterator to the first occurrence in the list, or an invalid iterator if not found.
	SListIterator<E> search(const E &e) {
		SListIterator<E> i;
		for (i = begin(); i.valid(); ++i)
			if (*i == e) return i;
		return i;
	}

	//! Scans the list for the specified element (using the user-defined comparer) and returns an iterator to the first occurrence in the list, or an invalid iterator if not found.
	template<class COMPARER>
	SListConstIterator<E> search(const E &e, const COMPARER &comp) const {
		SListConstIterator<E> i;
		for (i = begin(); i.valid(); ++i)
			if (comp.equal(*i, e)) return i;
		return i;
	}

	//! Scans the list for the specified element (using the user-defined comparer) and returns an iterator to the first occurrence in the list, or an invalid iterator if not found.
	template<class COMPARER>
	SListIterator<E> search(const E &e, const COMPARER &comp) {
		SListIterator<E> i;
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

	//! Sorts the list using bucket sort.
	void bucketSort(BucketFunc<E> &f);


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
		permute(size(), rng);
	}

	//! Randomly permutes the elements in the list using random number generator \a rng.
	template<class RNG>
	void permute(RNG &rng) {
		permute(size(), rng);
	}

	//@}

protected:
	void copy(const SListPure<E> &L) {
		for(SListElement<E> *pX = L.m_head; pX != nullptr; pX = pX->m_next)
			pushBack(pX->m_x);
	}

	template<class RNG>
	void permute(const int n, RNG &rng);

	OGDF_NEW_DELETE
}; // class SListPure



//! Singly linked lists (maintaining the length of the list).
/**
 * @ingroup containers
 *
 * Elements of the list are instances of type SListElement.
 * Use SListConstIterator<E> or SListIterator in order to iterate over the list.
 * In contrast to SListPure, instances of SList store the length of the list
 * and thus allow constant time access to the length.
 *
 * @tparam E is the data type stored in list elements.
 */

template<class E>
class SList : private SListPure<E> {

	int m_count; //!< The length of the list.

public:
	//! Represents the data type stored in a list element.
	typedef E value_type;
	//! Provides a reference to an element stored in a list.
	typedef E &reference;
	//! Provides a reference to a const element stored in a list for reading and performing const operations.
	typedef const E &const_reference;
	//! Provides a bidirectional iterator that can read a const element in a list.
	typedef SListConstIterator<E> const_iterator;
	//! Provides a bidirectional iterator that can read or modify any element in a list.
	typedef SListIterator<E> iterator;

	//! Constructs an empty singly linked list.
	SList() : m_count(0) { }

	//! Constructs a singly linked list containing the elements in \a init.
	SList(std::initializer_list<E> init) : SListPure<E>(init), m_count((int) init.size()) { }

	//! Constructs a singly linked list that is a copy of \a L.
	SList(const SList<E> &L) : SListPure<E>(L), m_count(L.m_count) { }

	//! Constructs a singly linked list containing the elements of \a L (move semantics).
	/**
	 * The list \a L is empty afterwards.
	 */
	SList(SList<E> &&L) : SListPure<E>(std::move(L)), m_count(L.m_count) {
		L.m_count = 0;
	}

	//! Destructor.
	~SList() { }


	/**
	 * @name Access methods
	 * These methods provide simple access without changing the list.
	 */
	//@{

	//! Returns true iff the list is empty.
	bool empty() const { return SListPure<E>::empty(); }

	//! Returns the number of elements in the list.
	/**
	 * This method has constant runtime (in contrast to SListPure::size()), since the list maintains the current size.
	 */
	int size() const { return m_count; }

	//! Returns a reference to the first element.
	/**
	 * \pre The list is not empty!
	 */
	const_reference front() const { return SListPure<E>::front(); }

	//! Returns a reference to the first element.
	/**
	 * \pre The list is not empty!
	 */
	reference front() { return SListPure<E>::front(); }

	//! Returns a reference to the last element.
	/**
	 * \pre The list is not empty!
	 */
	const_reference back() const { return SListPure<E>::back(); }

	//! Returns a reference to the last element.
	/**
	 * \pre The list is not empty!
	 */
	reference back() { return SListPure<E>::back(); }

	//! Returns an iterator pointing to the element at position \a pos.
	/**
	 * The running time of this method is linear in \a pos.
	 */
	const_iterator get(int pos) const {
		return SListPure<E>::get(pos);
	}

	//! Returns an iterator pointing to the element at position \a pos.
	/**
	 * The running time of this method is linear in \a pos.
	 */
	iterator get(int pos) {
		return SListPure<E>::get(pos);
	}

	//! Returns the position (starting with 0) of \a it in the list.
	/**
	 * Positions are numbered 0,1,...
	 * \pre \a it is an iterator pointing to an element in this list.
	 */
	int pos(const_iterator it) const {
		return SListPure<E>::pos(it);
	}

	//! Conversion to const SListPure.
	const SListPure<E> &getSListPure() const { return *this; }


	//@}
	/**
	 * @name Iterators
	 * These methods return forward iterators to elements in the list and allow to iterate over the elements in linear or cyclic order.
	 */
	//@{

	//! Returns an iterator to the first element of the list.
	/**
	 * If the list is empty, a null pointer iterator is returned.
	 */
	iterator begin() { return SListPure<E>::begin(); }

	//! Returns a const iterator to the first element of the list.
	/**
	 * If the list is empty, a null pointer iterator is returned.
	 */
	const_iterator begin() const { return SListPure<E>::begin(); }

	//! Returns a const iterator to the first element of the list.
	/**
	 * If the list is empty, a null pointer iterator is returned.
	 */
	const_iterator cbegin() const { return SListPure<E>::begin(); }

	//! Returns an iterator to one-past-last element of the list.
	/**
	 * This is always a null pointer iterator.
	 */
	iterator end() { return SListIterator<E>(); }

	//! Returns a const iterator to one-past-last element of the list.
	/**
	 * This is always a null pointer iterator.
	 */
	const_iterator end() const { return SListConstIterator<E>(); }

	//! Returns a const iterator to one-past-last element of the list.
	/**
	 * This is always a null pointer iterator.
	 */
	const_iterator cend() const { return SListConstIterator<E>(); }

	//! Returns an iterator to the last element of the list.
	/**
	 * If the list is empty, a null pointer iterator is returned.
	 */
	iterator rbegin() { return SListPure<E>::rbegin(); }

	//! Returns a const iterator to the last element of the list.
	/**
	 * If the list is empty, a null pointer iterator is returned.
	 */
	const_iterator rbegin() const { return SListPure<E>::rbegin(); }

	//! Returns a const iterator to the last element of the list.
	/**
	 * If the list is empty, a null pointer iterator is returned.
	 */
	const_iterator crbegin() const { return SListPure<E>::rbegin(); }

	//! Returns an iterator to the cyclic successor of \a it.
	/**
	 * \pre \a it points to an element in this list!
	 */
	const_iterator cyclicSucc(const_iterator it) const {
		return SListPure<E>::cyclicSucc(it);
	}

	//! Returns an iterator to the cyclic successor of \a it.
	/**
	 * \pre \a it points to an element in this list!
	 */
	iterator cyclicSucc(iterator it) {
		return SListPure<E>::cyclicSucc(it);
	}


	//@}
	/**
	 * @name Operators
	 * The following operators are provided by lists.
	 */
	//@{

	//! Assignment operator.
	SList<E> &operator=(const SList<E> &L) {
		SListPure<E>::operator=(L);
		m_count = L.m_count;
		return *this;
	}

	//! Assignment operator (move semantics)
	/**
	 * The list \a L is empty afterwards.
	 */
	SList<E> &operator=(SList<E> &&L) {
		m_count = L.m_count;
		SListPure<E>::operator=(std::move(L));
		L.m_count = 0;
		return *this;
	}

	//! Equality operator.
	bool operator==(const SList<E> &L) const {
		return (m_count == L.m_count) && SListPure<E>::operator==(L);
	}

	//! Inequality operator.
	bool operator!=(const SList<E> &L) const {
		return !operator==(L);
	}


	//@}
	/**
	 * @name Adding elements
	 * These method add elements to the list.
	 */
	//@{

	//! Adds element \a x at the begin of the list.
	SListIterator<E> pushFront(const E &x) {
		++m_count;
		return SListPure<E>::pushFront(x);
	}

	//! Adds a new element at the beginning of the list.
	/**
	* The element is constructed in-place with exactly the same arguments \a args as supplied to the function.
	*/
	template<class ... Args>
	iterator emplaceFront(Args && ... args) {
		++m_count;
		return SListPure<E>::emplaceFront(std::forward<Args>(args)...);
	}

	//! Adds element \a x at the end of the list.
	SListIterator<E> pushBack(const E &x) {
		++m_count;
		return SListPure<E>::pushBack(x);
	}

	//! Adds a new element at the end of the list.
	/**
	* The element is constructed in-place with exactly the same arguments \a args as supplied to the function.
	*/
	template<class ... Args>
	iterator emplaceBack(Args && ... args) {
		++m_count;
		return SListPure<E>::emplaceBack(std::forward<Args>(args)...);
	}

	//! Inserts element \a x after \a pBefore.
	/**
	 * \pre \a pBefore points to an element in this list.
	 */
	SListIterator<E> insertAfter(const E &x, SListIterator<E> itBefore) {
		++m_count;
		return SListPure<E>::insertAfter(x, itBefore);
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
		SListPure<E>::popFront();
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

	//! Removes the succesor of \a pBefore.
	/**
	 * \pre \a pBefore points to an element in this list.
	 */
	void delSucc(SListIterator<E> itBefore) {
		--m_count;
		SListPure<E>::delSucc(itBefore);
	}

	//! Removes all elements from the list.
	void clear() {
		m_count = 0;
		SListPure<E>::clear();
	}


	//@}
	/**
	 * @name Moving elements
	 * The method allow to change the order of elements within the list, or to move elements to another list.
	 */
	//@{

	//! Moves the first element of this list to the begin of list \a L2.
	void moveFrontToFront(SList<E> &L2) {
		SListPure<E>::moveFrontToFront(L2);
		--m_count; ++L2.m_count;
	}

	//! Moves the first element of this list to the end of list \a L2.
	void moveFrontToBack(SList<E> &L2) {
		SListPure<E>::moveFrontToBack(L2);
		--m_count; ++L2.m_count;
	}

	//! Moves the first element of this list to list \a L2 inserted after \a itBefore.
	/**
	 * \pre \a itBefore points to an element in \a L2.
	 */
	void moveFrontToSucc(SList<E> &L2, SListIterator<E> itBefore) {
		SListPure<E>::moveFrontToSucc(L2,itBefore);
		--m_count; ++L2.m_count;
	}

	//! Appends \a L2 to this list and makes \a L2 empty.
	void conc(SList<E> &L2) {
		SListPure<E>::conc(L2);
		m_count += L2.m_count;
		L2.m_count = 0;
	}

	//! Reverses the order of the list elements.
	void reverse() {
		SListPure<E>::reverse();
	}


	//@}
	/**
	 * @name Searching and sorting
	 * These methods provide searching for values and sorting the list.
	 */
	//@{

	//! Scans the list for the specified element and returns an iterator to the first occurrence in the list, or an invalid iterator if not found.
	SListConstIterator<E> search(const E &e) const {
		return SListPure<E>::search(e);
	}

	//! Scans the list for the specified element and returns an iterator to the first occurrence in the list, or an invalid iterator if not found.
	SListIterator<E> search(const E &e) {
		return SListPure<E>::search(e);
	}

	//! Scans the list for the specified element (using the user-defined comparer) and returns an iterator to the first occurrence in the list, or an invalid iterator if not found.
	template<class COMPARER>
	SListConstIterator<E> search(const E &e, const COMPARER &comp) const {
		return SListPure<E>::search(e, comp);
	}

	//! Scans the list for the specified element (using the user-defined comparer) and returns an iterator to the first occurrence in the list, or an invalid iterator if not found.
	template<class COMPARER>
	SListIterator<E> search(const E &e, const COMPARER &comp) {
		return SListPure<E>::search(e, comp);
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
		SListPure<E>::bucketSort(l,h,f);
	}

	//! Sorts the list using bucket sort.
	void bucketSort(BucketFunc<E> &f) {
		SListPure<E>::bucketSort(f);
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
		SListPure<E>::permute(m_count,rng);
	}

	//@}

	OGDF_NEW_DELETE
}; // class SList




// sorts list L using bucket sort
// computes l and h value
template<class E>
void SListPure<E>::bucketSort(BucketFunc<E> &f)
{
	// if less than two elements, nothing to do
	if (m_head == m_tail) return;

	int l, h;
	l = h = f.getBucket(m_head->m_x);

	SListElement<E> *pX;
	for(pX = m_head->m_next; pX; pX = pX->m_next)
	{
		int i = f.getBucket(pX->m_x);
		if (i < l) l = i;
		if (i > h) h = i;
	}

	bucketSort(l,h,f);
}


// sorts list L using bucket sort
template<class E>
void SListPure<E>::bucketSort(int l, int h, BucketFunc<E> &f)
{
	// if less than two elements, nothing to do
	if (m_head == m_tail) return;

	Array<SListElement<E> *> head(l,h,nullptr), tail(l,h);

	SListElement<E> *pX;
	for (pX = m_head; pX; pX = pX->m_next) {
		int i = f.getBucket(pX->m_x);
		if (head[i])
			tail[i] = (tail[i]->m_next = pX);
		else
			head[i] = tail[i] = pX;
	}

	SListElement<E> *pY = nullptr;
	for (int i = l; i <= h; i++) {
		pX = head[i];
		if (pX) {
			if (pY)
				pY->m_next = pX;
			else
				m_head = pX;
			pY = tail[i];
		}
	}

	m_tail = pY;
	pY->m_next = nullptr;
}


// permutes elements in list randomly; n is the length of the list
template<class E>
template<class RNG>
void SListPure<E>::permute(const int n, RNG &rng)
{
	Array<SListElement<E> *> A(n+1);
	A[n] = nullptr;

	int i = 0;
	SListElement<E> *pX;
	for (pX = m_head; pX; pX = pX->m_next)
		A[i++] = pX;

	A.permute(0,n-1,rng);

	for (i = 0; i < n; i++) {
		A[i]->m_next = A[i+1];
	}

	m_head = A[0];
	m_tail = A[n-1];
}

// prints list to output stream os using delimiter delim
template<class E>
void print(ostream &os, const SListPure<E> &L, char delim = ' ')
{
	SListConstIterator<E> pX = L.begin();
	if (pX.valid()) {
		os << *pX;
		for(++pX; pX.valid(); ++pX)
			os << delim << *pX;
	}
}

// prints list to output stream os using delimiter delim
template<class E>
void print(ostream &os, const SList<E> &L, char delim = ' ')
{
	print(os, L.getSListPure(), delim);
}

// output operator
template<class E>
ostream &operator<<(ostream &os, const SListPure<E> &L)
{
	print(os,L);
	return os;
}

template<class E>
ostream &operator<<(ostream &os, const SList<E> &L)
{
	return operator<<(os,L.getSListPure());
}


// sort array using bucket sort and bucket object f;
// the values of f must be in the interval [min,max]
template<class E>
void bucketSort(Array<E> &a, int min, int max, BucketFunc<E> &f)
{
	if (a.low() >= a.high()) return;

	Array<SListPure<E> > bucket(min,max);

	int i;
	for(i = a.low(); i <= a.high(); ++i)
		bucket[f.getBucket(a[i])].pushBack(a[i]);

	i = a.low();
	for(int j = min; j <= max; ++j) {
		SListConstIterator<E> it = bucket[j].begin();
		for(; it.valid(); ++it)
			a[i++] = *it;
	}
}




} // namespace ogdf
