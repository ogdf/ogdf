/*
 * $Revision: 3533 $
 *
 * last checkin:
 *   $Author: beyer $
 *   $Date: 2013-06-03 18:22:41 +0200 (Mo, 03. Jun 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration and implementation of ArrayBuffer class.
 *
 * \author Markus Chimani
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

#ifndef OGDF_ARRAY_BUFFER_H
#define OGDF_ARRAY_BUFFER_H

#include <ogdf/basic/Array.h>
#include <cstring>


namespace ogdf {

//! An array that keeps track of the number of inserted elements; also usable as an efficient stack.
/**
 * This is a (by default automatically growable) array (with some initial size \a s) which starts out being empty. Using
 * stack functions you can put elements into and out of it. The initial array size is automatically
 * expanded if neccessary (unless growing is forbidden), but never automatically shrunken. You may also access the elements it
 * contains using the []-operator. The valid indices are 0..(\a s - 1).
 *
 * @tparam E     denotes the element type.
 * @tparam INDEX denotes the index type. The index type must be chosen such that it can
 *               express the whole index range of the array instance, as well as its size.
 *               The default index type is \c int, other possible types are \c short and
 *               <code>long long</code> (on 64-bit systems).
 */
template<class E, class INDEX = int>
class ArrayBuffer : private Array<E, INDEX> {
	INDEX num; //!< The number of elements in the buffer
	bool growable;
public:
	//! Creates an empty array buffer, without initial memory allocation.
	ArrayBuffer() : Array<E,INDEX>(), num(0), growable(true) {}

	//! Creates an empty array buffer, allocating memory for up to \a size elements; you may specify that the array should not grow automatically.
	explicit ArrayBuffer(INDEX size, bool autogrow = true) : Array<E,INDEX>(size), num(0), growable(autogrow) {}

	//! Creates an array buffer, initialized by the given array; you may specify that the array should not grow.
	explicit ArrayBuffer(const Array<E,INDEX>& source, bool autogrow = true) : Array<E,INDEX>(source), num(0), growable(autogrow) {}

	//! Creates an array buffer that is a copy of \a buffer.
	ArrayBuffer(const ArrayBuffer<E,INDEX> &buffer) : Array<E,INDEX>(buffer), num(buffer.num), growable(buffer.growable) { }

	//! Reinitializes the array, clearing it, and without initial memory allocation.
	void init() { Array<E,INDEX>::init(); }
	//! Reinitializes the array, clearing it, and allocating memory for up to \a size elements.
	void init(INDEX size) { Array<E,INDEX>::init(size); }

	//! Clears the buffer
	void clear() { num = 0; }

	//! Returns the newest element of the buffer.
	const E &top() const { OGDF_ASSERT(num>0); return Array<E,INDEX>::operator[](num-1); }
	//! Returns the newest element of the buffer.
	E &top() { OGDF_ASSERT(num>0); return Array<E,INDEX>::operator[](num-1); }

	//! Puts a new element in the buffer.
	void push(E e) {
		if (num == Array<E,INDEX>::size()) {
			if (growable)
				Array<E,INDEX>::grow(max(num,1)); // double the size
			else
				OGDF_THROW_PARAM(PreconditionViolatedException, pvcFull);
		}
		Array<E,INDEX>::operator[](num++) = e;
	}

	//! Removes the newest element from the buffer.
	void pop() { OGDF_ASSERT(num>0); --num; }
	//! Removes the newest element from the buffer and returns it.
	E popRet() { OGDF_ASSERT(num>0); return Array<E,INDEX>::operator[](--num); }

	//! Returns true if the buffer is empty, false otherwise.
	bool empty() const { return !num; }

	//! Returns true iff the buffer is non-growable and filled.
	bool full() const { return (!growable) && (num == Array<E,INDEX>::size()); }

	//! Returns number of elements in the buffer.
	INDEX size() const { return num; }

	//! Returns the current capacity of the datastructure. Note that this value is rather irrelevant if the array is growable.
	INDEX capacity() const { return Array<E,INDEX>::size(); }

	//! Returns whether the buffer will automatically expand if the initial size is insufficient
	bool isGrowable() const { return growable; }

	//! Sets the flag whether the buffer will automatically expand if the initial size is insufficient
	void setGrowable(bool _growable) { growable = _growable; }

	//! Returns a pointer to the first element.
	E *begin() { return Array<E, INDEX>::begin(); }

	//! Returns a pointer to the first element.
	const E *begin() const { return Array<E, INDEX>::begin(); }

	//! Returns a pointer to one past the last element.
	E *end() { return Array<E, INDEX>::begin()+num; }

	//! Returns a pointer to one past the last element.
	const E *end() const { return Array<E, INDEX>::begin()+num; }

	//! Returns a pointer to the last element.
	E *rbegin() { return Array<E, INDEX>::begin()+(num-1); }

	//! Returns a pointer to the last element.
	const E *rbegin() const { return Array<E, INDEX>::begin()+(num-1); }

	//! Returns a pointer to one before the first element.
	E *rend() { return Array<E, INDEX>::rend(); }

	//! Returns a pointer to one before the first element.
	const E *rend() const { return Array<E, INDEX>::rend(); }

	//! Returns a reference to the element at position \a i.
	const E &operator[](INDEX i) const {
		OGDF_ASSERT(0 <= i && i < num)
		return Array<E,INDEX>::operator[](i);
	}
	//! Returns a reference to the element at position \a i.
	E &operator[](INDEX i) {
		OGDF_ASSERT(0 <= i && i < num)
		return Array<E,INDEX>::operator[](i);
	}

	//! Assignment operator.
	ArrayBuffer<E,INDEX> &operator=(const ArrayBuffer<E,INDEX> &buffer) {
		Array<E,INDEX>::operator=(buffer);
		num      = buffer.num;
		growable = buffer.growable;
		return *this;
	}

	//! Generates a compact copy holding the current elements.
	/**
	 * Creates a copy of the ArrayBuffer and stores it into
	 * the given Array \a A.
	 * \a A has exactly the neccessary size to hold all
	 * elements in the buffer.
	 *
	 * This method uses an elementwise operator=.
	 * If you need a bitcopy of the buffer, use compactMemcpy()
	 * instead; if you need a traditional array copy (using the Array's
	 * copy-constructor) use compactCpycon() instead.
	 */
	void compactCopy(Array<E,INDEX>& A2) const {
		OGDF_ASSERT(this != &A2);
		if(num) {
			A2.init(num);
			for(INDEX i = num; i-->0;)
				A2[i] = (*this)[i];
		} else
			A2.init(0);
	}

	//! Generates a compact copy holding the current elements.
	/**
	 * Creates a copy of the ArrayBuffer and stores it into
	 * the given Array \a A.
	 * \a A has exactly the neccessary size to hold all
	 * elements in the buffer
	 *
	 * This method uses the Array's copy constructur. If you
	 * need a bitcopy of the buffer, use compactMemcpy()
	 * instead; if you need a elementwise operator=-copy, use
	 * compactCopy() instead.
	 */
	void compactCpycon(Array<E,INDEX>& A2) const {
		OGDF_ASSERT(this != &A2);
		if(num) {
			INDEX tmp = Array<E,INDEX>::m_high; // thank god i'm a friend of Array
			Array<E,INDEX>::m_high = num-1; // fake smaller size
			A2.copy(*this); // copy
			Array<E,INDEX>::m_high = tmp;
		} else
			A2.init(0);
	}

	//! Generates a compact copy holding the current elements.
	/**
	 * Creates a copy of the ArrayBuffer and stores it into
	 * the given Array \a A.
	 * \a A has exactly the neccessary size to hold all
	 * elements in the buffer.
	 *
	 * This method uses memcpy. If you need a traditional
	 * arraycopy using a copy constructur, use compactCoycon()
	 * instead; if you neeed a elementwise operator=-copy, use
	 * compactCopy() instead.
	 */
	void compactMemcpy(Array<E,INDEX>& A2) const {
		OGDF_ASSERT(this != &A2);
		if(num) {
			A2.init(num);
			memcpy(A2.m_pStart,this->m_pStart,sizeof(E)*num);
		} else
			A2.init(0);
	}

	//! Performs a linear search for element \a x.
	/**
	 * Warning: linear running time!
	 * Note that the linear search runs from back to front.
	 * \return the index of the found element, and low()-1 if not found.
	 */
	INDEX linearSearch (const E& x) const {
		INDEX i;
		for(i = num; i-->0;)
			if(x == Array<E,INDEX>::m_vpStart[i]) break;
		return i;
	}

	//! Performs a linear search for element \a x with comparer \a comp.
	/**
	 * Warning: linear running time!
	 * Note that the linear search runs from back to front.
	 * \return the index of the found element, and low()-1 if not found.
	 */
	template<class COMPARER>
	INDEX linearSearch (const E& x, const COMPARER &comp) const {
		INDEX i;
		for(i = num; i-->0;)
			if(comp.equal(x, Array<E,INDEX>::m_vpStart[i])) break;
		return i;
	}

	//! Sorts buffer using Quicksort.
	inline void quicksort() {
		Array<E,INDEX>::quicksort(0,num-1,StdComparer<E>());
	}

	//! Sorts buffer using Quicksort and a user-defined comparer \a comp.
	/**
	 * @param comp is a user-defined comparer; \a C must be a class providing a \a less(x,y) method.
	 */
	template<class COMPARER>
	inline void quicksort(const COMPARER &comp) {
		Array<E,INDEX>::quicksort(0,num-1,comp);
	}

	//! Performs a binary search for element \a x.
	/**
	 * \pre The buffer must be sorted!
	 * \return the index of the found element, and low()-1 if not found.
	 */
	inline INDEX binarySearch (const E& e) const {
		return Array<E,INDEX>::binarySearch(0, num-1, e, StdComparer<E>());
	}

	//! Performs a binary search for element \a x with comparer \a comp.
	/**
	 * \pre The buffer must be sorted according to \a comp!
	 * \return the index of the found element, and low()-1 if not found.
	 */
	template<class COMPARER>
	inline INDEX binarySearch(const E& e, const COMPARER &comp) const {
		return Array<E,INDEX>::binarySearch(0, num-1, e, comp);
	}

	//! Randomly permutes the array.
	void permute() {
		Array<E,INDEX>::permute(0, num-1);
	}

	//! Removes the components listed in the buffer \a ind by shifting the remaining components to the left.
	/**
	 * The values stored in \a ind have to be upward sorted.
	 * Memory management of the removed components must be
	 * carefully implemented by the user of this function to avoid
	 * memory leaks.
	 *
	 * If this function is compiled with <tt>OGDF_DEBUG</tt>
	 * then it is checked if each value of \a ind is in the
	 * range 0,..., \a number()-1.
	 *
	 * \param ind The numbers of the components being removed.
	 */
	void leftShift(ArrayBuffer<INDEX, INDEX> &ind) {
		const INDEX nInd = ind.size();
		if (nInd == 0) return;

		//! shift all items up to the last element of \a ind to the left
	#ifdef OGDF_DEBUG
		if(ind[0] < 0 || ind[0] >= num)
			OGDF_THROW_PARAM(AlgorithmFailureException, afcIndexOutOfBounds);
	#endif

		INDEX j, current = ind[0];
		for (INDEX i = 0; i < nInd - 1; i++) {
	#ifdef OGDF_DEBUG
			if(ind[i+1] < 0 || ind[i+1] >= num)
				OGDF_THROW_PARAM(AlgorithmFailureException, afcIndexOutOfBounds);
	#endif

			const INDEX last = ind[i+1];
			for(j = ind[i]+1; j < last; j++)
				operator[](current++) = operator[](j);
		}

		//! copy the rest of the buffer
		for (j = ind[nInd - 1]+1; j < size(); j++)
			operator[](current++) = operator[](j);

		num -= nInd;
	}

	//! Changes the capacity of the buffer (independent whether the buffer is growable of not).
	/**
	 * If the new capacity if smaller that the currently stored elements, only the first elements (as many as fit) are
	 * retained in the buffer. The user is responsible that no memory leaks occur.
	 */
	void setCapacity(INDEX newCapacity) {
		Array<E, INDEX>::resize(newCapacity);
	}

	OGDF_NEW_DELETE
};


// prints array a to output stream os using delimiter delim
template<class E, class INDEX>
void print(ostream &os, const ArrayBuffer<E,INDEX> &a, char delim = ' ')
{
	for (int i = 0; i < a.size(); i++) {
		if (i > 0) os << delim;
		os << a[i];
	}
}


// output operator
template<class E, class INDEX>
ostream &operator<<(ostream &os, const ogdf::ArrayBuffer<E,INDEX> &a)
{
	print(os,a);
	return os;
}

} // end namespace ogdf

#endif
