/*
 * $Revision: 3373 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-05 14:02:17 +0200 (Fr, 05. Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Decralation of GraphElement and GraphList classes.
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

#ifdef _MSC_VER
#pragma once
#endif

#ifndef OGDF_GRAPH_LIST_H
#define OGDF_GRAPH_LIST_H

#include <ogdf/basic/List.h>

namespace ogdf {

class OGDF_EXPORT GraphListBase;

//! The base class for objects used by (hyper)graphs.
/**
 * Such graph objects are maintained in list (see GraphList<T>),
 * and \a GraphElement basically provides a next and previous pointer
 * for these objects.
 */
class OGDF_EXPORT GraphElement {

	friend class Graph;
	friend class GraphListBase;

protected:

	GraphElement *m_next; //!< The successor in the list.
	GraphElement *m_prev; //!< The predecessor in the list.

	OGDF_NEW_DELETE

}; // class GraphElement

//! Base class for GraphElement lists.
class OGDF_EXPORT GraphListBase {

protected:

	int m_size;           //!< The size of the list.
	GraphElement *m_head; //!< Pointer to the first element in the list.
	GraphElement *m_tail; //!< Pointer to the last element in the list.

public:

	//! Constructs an empty list.
	GraphListBase() { m_head = m_tail = 0; m_size = 0; }

	// destruction
	~GraphListBase() { }

	//! Retuns the size of the list.
	int size() const
	{
		return m_size;
	}

	//! Adds element \a pX at the end of the list.
	void pushBack(GraphElement *pX) {
		pX->m_next = 0;
		pX->m_prev = m_tail;
		if (m_head)
			m_tail = m_tail->m_next = pX;
		else
			m_tail = m_head = pX;
		++m_size;
	}

	//! Inserts element \a pX after element \a pY.
	void insertAfter(GraphElement *pX, GraphElement *pY) {
		pX->m_prev = pY;
		GraphElement *pYnext = pX->m_next = pY->m_next;
		pY->m_next = pX;
		if (pYnext) pYnext->m_prev = pX;
		else m_tail = pX;
		++m_size;
	}

	//! Inserts element \a pX before element \a pY.
	void insertBefore(GraphElement *pX, GraphElement *pY) {
		pX->m_next = pY;
		GraphElement *pYprev = pX->m_prev = pY->m_prev;
		pY->m_prev = pX;
		if (pYprev) pYprev->m_next = pX;
		else m_head = pX;
		++m_size;
	}

	//! Removes element \a pX from the list.
	void del(GraphElement *pX) {
		GraphElement *pxPrev = pX->m_prev, *pxNext = pX->m_next;

		if (pxPrev)
			pxPrev->m_next = pxNext;
		else
			m_head = pxNext;
		if (pxNext)
			pxNext->m_prev = pxPrev;
		else
			m_tail = pxPrev;
		m_size--;
	}

	//! Sorts the list according to \a newOrder.
	template<class LIST>
	void sort(const LIST &newOrder) {
		GraphElement *pPred = 0;
		typename LIST::const_iterator it = newOrder.begin();
		if (!it.valid()) return;

		m_head = *it;
		for(; it.valid(); ++it) {
			GraphElement *p = *it;
			if ((p->m_prev = pPred) != 0) pPred->m_next = p;
			pPred = p;
		}

		(m_tail = pPred)->m_next = 0;
	}

	//! Reverses the order of the list elements.
	void reverse() {
		GraphElement *pX = m_head;
		m_head = m_tail;
		m_tail = pX;
		while(pX) {
			GraphElement *pY = pX->m_next;
			pX->m_next = pX->m_prev;
			pX = pX->m_prev = pY;
		}
	}

	//! Exchanges the positions of \a pX and \a pY in the list.
	void swap(GraphElement *pX, GraphElement *pY) {
		if (pX->m_next == pY) {
			pX->m_next = pY->m_next;
			pY->m_prev = pX->m_prev;
			pY->m_next = pX;
			pX->m_prev = pY;

		} else if(pY->m_next == pX) {
			pY->m_next = pX->m_next;
			pX->m_prev = pY->m_prev;
			pX->m_next = pY;
			pY->m_prev = pX;

		} else {
			::swap(pX->m_next,pY->m_next);
			::swap(pX->m_prev,pY->m_prev);
		}

		if(pX->m_prev)
			pX->m_prev->m_next = pX;
		else
			m_head = pX;
		if(pX->m_next)
			pX->m_next->m_prev = pX;
		else
			m_tail = pX;

		if(pY->m_prev)
			pY->m_prev->m_next = pY;
		else
			m_head = pY;
		if(pY->m_next)
			pY->m_next->m_prev = pY;
		else
			m_tail = pY;

		OGDF_ASSERT(consistencyCheck());
	}

	//! Checks consistency of graph list.
	bool consistencyCheck() {
		if (m_head == 0) {
			return (m_tail == 0);

		} else if (m_tail == 0) {
			return false;

		} else {
			if (m_head->m_prev != 0)
				return false;
			if (m_tail->m_next != 0)
				return false;

			GraphElement *pX = m_head;
			for(; pX; pX = pX->m_next) {
				if (pX->m_prev) {
					if (pX->m_prev->m_next != pX)
						return false;
				} else if(pX != m_head)
					return false;

				if (pX->m_next) {
					if (pX->m_next->m_prev != pX)
						return false;
				} else if (pX != m_tail)
					return false;
			}
		}

		return true;
	}

	OGDF_NEW_DELETE

}; // class GraphListBase

//! Lists of graph objects (like nodes, edges, etc.).
/**
 * The template type \a T must be a class derived from GraphElement.
 */
template<class T> class GraphList : protected GraphListBase {

public:

	//! Constructs an empty list.
	GraphList() { }

	// destruction (deletes all elements)
	~GraphList() {
		if (m_head)
			OGDF_ALLOCATOR::deallocateList(sizeof(T), m_head,m_tail);
	}

	//! Returns the size of the list.
	int size () const { return m_size; }

	//! Returns the first element in the list.
	T *begin () const { return (T *)m_head; }

	//! Returns the last element in the list.
	T *rbegin() const { return (T *)m_tail; }

	//! Returns true iff the list is empty.
	bool empty() { return m_head; }

	//! Adds element \a pX at the end of the list.
	void pushBack(T *pX) {
		GraphListBase::pushBack(pX);
	}

	//! Inserts element \a pX after element \a pY.
	void insertAfter(T *pX, T *pY) {
		GraphListBase::insertAfter(pX,pY);
	}

	//! Inserts element \a pX before element \a pY.
	void insertBefore(T *pX, T *pY) {
		GraphListBase::insertBefore(pX,pY);
	}

	//! Moves element \a pX to list \a L and inserts it before or after \a pY.
	void move(T *pX, GraphList<T> &L, T *pY, Direction dir) {
		GraphListBase::del(pX);
		if (dir == after)
			L.insertAfter(pX,pY);
		else
			L.insertBefore(pX,pY);
	}

	//! Moves element \a pX to list \a L and inserts it at the end.
	void move(T *pX, GraphList<T> &L) {
		GraphListBase::del(pX);
		L.pushBack(pX);
	}

	//! Moves element \a pX from its current position to a position after \a pY.
	void moveAfter(T *pX, T *pY){
		GraphListBase::del(pX);
		insertAfter(pX,pY);
	}

	//! Moves element \a pX from its current position to a position before \a pY.
	void moveBefore(T *pX, T *pY){
		GraphListBase::del(pX);
		insertBefore(pX,pY);
	}

	//! Removes element \a pX from the list and deletes it.
	void del(T *pX) {
		GraphListBase::del(pX);
		delete pX;
	}

	//! Only removes element \a pX from the list; does not delete it.
	void delPure(T *pX) {
		GraphListBase::del(pX);
	}

	//! Removes all elements from the list and deletes them.
	void clear() {
		if (m_head) {
			OGDF_ALLOCATOR::deallocateList(sizeof(T),m_head,m_tail);
			m_head = m_tail = 0;
			m_size = 0;
		}
	}

	//! Sorts all elements according to \a newOrder.
	template<class T_LIST>
	void sort(const T_LIST &newOrder) {
		GraphListBase::sort(newOrder);
	}


	//! Reverses the order of the list elements.
	void reverse() {
		GraphListBase::reverse();
	}

	//! Exchanges the positions of \a pX and \a pY in the list.
	void swap(T *pX, T *pY) {
		GraphListBase::swap(pX,pY);
	}


	//! Checks consistency of graph list; returns true if ok.
	bool consistencyCheck() {
		return GraphListBase::consistencyCheck();
	}

	OGDF_NEW_DELETE

}; // class GraphList<T>

} //namespace

#endif

