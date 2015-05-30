/*
 * $Revision: 3951 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2014-03-03 13:57:46 +0100 (Mon, 03 Mar 2014) $
 ***************************************************************/

/** \file
 * \brief declaration and implementation of class FaceSetSimple,
 * FaceSetPure and FaceSet
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

#ifndef OGDF_FACE_SET_H
#define OGDF_FACE_SET_H


#include <ogdf/basic/FaceArray.h>
#include <ogdf/basic/SList.h>



namespace ogdf {


//! Simple face sets.
/**
 * A face set maintains a subset \a S of the faces contained in an associated
 * combinatorial embedding.
 * This kind of face set only provides efficient operation for testing membership,
 * insertion, and clearing the set.
 *
 * \sa
 *   - FaceSet, FaceSetPure
 *   - NodeSet, NodeSetPure, NodeSetSimple
 */
class OGDF_EXPORT FaceSetSimple {
public:
	//! Creates an empty face set associated with combinatorial embedding \a E.
	FaceSetSimple(const CombinatorialEmbedding &E) : m_isContained(E,false) { }

	// destructor
	~FaceSetSimple() { }

	//! Inserts face \a f into set \a S.
	/**
	 * This operation has constant runtime.
	 *
	 * \pre \a f is a face in the associated combinatorial embedding.
	 */
	void insert(face f) {
		OGDF_ASSERT(f->embeddingOf() == m_isContained.embeddingOf());
		bool &isContained = m_isContained[f];
		if (isContained == false) {
			isContained = true;
			m_faces.pushFront(f);
		}
	}


	//! Removes all faces from \a S.
	/**
	 * After this operation, the \a S is empty and still associated with the same combinatorial embedding.
	 * The runtime of this operations is O(k), where k is the number of faces in \a S before
	 * this operation.
	 */
	void clear() {
		SListIterator<face> it;
		for(it = m_faces.begin(); it.valid(); ++it) {
			m_isContained[*it] = false;
		}
		m_faces.clear();
	}


	//! Returns true if face \a f is contained in \a S, false otherwise.
	/**
	 * This operation has constant runtime.
	 *
	 * \pre \a f is a face in the associated combinatorial embedding.
	 */
	bool isMember(face f) const {
		OGDF_ASSERT(f->embeddingOf() == m_isContained.embeddingOf());
		return m_isContained[f];
	}

	//! Returns a reference to the list of faces contained in \a S.
	/**
	 * This list can be used for iterating over all faces in \a S.
	 */
	const SListPure<face> &faces() const {
		return m_faces;
	}

private:
	//! m_isContained[f] is true iff \a f is contained in \a S.
	FaceArray<bool> m_isContained;

	//! The list of faces contained in \a S.
	SListPure<face> m_faces;
};



//! Face sets.
/**
 * A face set maintains a subset \a S of the faces contained in an associated
 * combinatorial embedding. This kind of face set provides efficient operations for testing
 * membership, insertion and deletion of elements, and clearing the set.
 *
 * In contrast to FaceSet, a FaceSetPure does not provide efficient access
 * to the number of faces stored in the set.
 *
 * \sa
 *   - FaceSet, FaceSetSimple
 *   - NodeSet, NodeSetPure, NodeSetSimple
 */
class OGDF_EXPORT FaceSetPure {
public:
	//! Creates an empty node set associated with combinatorial embedding \a E.
	FaceSetPure(const CombinatorialEmbedding &E) : m_it(E,ListIterator<face>()) { }

	// destructor
	~FaceSetPure() { }

	//! Inserts face \a f into \a S.
	/**
	 * This operation has constant runtime.
	 *
	 * \pre \a f is a face in the associated combinatorial embedding.
	 */
	void insert(face f) {
		OGDF_ASSERT(f->embeddingOf() == m_it.embeddingOf());
		ListIterator<face> &itF = m_it[f];
		if (!itF.valid())
			itF = m_faces.pushBack(f);
	}

	//! Removes face \a f from \a S.
	/**
	 * This operation has constant runtime.
	 *
	 * \pre \a f is a face in the associated combinatorial embedding.
	 */
	void remove(face f) {
		OGDF_ASSERT(f->embeddingOf() == m_it.embeddingOf());
		ListIterator<face> &itF = m_it[f];
		if (itF.valid()) {
			m_faces.del(itF);
			itF = ListIterator<face>();
		}
	}


	//! Removes all faces from \a S.
	/**
	 * After this operation, \a S is empty and still associated with the same combinatorial embedding.
	 * The runtime of this operations is O(k), where k is the number of faces in \a S
	 * before this operation.
	 */
	void clear() {
		ListIterator<face> it;
		for(it = m_faces.begin(); it.valid(); ++it) {
			m_it[*it] = ListIterator<face>();
		}
		m_faces.clear();
	}


	//! Returns true if face \a f is contained in \a S, false otherwise.
	/**
	 * This operation has constant runtime.
	 *
	 * \pre \a f is a face in the associated combinatorial embedding.
	 */
	bool isMember(face f) const {
		OGDF_ASSERT(f->embeddingOf() == m_it.embeddingOf());
		return m_it[f].valid();
	}

	//! Returns a reference to the list of faces contained in \a S.
	/**
	 * This list can be used for iterating over all faces in \a S.
	 */
	const ListPure<face> &faces() const {
		return m_faces;
	}

private:
	//! m_it[f] contains the list iterator pointing to \a f if \a f is contained in S,
	//! an invalid list iterator otherwise.
	FaceArray<ListIterator<face> > m_it;

	//! The list of faces contained in \a S.
	ListPure<face> m_faces;
};



//! Face sets.
/**
 * A face set maintains a subset \a S of the faces contained in an associated
 * combinatorial embedding. This kind of face set provides efficient operations for testing
 * membership, insertion and deletion of elements, and clearing the set.
 *
 * In contrast to FaceSetPure, a FaceSet provides efficient access
 * to the number of elements stored in the set.
 *
 * \sa
 *   - FaceSetPure, FaceSetSimple
 *   - NodeSet, NodeSetPure, NodeSetSimple
 */
class OGDF_EXPORT FaceSet {
public:
	//! Creates an empty node set associated with combinatorial embedding \a E.
	FaceSet(const CombinatorialEmbedding &E) : m_it(E,ListIterator<face>()) { }

	// destructor
	~FaceSet() { }

	//! Inserts face \a f into \a S.
	/**
	 * This operation has constant runtime.
	 *
	 * \pre \a f is a face in the associated combinatorial embedding.
	 */
	void insert(face f) {
		OGDF_ASSERT(f->embeddingOf() == m_it.embeddingOf());
		ListIterator<face> &itF = m_it[f];
		if (!itF.valid())
			itF = m_faces.pushBack(f);
	}

	//! Removes face \a f from \a S.
	/**
	 * This operation has constant runtime.
	 *
	 * \pre \a f is a face in the associated combinatorial embedding.
	 */
	void remove(face f) {
		OGDF_ASSERT(f->embeddingOf() == m_it.embeddingOf());
		ListIterator<face> &itF = m_it[f];
		if (itF.valid()) {
			m_faces.del(itF);
			itF = ListIterator<face>();
		}
	}


	//! Removes all faces from \a S.
	/**
	 * After this operation, \a S is empty and still associated with the same combinatorial embedding.
	 * The runtime of this operations is O(k), where k is the number of faces in \a S
	 * before this operation.
	 */
	void clear() {
		ListIterator<face> it;
		for(it = m_faces.begin(); it.valid(); ++it) {
			m_it[*it] = ListIterator<face>();
		}
		m_faces.clear();
	}


	//! Returns true if face \a f is contained in \a S, false otherwise.
	/**
	 * This operation has constant runtime.
	 *
	 * \pre \a f is a face in the associated combinatorial embedding.
	 */
	bool isMember(face f) const {
		OGDF_ASSERT(f->embeddingOf() == m_it.embeddingOf());
		return m_it[f].valid();
	}

	//! Returns the size of \a S.
	/**
	 * This operation has constant runtime.
	 */
	int size() const {
		return m_faces.size();
	}

	//! Returns a reference to the list of faces contained in \a S.
	/**
	 * This list can be used for iterating over all faces in \a S.
	 */
	const List<face> &faces() const {
		return m_faces;
	}

private:
	//! m_it[f] contains the list iterator pointing to \a f if \a f is contained in S,
	//! an invalid list iterator otherwise.
	FaceArray<ListIterator<face> > m_it;

	//! The list of faces contained in \a S.
	List<face> m_faces;
};


} // end namespace ogdf


#endif

