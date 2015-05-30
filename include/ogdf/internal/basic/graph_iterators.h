/** \file
 * \brief Decralation of graph iterators.
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

#ifndef OGDF_GRAPH_ITERATORS_H
#define OGDF_GRAPH_ITERATORS_H


namespace ogdf {

	namespace internal {

		template<class GraphObjectPtr>
		class GraphIterator {
			GraphObjectPtr m_ptr;

		public:
			GraphIterator() : m_ptr(nullptr) { }
			GraphIterator(GraphObjectPtr ptr) : m_ptr(ptr) { }

			bool operator==(const GraphIterator<GraphObjectPtr> &other) const {
				return (m_ptr == other.m_ptr);
			}

			bool operator!=(const GraphIterator<GraphObjectPtr> &other) const {
				return (m_ptr != other.m_ptr);
			}

			GraphObjectPtr operator*() {
				return m_ptr;
			}

			//! Increment operator (prefix).
			GraphIterator<GraphObjectPtr> &operator++() {
				OGDF_ASSERT(m_ptr != nullptr);
				m_ptr = m_ptr->succ();
				return *this;
			}

			//! Increment operator (postfix).
			GraphIterator<GraphObjectPtr> operator++(int) {
				OGDF_ASSERT(m_ptr != nullptr);
				GraphObjectPtr ptr = m_ptr;
				m_ptr = m_ptr->succ();
				return ptr;
			}

			//! Decrement operator (prefix).
			GraphIterator<GraphObjectPtr> &operator--() {
				OGDF_ASSERT(m_ptr != nullptr);
				m_ptr = m_ptr->pred();
				return *this;
			}

			//! Decrement operator (postfix).
			GraphIterator<GraphObjectPtr> operator--(int) {
				OGDF_ASSERT(m_ptr != nullptr);
				GraphObjectPtr ptr = m_ptr;
				m_ptr = m_ptr->pred();
				return ptr;
			}
		};


		template<class ArrayType> class GraphArrayConstIterator;


		template<class ArrayType>
		class GraphArrayIterator {

			friend class GraphArrayConstIterator<ArrayType>;

			typedef typename ArrayType::key_type    key_type;
			typedef typename ArrayType::value_type  value_type;

			key_type   m_key;
			ArrayType *m_array;

		public:
			GraphArrayIterator() : m_key(nullptr), m_array(nullptr) { }
			GraphArrayIterator(const GraphArrayIterator<ArrayType> &iter) : m_key(iter.m_key), m_array(iter.m_array) { }

			GraphArrayIterator(key_type key, ArrayType *a) : m_key(key), m_array(a) { }

			GraphArrayIterator<ArrayType> &operator=(const GraphArrayIterator<ArrayType> &iter) {
				m_key = iter.m_key;
				m_array = iter.m_array;
				return *this;
			}

			bool operator==(const GraphArrayIterator<ArrayType> &iter) const {
				return m_key == iter.m_key && m_array == iter.m_array;
			}

			bool operator!=(const GraphArrayIterator<ArrayType> &iter) const {
				return !operator==(iter);
			}

			key_type key() const { return m_key; }

			value_type &value() const { return (*m_array)[m_key]; }

			value_type &operator*() const { return (*m_array)[m_key]; }

			//! Increment operator (prefix).
			GraphArrayIterator<ArrayType> &operator++() {
				OGDF_ASSERT(m_key != nullptr);
				//m_key = m_key->succ();
				m_key = ArrayType::findSuccKey(m_key);
				return *this;
			}

			//! Increment operator (postfix).
			GraphArrayIterator<ArrayType> operator++(int) {
				GraphArrayIterator<ArrayType> iter = *this;
				//m_key = m_key->succ();
				m_key = ArrayType::findSuccKey(m_key);
				return iter;
			}

			//! Decrement operator (prefix).
			GraphArrayIterator<ArrayType> &operator--() {
				OGDF_ASSERT(m_key != nullptr);
				//m_key = m_key->pred();
				m_key = ArrayType::findPredKey(m_key);
				return *this;
			}

			//! Decrement operator (postfix).
			GraphArrayIterator<ArrayType> operator--(int) {
				GraphArrayIterator<ArrayType> iter = *this;
				//m_key = m_key->pred();
				m_key = ArrayType::findPredKey(m_key);
				return iter;
			}
		};

		template<class ArrayType>
		class GraphArrayConstIterator {

			typedef typename ArrayType::key_type    key_type;
			typedef typename ArrayType::value_type  value_type;

			key_type         m_key;
			const ArrayType *m_array;

		public:
			GraphArrayConstIterator() : m_key(nullptr), m_array(nullptr) { }
			GraphArrayConstIterator(const GraphArrayConstIterator<ArrayType> &iter) : m_key(iter.m_key), m_array(iter.m_array) { }
			GraphArrayConstIterator(const GraphArrayIterator<ArrayType> &iter) : m_key(iter.m_key), m_array(iter.m_array) { }

			GraphArrayConstIterator(key_type key, const ArrayType *a) : m_key(key), m_array(a) { }

			GraphArrayConstIterator<ArrayType> &operator=(const GraphArrayConstIterator<ArrayType> &iter) {
				m_key = iter.m_key;
				m_array = iter.m_array;
				return *this;
			}

			bool operator==(const GraphArrayConstIterator<ArrayType> &iter) const {
				return m_key == iter.m_key && m_array == iter.m_array;
			}

			bool operator!=(const GraphArrayConstIterator<ArrayType> &iter) const {
				return !operator==(iter);
			}

			key_type key() const { return m_key; }

			const value_type &value() const { return (*m_array)[m_key]; }

			const value_type &operator*() const { return (*m_array)[m_key]; }

			GraphArrayConstIterator<ArrayType> &operator++() {
				OGDF_ASSERT(m_key != nullptr);
				m_key = m_key->succ();
				return *this;
			}

			//! Increment operator (postfix).
			GraphArrayConstIterator<ArrayType> operator++(int) {
				GraphArrayConstIterator<ArrayType> iter = *this;
				m_key = m_key->succ();
				return iter;
			}

			GraphArrayConstIterator<ArrayType> &operator--() {
				OGDF_ASSERT(m_key != nullptr);
				m_key = m_key->pred();
				return *this;
			}

			//! Decrement operator (postfix).
			GraphArrayConstIterator<ArrayType> operator--(int) {
				GraphArrayConstIterator<ArrayType> iter = *this;
				m_key = m_key->pred();
				return iter;
			}
		};

	}
}


#endif
