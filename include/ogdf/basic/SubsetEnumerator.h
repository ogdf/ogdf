/*
 * $Revision: 3624 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-07-02 13:57:06 +0200 (Di, 02. Jul 2013) $
 ***************************************************************/

/** \file
 * \brief A class that allows to enumerate k-subsets of lists.
 *
 * \author Stephan Beyer
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

#ifndef OGDF_SUBSET_ENUMERATOR_H
#define OGDF_SUBSET_ENUMERATOR_H

#include <ogdf/basic/List.h>

namespace ogdf {

//! Enumerator for k-subsets of a given type.
/**
 * <H3>Usage examples</H3>
 * <ul>
 *   <li> Enumerate all subsets of edges with cardinality 3:
 *    \code
 *     List<edge> edges;
 *
 *     do_something_eg_fill_edges();
 *
 *     SubsetEnumerator<edge> edgeSubset(edges);
 *
 *     for (subset.begin(3); subset.valid(); subset.next()) {
 *       do_something_with(subset[0], subset[1], subset[2]);
 *     }
 *    \endcode
 *   </li>
 *   <li> Enumerate all subsets of edges:
 *    \code
 *     SubsetEnumerator<edge> edgeSubset(edges);
 *
 *     for (int k = 1; k <= edges.size(); ++k) {
 *       for (subset.begin(3); subset.valid(); subset.next()) {
 *         for (int i = 0; i < subset.size(); ++i) {
 *           do_something_with(subset[i]);
 *         }
 *         do_stuff();
 *       }
 *     }
 *    \endcode
 *   </li>
 * </ul>
 *
 * Please note that the internal data structures of SubsetEnumerator do not use references of the type T.
 * Hence, T should either be a simple type or a pointer to a complex type (which is also only sane for Lists, too).
 * Otherwise the data structure will slow down due to extensive implicit copying.
 */
template<typename T>
class SubsetEnumerator {
protected:
	const List<T> &m_set;
	bool m_valid;
	Array<T> m_subset;
	Array<int> m_index;

public:
	//! \brief Constructor.
	//! @param set The list of elements we want to enumerate subsets for.
	SubsetEnumerator(const List<T> &set)
	  : m_set(set)
	  , m_valid(false)
	  , m_subset(set.size())
	{
		int i = 0;
		forall_listiterators(T, it, m_set) {
			m_subset[i++] = *it;
		}
	}

	//! Set the SubsetEnumerator to the first subset of given cardinality.
	void begin(int card)
	{
		if (card > 0) {
			m_index.init(card);
			for (int i = 0; i < card; ++i) {
				m_index[i] = i;
			}
			m_valid = true;
		}
	}

	//! Return the cardinality of the subset.
	int size() const
	{
		return m_index.size();
	}

	//! Is the current subset valid? If not, all subsets have already been enumerated.
	bool valid() const
	{
		return m_valid;
	}

	//! Check in O(subset cardinality) whether the argument is a member of the subset.
	bool hasMember(const T &element) const
	{
		for (int i = 0; i < m_index.size(); ++i) {
			if (element == m_subset[m_index[i]]) {
				return true;
			}
		}
		return false;
	}

	//! Get element of subset by index (starting from 0).
	T operator[](int i) const
	{
		return m_subset[m_index[i]];
	}

	//! Obtain the next subset if possible. The result should be checked using the valid() method.
	void next()
	{
		if (m_valid) {
			const int n = m_subset.size();
			const int t = m_index.size();
			int i;
			for (i = t - 1; m_index[i] == i + n - t; --i) {
				if (i == 0) {
					m_valid = false; // the last subset has been found
					return;
				}
			}
			for (++m_index[i]; i < t - 1; ++i) {
				m_index[i + 1] = m_index[i] + 1;
			}
		}
	}
};

}
#endif // OGDF_SUBSET_ENUMERATOR_H
