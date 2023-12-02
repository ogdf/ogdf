/** \file
 * \brief Declaration and implementation of EdgeArray class.
 *
 * \author Carsten Gutwenger
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

#include <ogdf/basic/Graph_d.h>
#include <ogdf/basic/RegisteredArray.h>

namespace ogdf {
//! Bucket function for edges.
/**
 * The bucket of an edge is stored in an edge array which is passed
 * by the user at construction; only a pointer is stored to that array.
 */
class OGDF_EXPORT BucketEdgeArray : public BucketFunc<edge> {
	const EdgeArray<int>* m_pEdgeArray; //!< Pointer to edge array.

public:
	//! Constructs a bucket function.
	/**
	 * @param edgeArray contains the buckets for the edges. May not be deleted
	 *        as long as the bucket function is used.
	 */
	explicit BucketEdgeArray(const EdgeArray<int>& edgeArray) : m_pEdgeArray(&edgeArray) { }

	//! Returns bucket of edge \p e.
	int getBucket(const edge& e) override { return (*m_pEdgeArray)[e]; }
};
}
