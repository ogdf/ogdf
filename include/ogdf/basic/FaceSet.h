/** \file
 * \brief Declaration and implementation of ogdf::FaceSet.
 *
 * \author Carsten Gutwenger, Tilo Wiedera
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

#include <ogdf/basic/FaceArray.h>
#include <ogdf/basic/RegisteredSet.h>

namespace ogdf {

//! Face sets.
/**
 * @ingroup graph-containers
 *
 * Maintains a subset of faces contained in an associated combinatorial embedding.
 * Provides efficient operations for testing membership,
 * iteration, insertion and deletion of elements, as well as clearing the set.
 *
 * \tparam SupportFastSizeQuery Whether this set supports querying it's #size in
 * constant instead of linear time (in the size).
 *
 * \sa NodeSet
 */
template<bool SupportFastSizeQuery = true>
class FaceSet : public RegisteredSet<ConstCombinatorialEmbedding, SupportFastSizeQuery> {
	using RS = RegisteredSet<ConstCombinatorialEmbedding, SupportFastSizeQuery>;

public:
	using RS::RS;

	//! Returns a reference to the list of faces contained in this set.
	const typename RS::list_type& faces() const { return RS::elements(); }

	//! Returns the associated combinatorial embedding
	const ConstCombinatorialEmbedding& embeddingOf() const {
		OGDF_ASSERT(RS::registeredAt());
		return *RS::registeredAt();
	}
};

}
