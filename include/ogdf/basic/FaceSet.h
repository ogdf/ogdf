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

#include <ogdf/basic/CombinatorialEmbedding.h>
#include <ogdf/basic/RegisteredSet.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/internal/copy_move.h>

namespace ogdf {

//! Face sets.
/**
 * @ingroup graph-containers
 *
 * Maintains a subset of faces contained in an associated combinatorial embedding.
 * Provides efficient operations for testing membership,
 * iteration, insertion and deletion of elements, as well as clearing the set.
 *
 * \sa NodeSet
 */
class OGDF_EXPORT FaceSet : public RegisteredSet<ConstCombinatorialEmbedding> {
public:
	using RegisteredSet::RegisteredSet;
	OGDF_DEFAULT_COPY(FaceSet);

	//! Returns a reference to the list of faces contained in this set.
	const list_type& faces() const { return elements(); }

	//! Returns the associated combinatorial embedding
	const ConstCombinatorialEmbedding& embeddingOf() const {
		OGDF_ASSERT(RegisteredSet::registeredAt());
		return *registeredAt();
	}
};

}
