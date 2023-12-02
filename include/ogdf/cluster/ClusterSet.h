/** \file
 * \brief Declaration and implementation of class ClusterSet
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

#include <ogdf/basic/RegisteredSet.h>
#include <ogdf/cluster/ClusterArray.h>

namespace ogdf {

//! Cluster sets.
/**
 * @ingroup graph-containers
 *
 * A cluster set maintains a subset \a S of the clusters contained in an associated
 * clustered graph. This kind of cluster set provides efficient operations for testing
 * membership, insertion and deletion of elements, and clearing the set.
 *
 * @tparam SupportFastSizeQuery Whether this set supports querying its #size in
 * constant instead of linear time (in the size).
 */
template<bool SupportFastSizeQuery = true>
class ClusterSet : public RegisteredSet<ClusterGraph, SupportFastSizeQuery> {
	using RS = RegisteredSet<ClusterGraph, SupportFastSizeQuery>;

public:
	using RS::RS;

	//! Returns a reference to the list of clusters contained in \a S.
	/**
	 * This list can be used for iterating over all clusters in \a S.
	 */
	const typename RS::list_type& clusters() const { return RS::elements(); }
};

}
