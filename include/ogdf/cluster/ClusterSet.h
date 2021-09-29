/** \file
 * \brief Declaration and implementation of class ClusterSetSimple,
 * ClusterSetPure and ClusterSet
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
 * In contrast to ClusterSet, a ClusterSetPure does not provide efficient access
 * to the number of clusters stored in the set.
 *
 * \sa ClusterSet
 */
class OGDF_EXPORT ClusterSetPure : public RegisteredSet<ClusterGraph, false> {
	using RS = RegisteredSet<ClusterGraph, false>;

public:
	using RS::RS;

	//! Returns a reference to the list of clusters contained in \a S.
	/**
	 * This list can be used for iterating over all clusters in \a S.
	 */
	const ListPure<cluster>& clusters() const { return RS::elements(); }
};

//! Cluster sets.
/**
 * @ingroup graph-containers
 *
 * A cluster set maintains a subset \a S of the clusters contained in an associated
 * clustered graph. This kind of cluster set provides efficient operations for testing
 * membership, insertion and deletion of elements, and clearing the set.
 *
 * In contrast to ClusterSetPure, a ClusterSet provides efficient access
 * to the number of clusters stored in the set.
 *
 * \sa - ClusterSetPure
 */
class OGDF_EXPORT ClusterSet : public RegisteredSet<ClusterGraph, true> {
	using RS = RegisteredSet<ClusterGraph, true>;

public:
	using RS::RS;

	//! Returns a reference to the list of clusters contained in \a S.
	/**
	 * This list can be used for iterating over all clusters in \a S.
	 */
	const List<cluster>& clusters() const { return RS::elements(); }
};

using ClusterSetSimple = ClusterSetPure;

}
