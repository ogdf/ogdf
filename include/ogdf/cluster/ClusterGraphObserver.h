/** \file
 * \brief Abstract base class for structures on graphs, that need
 *        to be informed about cluster graph changes.
 *
 * Follows the observer pattern: cluster graphs are observable
 * objects that can inform observers on changes made to their
 * structure.
 *
 * \author Martin Gronemann
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

#include <ogdf/basic/List.h>
#include <ogdf/basic/Observer.h>
#include <ogdf/cluster/ClusterGraph.h>

namespace ogdf {

/**
 * Abstract base class for cluster graph observers.
 *
 * @ingroup graphs
 *
 * If a class needs to keep track of changes in a clustered graph like addition or deletion
 * of clusters, you can derive it from ClusterGraphObserver and override the
 * notification methods clusterDeleted, clusterAdded.
 */
class OGDF_EXPORT ClusterGraphObserver : public Observer<ClusterGraph, ClusterGraphObserver> {
public:
	ClusterGraphObserver() : Observer(nullptr) { }

	explicit ClusterGraphObserver(const ClusterGraph* CG) : Observer(CG) { }

	virtual void clusterDeleted(cluster v) = 0;
	virtual void clusterAdded(cluster v) = 0;

	// virtual void cleared() = 0;

	const ClusterGraph* getGraph() const { return getObserved(); }
};

}
