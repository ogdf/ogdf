/** \file
 * \brief Abstract base class for observers on hypergraphs, that
 *        need to be informed about hypergraph changes
 *        (e.g. associated graph edge standard representation).
 *
 * Follows the observer pattern: hypergraphs are observable
 * objects that can inform observers on changes made to their
 * structure.
 *
 * Based on the original GraphObserver class written by Karsten Klein.
 *
 * \author Ondrej Moris
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
#include <ogdf/hypergraph/Hypergraph.h>

namespace ogdf {

// HypergraphObserver
class OGDF_EXPORT HypergraphObserver : public Observer<Hypergraph, HypergraphObserver> {
public:
	//! Constructor.
	HypergraphObserver() = default;

	//! Constructor assigning \p pH hypergraph to the observer.
	explicit HypergraphObserver(const Hypergraph* pH) { reregister(pH); }

	//! Associates an observer instance with hypergraph \p pH
	OGDF_DEPRECATED("reregister() should be used instead.")

	void init(const Hypergraph* pH) { reregister(pH); }

	//! Called by an observed hypergraph when a hypernode is deleted.
	virtual void hypernodeDeleted(hypernode v) = 0;

	//! Called by an observed hypergraph when a hypernode is added.
	virtual void hypernodeAdded(hypernode v) = 0;

	//! Called by an observed hypergraph when a hyperedge is deleted.
	virtual void hyperedgeDeleted(hyperedge e) = 0;

	//! Called by an observed hypergraph when a hyperedge is added.
	virtual void hyperedgeAdded(hyperedge e) = 0;

	//! Called by the observed hypergraph when it is cleared.
	virtual void cleared() = 0;

	//! Returns the observer hypergraph.
	const Hypergraph* hypergraph() const { return getObserved(); }
};

}
