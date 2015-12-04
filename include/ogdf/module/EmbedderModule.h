/** \file
 * \brief Declaration of interface for embedder for
 * graphs.
 *
 * \author Thorsten Kerkhof (thorsten.kerkhof@udo.edu)
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

#pragma once

#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/Module.h>
#include <ogdf/basic/Timeouter.h>

namespace ogdf {

/**
 * \brief Base class for embedder algorithms.
 *
 * An embedder algorithm computes a planar embedding of a planar
 * graph.
 *
 * \see PlanarizationLayout, PlanarizationGridLayout
 */
class OGDF_EXPORT EmbedderModule : public Module, public Timeouter {
public:
	//! Initializes an embedder module.
	EmbedderModule() { }

	virtual ~EmbedderModule() { }

	/**
	 * \brief Calls the embedder algorithm for graph \a G.
	 * \pre \a G is planar.
	 * \param G is the graph that shall be embedded.
	 * \param adjExternal is set (by the algorithm) to an adjacency entry on the
	 *        external face of \a G.
	 */
	void call(Graph& G, adjEntry& adjExternal) {
#ifdef OGDF_DEBUG
		if(!isPlanar(G)) {
			throw PreconditionViolatedException();
		}
#endif
		doCall(G, adjExternal);
	};

	//! Calls the embedder algorithm for planarized representation \a PG.
	void operator()(Graph& G, adjEntry& adjExternal) { call(G, adjExternal); }

	OGDF_MALLOC_NEW_DELETE
protected:

	/**
	 * \brief Calls the embedder algorithm for graph \a G.
	 * \a G is guaranteed to be planar.
	 * See #call .
	 */
	virtual void doCall(Graph& G, adjEntry& adjExternal) = 0;
};

} // end namespace ogdf
