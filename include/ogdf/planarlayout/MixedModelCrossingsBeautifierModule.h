/** \file
 * \brief Declaration of interface for mixed-model crossings
 * beautifier algorithms
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

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/memory.h>

namespace ogdf {
class GridLayout;
class PlanRep;
template<class E>
class List;

/**
 * \brief The base class for Mixed-Model crossings beautifier algorithms.
 *
 * The class MixedModelCrossingsBeautifierModule is the base class for
 * mixed model bend crossing modules. Such a module transforms an input
 * graph \a G into an output graph \a G' such that crossings of edges don't
 * look weird.
 *
 * <H3>Implementation of Mixed-Model Crossings Beautifier Algorithms</H3>
 *
 * An implementation of a Mixed-Model crossings beautifier module must override
 * the protected method doCall().
 */

class OGDF_EXPORT MixedModelCrossingsBeautifierModule {
public:
	//! Initializes the Mixed-Model crossings beautifier module.
	MixedModelCrossingsBeautifierModule() { }

	// destruction
	virtual ~MixedModelCrossingsBeautifierModule() { }

	/*
	 * \brief Calls the Mixed-Model crossings beautifier module for graph \p PG and grid layout \p gl.
	 *
	 * @param PG is the input graph.
	 * @param gl is the grid layout of \p PG.
	 */
	void call(const PlanRep& PG, GridLayout& gl);

	//! Returns the number of processed crossings.
	int numberOfCrossings() const { return m_nCrossings; }


protected:
	/**
	 * \brief Implements the crossings beautifier module.
	 *
	 * @param PG is the input graph.
	 * @param gl is the grid layout of \p PG.
	 * @param L is the list of crossing nodes.
	 */
	virtual void doCall(const PlanRep& PG, GridLayout& gl, const List<node>& L) = 0;

private:
	int m_nCrossings; //!< the number of processed crossings.

	OGDF_MALLOC_NEW_DELETE
};

//! Dummy implementation of Mixed-Model crossings beautifier.
/**
 * This implementation does no beautification at all and can thus be used
 * for obtaining the original Mixed-Model layout.
 */
class MMDummyCrossingsBeautifier : public MixedModelCrossingsBeautifierModule {
protected:
	//! Dummy implementation.
	virtual void doCall(const PlanRep&, GridLayout&, const List<node>&) override { }
};

}
