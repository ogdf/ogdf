/** \file
 * \brief Declaration of class SeparatorLiptonTarjanFC.
 *
 * \author Thomas Klein
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

#include <ogdf/graphalg/PlanarSeparatorModule.h>

#include <set>

namespace ogdf {

namespace planar_separators {

/**
 * Custom BFS Tree.
 */
class OGDF_EXPORT BFSTreeFC : public ArrayBFSTree {
public:
	/**
	 * Constructor.
	 *
	 * @param G the graph
	 * @param rootNode the node at which the tree should be rooted
	 */
	BFSTreeFC(GraphCopy& G, node rootNode) : ArrayBFSTree(G, rootNode) { construct(); }

	void construct();
};

/**
 * Triangulating BFS tree that operates on a non-triangulated graph and
 * constructs the triangulation together with the BFS, which should lead to broader trees.
 */
class OGDF_EXPORT TriangulatingBFSTree : public ArrayBFSTree {
public:
	/**
	 * Constructor.
	 *
	 * @param G the graph
	 * @param rootNode the node at which the tree should be rooted
	 */
	TriangulatingBFSTree(GraphCopy& G, node rootNode) : ArrayBFSTree(G, rootNode) { construct(); }

	void construct();

	void visit(node v, node parent, adjEntry adj, SListPure<node>& bfs);
};

}

//! Computes planar separators using Fundamental Cycles
/**
 * Computes planar separators using only the Fundamental Cycle Lemma (ie. the
 * second half of Lipton Tarjan).
 *
 * @ingroup ga-plansep
 */
class OGDF_EXPORT SeparatorLiptonTarjanFC : public PlanarSeparatorModule {
public:
	/**
	 * Constructor.
	 *
	 * @param useTriBFS whether to use triangulating BFS or not
	 */
	SeparatorLiptonTarjanFC(bool useTriBFS = false) : useTriangulatingBFS {useTriBFS} { }

	virtual double getMaxSeparatorSize(int n) const override { return -1; }

protected:
	virtual bool doSeparate(const Graph& G, List<node>& separator, List<node>& first,
			List<node>& second) override;

	bool useTriangulatingBFS;
	std::shared_ptr<ArrayBFSTree> tree;

	virtual std::string getSpecificName() const override {
		std::string name = "LTFC";
		if (useTriangulatingBFS) {
			name += "-triBFS";
		}
		return name;
	}

	/**
	 * Finds a cycle that works as a separator.
	 *
	 * @param separator the separator nodes
	 * @param first the first component
	 * @param second the second component
	 * @return true on success
	 */
	virtual bool findCycle(List<node>& separator, List<node>& first, List<node>& second);

	/**
	 * Randomly selects the initial edge for the first cycle.
	 *
	 * @return the first edge
	 */
	edge chooseEdge() const;
};

}
