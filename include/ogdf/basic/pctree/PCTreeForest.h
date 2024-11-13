/** \file
 * \brief PCTreeForests contain multiple PCTrees that can be merged with each other.
 *
 * \author Simon D. Fink <ogdf@niko.fink.bayern>
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

#include <ogdf/basic/DisjointSets.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/pctree/PCRegistry.h>

#include <cstddef>
#include <limits>
#include <vector>

namespace ogdf::pc_tree {
class PCNode;
class PCTree;
} // namespace ogdf::pc_tree

#define OGDF_PCTREE_REUSE_NODES

namespace ogdf::pc_tree {
using UnionFindIndex = size_t;

const UnionFindIndex UNIONFINDINDEX_EMPTY = std::numeric_limits<UnionFindIndex>::max();

/**
 * Multiple PCTrees can be created within the same PCTreeForest, which allows merging the trees later on by making one
 * a child of another. This is extensively used during planarity testing.
 * @sa PCTree:insertTree()
 */
class OGDF_EXPORT PCTreeForest {
	friend class PCNode;
	friend class PCTree;
	friend class PCTreeRegistry;

private:
	std::vector<PCTree*> m_trees;
	std::vector<PCNode*> m_cNodes;
	ogdf::DisjointSets<> m_parents {1 << 8};
	int m_nextNodeId = 0;
	size_t m_timestamp = 0;
	PCTreeRegistry m_nodeArrayRegistry;
	bool m_autodelete;

#ifdef OGDF_PCTREE_REUSE_NODES
	PCNode* m_reusableNodes = nullptr;
#endif

public:
	/**
	 * @param autodelete whether the trees created by makeTree() should be deleted automatically
	 *   on destruction of this forrest. Note that this does not affect PCTrees directly created by
	 *   calling PCTree::PCTree(PCTreeForest*).
	 */
	PCTreeForest(bool autodelete = true) : m_nodeArrayRegistry(this), m_autodelete(autodelete) {};

	virtual ~PCTreeForest();

	//! Create a new tree that may be automatically deleted when this forest is deleted.
	PCTree* makeTree(void);

	//! Delete all trees created by makeTree().
	void clear(void);

	operator const PCTreeRegistry&() const { return m_nodeArrayRegistry; }
};
}
