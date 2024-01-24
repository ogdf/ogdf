/** \file
 * \brief // TODO DESCRIBE WHAT IS IMPLEMENTED
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
#include <ogdf/basic/pctree/PCEnum.h>
#include <ogdf/basic/pctree/PCRegistry.h>

#include <cstdint>
#include <vector>

#define OGDF_PCTREE_REUSE_NODES

namespace pc_tree {
using UnionFindIndex = int;

const int UNIONFINDINDEX_EMPTY = -1;

class PCTreeForest {
	friend class PCNode;
	friend class PCTree;

	template<class Key>
	friend class PCTreeRegistry;

private:
	std::vector<PCTree*> trees;
	std::vector<PCNode*> cNodes;
	ogdf::DisjointSets<> parents {1 << 8};
	int nextNodeId = 0;
	int timestamp = 0;
	PCTreeRegistry<PCNode*> nodeArrayRegistry;
	bool autodelete;

#ifdef OGDF_PCTREE_REUSE_NODES
	// TODO: also reuse PCTrees?
	PCNode* reusableNodes = nullptr;
#endif

public:
	PCTreeForest(bool p_autodelete = true) : nodeArrayRegistry(this), autodelete(p_autodelete) {};

	virtual ~PCTreeForest();

	PCTree* makeTree(void);

	bool merge(PCTree* a, PCTree* b);

	void clear(void);

	operator const PCTreeRegistry<PCNode*>&() const { return nodeArrayRegistry; }
};
}
