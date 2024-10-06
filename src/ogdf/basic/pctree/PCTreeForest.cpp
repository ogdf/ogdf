/** \file
 * \brief Implementation for ogdf::pc_tree::PCTreeForest
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

#include <ogdf/basic/DisjointSets.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/pctree/PCNode.h>
#include <ogdf/basic/pctree/PCRegistry.h>
#include <ogdf/basic/pctree/PCTree.h>
#include <ogdf/basic/pctree/PCTreeForest.h>

#include <vector>

using namespace ogdf::pc_tree;

PCTreeForest::~PCTreeForest() {
	clear();
#ifdef OGDF_PCTREE_REUSE_NODES
	while (m_reusableNodes) {
		PCNode* tmp = m_reusableNodes;
		m_reusableNodes = m_reusableNodes->m_parentPNode;
		delete tmp;
	}
#endif
}

// forest auto deletes allocated trees when destructed
PCTree* PCTreeForest::makeTree() {
	PCTree* tree = new PCTree(this);
	m_trees.push_back(tree);

	return tree;
}

void PCTreeForest::clear() {
	if (m_autodelete) {
		for (auto* k : m_trees) {
			delete k;
		}
	}

	m_trees.clear();
	m_trees.shrink_to_fit();
	m_cNodes.clear();
	m_cNodes.shrink_to_fit();
	m_parents.init();
	m_nextNodeId = 0;
	m_timestamp = 0;
}

bool PCTreeRegistry::isKeyAssociated(PCNode* key) const {
#ifdef OGDF_DEBUG
	return key && key->getForest() == m_pForest;
#else
	return key;
#endif
}

int PCTreeRegistry::calculateArraySize(int add) const {
	return ogdf::calculateTableSize(m_pForest->m_nextNodeId + add);
}

int PCTreeRegistry::maxKeyIndex() const { return m_pForest->m_nextNodeId - 1; }
