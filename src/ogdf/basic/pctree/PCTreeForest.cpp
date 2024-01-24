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

#include <ogdf/basic/pctree/PCTree.h>
#include <ogdf/basic/pctree/PCTreeForest.h>

using namespace pc_tree;

PCTreeForest::~PCTreeForest() {
	clear();
#ifdef OGDF_PCTREE_REUSE_NODES
	while (reusableNodes) {
		PCNode* tmp = reusableNodes;
		reusableNodes = reusableNodes->parentPNode;
		delete tmp;
	}
#endif
}

// forest auto deletes allocated trees when destructed
PCTree* PCTreeForest::makeTree() {
	PCTree* tree = new PCTree(this);
	trees.push_back(tree);

	return tree;
}

// merges tree b into a at root node, tree b is invalidated
bool PCTreeForest::merge(PCTree* a, PCTree* b) {
	OGDF_ASSERT(a != nullptr && b != nullptr);
	OGDF_ASSERT(a != b);
	OGDF_ASSERT(a->forest == b->forest);
	OGDF_ASSERT(a->externalForest && b->externalForest);
	OGDF_ASSERT(a->getRootNode()->getNodeType() == PCNodeType::Leaf
			&& b->getRootNode()->getNodeType() == PCNodeType::Leaf);

	a->observers.splice(a->observers.end(), b->observers);
	a->leaves.splice(a->leaves.end(), b->leaves);

	a->pNodeCount += b->pNodeCount;
	a->cNodeCount += b->cNodeCount;

	OGDF_ASSERT(a->getRootNode()->getChild1()->getNodeType() == PCNodeType::PNode);

	PCNode* nodeDest = a->getRootNode()->getChild1();
	PCNode* nodeSrc = b->getRootNode()->getChild1();
	nodeSrc->detach();
	a->destroyNode(b->getRootNode());

	bool contractedPNode = false;

	if (nodeSrc->getChildCount() > 1) {
		if (nodeDest->getChildCount() > 0) {
			nodeDest->appendChild(nodeSrc);
		} else {
			nodeDest->detach();
			a->destroyNode(nodeDest);
			a->getRootNode()->appendChild(nodeSrc);
			contractedPNode = true;
		}
	} else if (nodeSrc->getChildCount() == 1) {
		PCNode* onlyChild = nodeSrc->getChild1();
		OGDF_ASSERT(onlyChild->isLeaf());
		onlyChild->detach();
		a->destroyNode(nodeSrc);
		nodeDest->appendChild(onlyChild);
	} else {
		OGDF_ASSERT(nodeSrc->getChildCount() == 0);
		a->destroyNode(nodeSrc);
	}

	b->rootNode = nullptr;
	b->pNodeCount = b->cNodeCount = 0;

	OGDF_ASSERT(a->checkValid());
	return contractedPNode;
}

void PCTreeForest::clear() {
	if (autodelete) {
		for (auto* k : trees) {
			delete k;
		}
	}

	trees.clear();
	trees.shrink_to_fit();
	cNodes.clear();
	cNodes.shrink_to_fit();
	parents.init();
	nextNodeId = 0;
	timestamp = 0;
}
