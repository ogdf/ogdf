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

#include <ogdf/basic/pctree/PCNode.h>
#include <ogdf/basic/pctree/PCTree.h>
#include <ogdf/basic/pctree/PCTreeForest.h>
#include <ogdf/basic/pctree/PCTreeIterators.h>

using namespace pc_tree;

void PCNode::appendChild(PCNode* p_node, bool begin) {
	OGDF_ASSERT(p_node != nullptr);
	OGDF_ASSERT(forest == p_node->forest);
	OGDF_ASSERT(this != p_node);
	OGDF_ASSERT(isValidNode());
	OGDF_ASSERT(p_node->isValidNode(forest));
	p_node->setParent(this);
	childCount++;
	if (child1 == nullptr) {
		// new child of node without other children
		OGDF_ASSERT(child2 == nullptr);
		child1 = child2 = p_node;
	} else {
		// append new child
		PCNode*& outerChild = begin ? child1 : child2;
		outerChild->replaceSibling(nullptr, p_node);
		p_node->replaceSibling(nullptr, outerChild);
		outerChild = p_node;
	}
	OGDF_ASSERT(isValidNode());
	OGDF_ASSERT(p_node->isValidNode(forest));
}

void PCNode::insertBetween(PCNode* sib1, PCNode* sib2) {
	if (sib1 == nullptr && sib2 != nullptr) {
		insertBetween(sib2, sib1);
		return;
	}
	OGDF_ASSERT(isValidNode());
	OGDF_ASSERT(sib1 != nullptr);
	OGDF_ASSERT(sib1->forest == forest);
	OGDF_ASSERT(sib1->isValidNode(forest));
	PCNode* parent = sib1->getParent();
	if (sib2 == nullptr) {
		// append at one end of the list
	} else if (sib1->isSiblingOf(sib2)) {
		// insert within one list
		OGDF_ASSERT(parent != nullptr);
		OGDF_ASSERT(parent->isParentOf(sib2));
		if (sib1->isSiblingAdjacent(sib2)) {
			// normal case, both nodes are adjacent children of the same parent
			OGDF_ASSERT(sib2->isValidNode(forest));
			OGDF_ASSERT(parent->isValidNode(forest));
			setParent(parent);
			parent->childCount++;
			sib1->replaceSibling(sib2, this);
			sib2->replaceSibling(sib1, this);
			this->replaceSibling(nullptr, sib1);
			this->replaceSibling(nullptr, sib2);
			OGDF_ASSERT(isValidNode());
			OGDF_ASSERT(parent->isValidNode(forest));
			OGDF_ASSERT(sib1->isValidNode(forest));
			OGDF_ASSERT(sib2->isValidNode(forest));
			return;
		} else {
			// wrap around if no parent is present and both nodes are outer nodes of the same parent
			OGDF_ASSERT(parent->isDetached());
			OGDF_ASSERT(parent->isChildOuter(sib1));
			OGDF_ASSERT(parent->isChildOuter(sib2));
		}
	} else if (parent != nullptr && sib2->isParentOf(parent)) {
		// sib2 is the parent of `parent`, sib1 is an outer child of `parent`
	} else {
		std::swap(sib1, sib2);
		parent = sib1->getParent();
		// now we should have the same situation as before
		OGDF_ASSERT(parent != nullptr);
		OGDF_ASSERT(sib2->isParentOf(parent));
	}

	OGDF_ASSERT(parent != nullptr);
	OGDF_ASSERT(parent->isValidNode(forest));
	OGDF_ASSERT(sib2 == nullptr || sib2->isValidNode(forest));
	setParent(parent);
	parent->childCount++;

	sib1->replaceSibling(nullptr, this);
	parent->replaceOuterChild(sib1, this);
	this->replaceSibling(nullptr, sib1);

	OGDF_ASSERT(isValidNode());
	OGDF_ASSERT(parent->isValidNode(forest));
	OGDF_ASSERT(sib1->isValidNode(forest));
	OGDF_ASSERT(sib2 == nullptr || sib2->isValidNode(forest));
}

void PCNode::detach() {
#ifdef OGDF_DEBUG
	OGDF_ASSERT(isValidNode());
	PCNode* parent = getParent();
	OGDF_ASSERT(parent == nullptr || parent->isValidNode(forest));
	forceDetach();
	OGDF_ASSERT(isValidNode());
	OGDF_ASSERT(parent == nullptr || parent->isValidNode(forest));
#else
	forceDetach();
#endif
}

void PCNode::forceDetach() {
	PCNode* parent = getParent();
	if (sibling1 != nullptr) {
		sibling1->replaceSibling(this, sibling2);
	} else if (parent != nullptr) {
		parent->replaceOuterChild(this, sibling2);
	}
	if (sibling2 != nullptr) {
		sibling2->replaceSibling(this, sibling1);
	} else if (parent != nullptr) {
		parent->replaceOuterChild(this, sibling1);
	}
	if (parent != nullptr) {
		OGDF_ASSERT(!parent->isChildOuter(this));
		parent->childCount--;
	}
	parentCNodeId = UNIONFINDINDEX_EMPTY;
	parentPNode = nullptr;
	sibling1 = sibling2 = nullptr;
}

void PCNode::replaceWith(PCNode* p_node) {
	OGDF_ASSERT(p_node != nullptr);
	OGDF_ASSERT(p_node != this);
	OGDF_ASSERT(forest == p_node->forest);
	OGDF_ASSERT(p_node->isDetached());
	OGDF_ASSERT(p_node->isValidNode(forest));
	OGDF_ASSERT(this != p_node);
	PCNode* parent = getParent();
	OGDF_ASSERT(parent == nullptr || parent->isValidNode(forest));
	p_node->parentCNodeId = parentCNodeId;
	p_node->parentPNode = parentPNode;
	p_node->sibling1 = sibling1;
	p_node->sibling2 = sibling2;
	if (p_node->sibling1 != nullptr) {
		p_node->sibling1->replaceSibling(this, p_node);
	}
	if (p_node->sibling2 != nullptr) {
		p_node->sibling2->replaceSibling(this, p_node);
	}
	while (parent != nullptr && parent->isChildOuter(this)) {
		parent->replaceOuterChild(this, p_node);
	}

	parentCNodeId = UNIONFINDINDEX_EMPTY;
	parentPNode = nullptr;
	sibling1 = sibling2 = nullptr;

	OGDF_ASSERT(isValidNode());
	OGDF_ASSERT(p_node->isValidNode());
	OGDF_ASSERT(parent == nullptr || parent->isValidNode(forest));
}

void PCNode::mergeIntoParent() {
	OGDF_ASSERT(nodeType == PCNodeType::CNode);
	OGDF_ASSERT(isValidNode());
	PCNode* p_parent = getParent();
	OGDF_ASSERT(p_parent->nodeType == PCNodeType::CNode);
	OGDF_ASSERT(p_parent->isValidNode(forest));

	UnionFindIndex pcid = forest->parents.link(nodeListIndex, p_parent->nodeListIndex);
	if (pcid == this->nodeListIndex) {
		std::swap(forest->cNodes[this->nodeListIndex], forest->cNodes[p_parent->nodeListIndex]);
		std::swap(this->nodeListIndex, p_parent->nodeListIndex);
	} else {
		OGDF_ASSERT(pcid == p_parent->nodeListIndex);
	}
	p_parent->childCount += childCount - 1;

	if (sibling1 != nullptr) {
		sibling1->replaceSibling(this, child1);
		child1->replaceSibling(nullptr, sibling1);
	} else {
		p_parent->replaceOuterChild(this, child1);
	}
	if (sibling2 != nullptr) {
		sibling2->replaceSibling(this, child2);
		child2->replaceSibling(nullptr, sibling2);
	} else {
		p_parent->replaceOuterChild(this, child2);
	}

	child1 = child2 = nullptr;
	sibling1 = sibling2 = nullptr;
	parentCNodeId = UNIONFINDINDEX_EMPTY;
	parentPNode = nullptr;
	childCount = 0;

	OGDF_ASSERT(p_parent->isValidNode(forest));
}

void PCNode::replaceSibling(PCNode* oldS, PCNode* newS) {
	OGDF_ASSERT((newS == nullptr) || (forest == newS->forest));
	OGDF_ASSERT(newS != this);
	if (oldS == sibling1) {
		sibling1 = newS;
	} else {
		OGDF_ASSERT(oldS == sibling2);
		sibling2 = newS;
	}
}

void PCNode::replaceOuterChild(PCNode* oldC, PCNode* newC) {
	OGDF_ASSERT((newC == nullptr) || (forest == newC->forest));
	OGDF_ASSERT(newC != this);
	if (oldC == child1) {
		child1 = newC;
	} else {
		OGDF_ASSERT(oldC == child2);
		child2 = newC;
	}
}

void pc_tree::proceedToNextSibling(PCNode*& pred, PCNode*& curr) {
	pred = curr->getNextSibling(pred);
	std::swap(pred, curr);
}

PCNode* PCNode::getNextSibling(const PCNode* pred) const {
	OGDF_ASSERT(pred == nullptr || isSiblingOf(pred));
	if (pred == sibling1) {
		return sibling2;
	} else {
		OGDF_ASSERT(pred == sibling2);
		return sibling1;
	}
}

PCNode* PCNode::getOtherOuterChild(const PCNode* child) const {
	OGDF_ASSERT(isParentOf(child));
	if (child == child1) {
		return child2;
	} else {
		OGDF_ASSERT(child == child2);
		return child1;
	}
}

void PCNode::proceedToNextNeighbor(PCNode*& pred, PCNode*& curr) const {
	pred = getNextNeighbor(pred, curr);
	std::swap(pred, curr);
}

PCNode* PCNode::getNextNeighbor(const PCNode* pred, const PCNode* curr) const {
	OGDF_ASSERT(curr != nullptr);

	PCNode* l_next;
	PCNode* l_parent = getParent();

	if (pred == nullptr) {
		if (curr == l_parent) {
			l_next = child1;
		} else {
			OGDF_ASSERT(isParentOf(curr));
			l_next = (curr->sibling1 != nullptr) ? curr->sibling1 : curr->sibling2;
		}
	} else {
		// find next sibling
		if (pred->isSiblingOf(curr)) {
			OGDF_ASSERT(isParentOf(curr));
			if (pred->isSiblingAdjacent(curr)) {
				l_next = curr->getNextSibling(pred);
			} else {
				OGDF_ASSERT(isDetached());
				OGDF_ASSERT(isChildOuter(pred));
				OGDF_ASSERT(isChildOuter(curr));
				l_next = curr->getNextSibling(nullptr);
			}

			// find second or second to last child
		} else if (pred == l_parent) {
			OGDF_ASSERT(isParentOf(curr));
			l_next = curr->getNextSibling(nullptr);

			// find first or last child
		} else {
			OGDF_ASSERT(curr == l_parent);
			OGDF_ASSERT(isParentOf(pred));
			l_next = getOtherOuterChild(pred);
		}
	}

	if (l_next == nullptr) {
		if (l_parent == nullptr) {
			return getOtherOuterChild(curr);
		} else {
			return l_parent;
		}
	} else {
		OGDF_ASSERT(isParentOf(l_next));
		return l_next;
	}
}

bool PCNode::areNeighborsAdjacent(const PCNode* neigh1, const PCNode* neigh2) const {
	OGDF_ASSERT(neigh1 != nullptr);
	OGDF_ASSERT(neigh2 != nullptr);
	OGDF_ASSERT(neigh1 != neigh2);

	if (isParentOf(neigh1) && isParentOf(neigh2)) {
		return neigh1->isSiblingAdjacent(neigh2)
				|| (isDetached() && isChildOuter(neigh1) && isChildOuter(neigh2));

	} else {
		PCNode* parent = getParent();
		if (neigh1 == parent) {
			OGDF_ASSERT(isParentOf(neigh2));
			return isChildOuter(neigh2);

		} else {
			OGDF_ASSERT(neigh2 == parent);
			OGDF_ASSERT(isParentOf(neigh1));
			return isChildOuter(neigh1);
		}
	}
}

bool PCNode::isValidNode(const PCTreeForest* ofForest) const {
	if (ofForest && forest != ofForest) {
		return false;
	}
	if (parentCNodeId == UNIONFINDINDEX_EMPTY && parentPNode == nullptr) {
		OGDF_ASSERT(sibling1 == nullptr && sibling2 == nullptr);
	} else {
		OGDF_ASSERT(parentCNodeId == UNIONFINDINDEX_EMPTY || parentPNode == nullptr);
		PCNode* parent = getParent();
		OGDF_ASSERT(parent->forest == forest);
		int null_sibs = 0;
		if (sibling1 == nullptr) {
			null_sibs++;
		} else {
			OGDF_ASSERT(sibling1->isSiblingAdjacent(this));
		}
		if (sibling2 == nullptr) {
			null_sibs++;
		} else {
			OGDF_ASSERT(sibling2->isSiblingAdjacent(this));
		}
		int outsides = 0;
		if (parent->child1 == this) {
			outsides++;
		}
		if (parent->child2 == this) {
			outsides++;
		}
		OGDF_ASSERT(null_sibs == outsides);
		OGDF_ASSERT((null_sibs > 0) == isOuterChild());
		if (isOuterChild()) {
			OGDF_ASSERT(parent->isChildOuter(this));
		}
	}

	if (childCount == 0) {
		OGDF_ASSERT(child1 == nullptr);
		OGDF_ASSERT(child2 == nullptr);
	} else if (childCount == 1) {
		OGDF_ASSERT(child1 != nullptr);
		OGDF_ASSERT(child2 != nullptr);
		OGDF_ASSERT(child1 == child2);
	} else {
		OGDF_ASSERT(childCount >= 2);
		OGDF_ASSERT(child1 != nullptr);
		OGDF_ASSERT(child2 != nullptr);
		OGDF_ASSERT(child1 != child2);
	}

	if (nodeType == PCNodeType::CNode) {
		OGDF_ASSERT(forest->cNodes.at(nodeListIndex) == this);
		return forest->parents.find(nodeListIndex) == nodeListIndex;
	} else if (nodeType == PCNodeType::Leaf) {
		OGDF_ASSERT(getDegree() <= 1);
		// OGDF_ASSERT(forest->leaves.at(nodeListIndex) == this);
		return true;
	} else {
		OGDF_ASSERT(nodeType == PCNodeType::PNode);
		return true;
	}
}

PCNode* PCNode::getParent() const {
	if (parentPNode != nullptr) {
		OGDF_ASSERT(parentCNodeId == UNIONFINDINDEX_EMPTY);
		OGDF_ASSERT(parentPNode->nodeType == PCNodeType::PNode
				|| parentPNode->nodeType == PCNodeType::Leaf);
		return parentPNode;
	} else if (parentCNodeId != UNIONFINDINDEX_EMPTY) {
		parentCNodeId = forest->parents.find(parentCNodeId);
		OGDF_ASSERT(forest->cNodes.at(parentCNodeId) != nullptr);
		PCNode* parent = forest->cNodes[parentCNodeId];
		OGDF_ASSERT(parent != this);
		OGDF_ASSERT(parent->nodeType == PCNodeType::CNode);
		OGDF_ASSERT(parent->nodeListIndex == parentCNodeId);
		return parent;
	} else {
		return nullptr;
	}
}

void PCNode::setParent(PCNode* parent) {
	OGDF_ASSERT(isDetached());
	OGDF_ASSERT(parent != nullptr);
	if (parent->nodeType == PCNodeType::CNode) {
		parentCNodeId = parent->nodeListIndex;
	} else {
		parentPNode = parent;
	}
}

// only allowed on root children
void PCNode::rotateChildOutside(bool p_child1) {
	PCNode* parent = getParent();
	OGDF_ASSERT(parent != nullptr);

	// parent must be root node
	OGDF_ASSERT(parent->getParent() == nullptr);

	/* rotate node to child1/2 of parent, children might be flipped but that's perfectly fine for a
     * CNode (and PNode as well) if parent is the root node */
	if (!isOuterChild()) {
		PCNode* left = parent->getChild1();
		PCNode* right = parent->getChild2();

		OGDF_ASSERT(right->isOuterChild());
		OGDF_ASSERT(left->isOuterChild());

		// create connection between outer children
		left->replaceSibling(nullptr, right);
		right->replaceSibling(nullptr, left);

		// find any existing sibling of node to later promote to child1/2
		if (p_child1) {
			parent->child1 = this;
			parent->child2 = (getSibling1() != nullptr) ? getSibling1() : getSibling2();

			OGDF_ASSERT(parent->child2 != nullptr);
		} else {
			parent->child2 = this;
			parent->child1 = (getSibling1() != nullptr) ? getSibling1() : getSibling2();

			OGDF_ASSERT(parent->child1 != nullptr);
		}

		// cut connection between node and one of its siblings and promote them to outer children
		parent->child1->replaceSibling(parent->child2, nullptr);
		parent->child2->replaceSibling(parent->child1, nullptr);
	} else if (p_child1 && parent->child1 != this) {
		PCNode* otherChild = parent->child1;
		parent->child1 = this;
		parent->child2 = otherChild;
	} else if (!p_child1 && parent->child2 != this) {
		PCNode* otherChild = parent->child2;
		parent->child2 = this;
		parent->child1 = otherChild;
	}

	OGDF_ASSERT((p_child1 && parent->child1 == this) || (!p_child1 && parent->child2 == this));
	OGDF_ASSERT(parent->child1->isOuterChild());
	OGDF_ASSERT(parent->child2->isOuterChild());
	OGDF_ASSERT(isOuterChild());
}

void PCNode::checkTimestamp() const {
	OGDF_ASSERT(isValidNode());
	if (forest->timestamp != timestamp) {
		OGDF_ASSERT(forest->timestamp > timestamp);
		label = NodeLabel::Unknown;
		timestamp = forest->timestamp;
		if (!isLeaf()) {
			temp.clear();
		}
	}
}

std::ostream& operator<<(std::ostream& os, const PCNode& node) { return os << &node; }

std::ostream& operator<<(std::ostream& os, const PCNode* node) {
	if (node == nullptr) {
		return os << "null-Node";
	}
	//    if (node->timestamp == node->tree->timestamp)
	//        os << node->temp.label << " ";
	os << node->nodeType << " " << node->index();
	os << " with " << node->childCount << " children";
	//    if (node->childCount == 1)
	//        os << " [" << node->child1->index() << "]";
	//    else if (node->childCount == 2)
	//        os << " [" << node->child1->index() << ", " << node->child2->index() << "]";
	//    else if (node->childCount > 2)
	//        os << " [" << node->child1->index() << ", ..., " << node->child2->index() << "]";
	os << " [";
	int c = 0;
	for (PCNode *pred = nullptr, *curr = node->child1; curr != nullptr;
			proceedToNextSibling(pred, curr)) {
		if (c > 0) {
			os << ", ";
		}
		os << curr->index();
		c++;
	}
	os << "]";
	PCNode* parent = node->getParent();
	if (parent != nullptr) {
		os << " and parent ";
		//        if (parent->timestamp == parent->tree->timestamp)
		//            os << parent->temp.label << " ";
		os << parent->nodeType << " " << parent->index();
	} else {
		os << " and no parent";
	}
	return os;
}

PCNodeChildrenIterable PCNode::children() { return PCNodeChildrenIterable(this); }

PCNodeNeighborsIterable PCNode::neighbors(PCNode* startWith) {
	return PCNodeNeighborsIterable(this, startWith);
}

PCNodeIterator& PCNodeIterator::operator++() {
	node->proceedToNextNeighbor(pred, curr);
	return *this;
}

PCNodeIterator PCNodeIterator::operator++(int) {
	PCNodeIterator before = *this;
	node->proceedToNextNeighbor(pred, curr);
	return before;
}

bool PCNodeIterator::isParent() {
	return node != nullptr && curr != nullptr && curr->isParentOf(node);
}

PCNodeIterator PCNodeChildrenIterable::begin() const noexcept {
	return PCNodeIterator(node, nullptr, node->child1);
}

PCNodeIterator PCNodeChildrenIterable::end() const noexcept {
	if (node->child1 != nullptr && !node->isDetached()) {
		return PCNodeIterator(node, node->child2, node->getParent());
	} else {
		return PCNodeIterator(node, node->child2, node->child1);
	}
}

unsigned long PCNodeChildrenIterable::count() const { return node->childCount; }

PCNodeIterator PCNodeNeighborsIterable::begin() const noexcept {
	return PCNodeIterator(node, nullptr, first);
}

PCNodeIterator PCNodeNeighborsIterable::end() const noexcept {
	PCNode* second = node->getNextNeighbor(nullptr, first);
	PCNode* last = node->getNextNeighbor(second, first);
	return PCNodeIterator(node, last, first);
}

unsigned long PCNodeNeighborsIterable::count() const { return node->getDegree(); }
