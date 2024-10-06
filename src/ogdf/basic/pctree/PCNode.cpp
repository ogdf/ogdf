/** \file
 * \brief Implementation for ogdf::pc_tree::PCNode
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
#include <ogdf/basic/basic.h>
#include <ogdf/basic/pctree/PCEnum.h>
#include <ogdf/basic/pctree/PCNode.h>
#include <ogdf/basic/pctree/PCTree.h>
#include <ogdf/basic/pctree/PCTreeForest.h>
#include <ogdf/basic/pctree/PCTreeIterators.h>

#include <cstddef>
#include <ostream>
#include <utility>
#include <vector>

using namespace ogdf::pc_tree;

void PCNode::appendChild(PCNode* child, bool begin) {
	OGDF_ASSERT(child != nullptr);
	OGDF_ASSERT(m_forest == child->m_forest);
	OGDF_ASSERT(this != child);
	OGDF_ASSERT(isValidNode());
	OGDF_ASSERT(child->isValidNode(m_forest));
	child->setParent(this);
	m_childCount++;
	if (m_child1 == nullptr) {
		// new child of node without other children
		OGDF_ASSERT(m_child2 == nullptr);
		m_child1 = m_child2 = child;
	} else {
		// append new child
		PCNode*& outerChild = begin ? m_child1 : m_child2;
		outerChild->replaceSibling(nullptr, child);
		child->replaceSibling(nullptr, outerChild);
		outerChild = child;
	}
	OGDF_ASSERT(isValidNode());
	OGDF_ASSERT(child->isValidNode(m_forest));
}

void PCNode::insertBetween(PCNode* sib1, PCNode* sib2) {
	if (sib1 == nullptr && sib2 != nullptr) {
		insertBetween(sib2, sib1);
		return;
	}
	OGDF_ASSERT(isValidNode());
	OGDF_ASSERT(sib1 != nullptr);
	OGDF_ASSERT(sib1->m_forest == m_forest);
	OGDF_ASSERT(sib1->isValidNode(m_forest));
	PCNode* parent = sib1->getParent();
	if (sib2 == nullptr) {
		// append at one end of the list
	} else if (sib1->isSiblingOf(sib2)) {
		// insert within one list
		OGDF_ASSERT(parent != nullptr);
		OGDF_ASSERT(parent->isParentOf(sib2));
		if (sib1->isSiblingAdjacent(sib2)) {
			// normal case, both nodes are adjacent children of the same parent
			OGDF_ASSERT(sib2->isValidNode(m_forest));
			OGDF_ASSERT(parent->isValidNode(m_forest));
			setParent(parent);
			parent->m_childCount++;
			sib1->replaceSibling(sib2, this);
			sib2->replaceSibling(sib1, this);
			this->replaceSibling(nullptr, sib1);
			this->replaceSibling(nullptr, sib2);
			OGDF_ASSERT(isValidNode());
			OGDF_ASSERT(parent->isValidNode(m_forest));
			OGDF_ASSERT(sib1->isValidNode(m_forest));
			OGDF_ASSERT(sib2->isValidNode(m_forest));
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
	OGDF_ASSERT(parent->isValidNode(m_forest));
	OGDF_ASSERT(sib2 == nullptr || sib2->isValidNode(m_forest));
	setParent(parent);
	parent->m_childCount++;

	sib1->replaceSibling(nullptr, this);
	parent->replaceOuterChild(sib1, this);
	this->replaceSibling(nullptr, sib1);

	OGDF_ASSERT(isValidNode());
	OGDF_ASSERT(parent->isValidNode(m_forest));
	OGDF_ASSERT(sib1->isValidNode(m_forest));
	OGDF_ASSERT(sib2 == nullptr || sib2->isValidNode(m_forest));
}

void PCNode::detach() {
#ifdef OGDF_DEBUG
	OGDF_ASSERT(isValidNode());
	PCNode* parent = getParent();
	OGDF_ASSERT(parent == nullptr || parent->isValidNode(m_forest));
	forceDetach();
	OGDF_ASSERT(isValidNode());
	OGDF_ASSERT(parent == nullptr || parent->isValidNode(m_forest));
#else
	forceDetach();
#endif
}

void PCNode::forceDetach() {
	PCNode* parent = getParent();
	if (m_sibling1 != nullptr) {
		m_sibling1->replaceSibling(this, m_sibling2);
	} else if (parent != nullptr) {
		parent->replaceOuterChild(this, m_sibling2);
	}
	if (m_sibling2 != nullptr) {
		m_sibling2->replaceSibling(this, m_sibling1);
	} else if (parent != nullptr) {
		parent->replaceOuterChild(this, m_sibling1);
	}
	if (parent != nullptr) {
		OGDF_ASSERT(!parent->isChildOuter(this));
		parent->m_childCount--;
	}
	m_parentCNodeId = UNIONFINDINDEX_EMPTY;
	m_parentPNode = nullptr;
	m_sibling1 = m_sibling2 = nullptr;
}

void PCNode::replaceWith(PCNode* repl) {
	OGDF_ASSERT(repl != nullptr);
	OGDF_ASSERT(repl != this);
	OGDF_ASSERT(m_forest == repl->m_forest);
	OGDF_ASSERT(repl->isDetached());
	OGDF_ASSERT(repl->isValidNode(m_forest));
	OGDF_ASSERT(this != repl);
	PCNode* parent = getParent();
	OGDF_ASSERT(parent == nullptr || parent->isValidNode(m_forest));
	repl->m_parentCNodeId = m_parentCNodeId;
	repl->m_parentPNode = m_parentPNode;
	repl->m_sibling1 = m_sibling1;
	repl->m_sibling2 = m_sibling2;
	if (repl->m_sibling1 != nullptr) {
		repl->m_sibling1->replaceSibling(this, repl);
	}
	if (repl->m_sibling2 != nullptr) {
		repl->m_sibling2->replaceSibling(this, repl);
	}
	while (parent != nullptr && parent->isChildOuter(this)) {
		parent->replaceOuterChild(this, repl);
	}

	m_parentCNodeId = UNIONFINDINDEX_EMPTY;
	m_parentPNode = nullptr;
	m_sibling1 = m_sibling2 = nullptr;

	OGDF_ASSERT(isValidNode());
	OGDF_ASSERT(repl->isValidNode());
	OGDF_ASSERT(parent == nullptr || parent->isValidNode(m_forest));
}

void PCNode::mergeIntoParent() {
	OGDF_ASSERT(m_nodeType == PCNodeType::CNode);
	OGDF_ASSERT(isValidNode());
	PCNode* parent = getParent();
	OGDF_ASSERT(parent->m_nodeType == PCNodeType::CNode);
	OGDF_ASSERT(parent->isValidNode(m_forest));

	UnionFindIndex pcid = m_forest->m_parents.link(m_nodeListIndex, parent->m_nodeListIndex);
	if (pcid == this->m_nodeListIndex) {
		std::swap(m_forest->m_cNodes[this->m_nodeListIndex],
				m_forest->m_cNodes[parent->m_nodeListIndex]);
		std::swap(this->m_nodeListIndex, parent->m_nodeListIndex);
	} else {
		OGDF_ASSERT(pcid == parent->m_nodeListIndex);
	}
	parent->m_childCount += m_childCount - 1;

	if (m_sibling1 != nullptr) {
		m_sibling1->replaceSibling(this, m_child1);
		m_child1->replaceSibling(nullptr, m_sibling1);
	} else {
		parent->replaceOuterChild(this, m_child1);
	}
	if (m_sibling2 != nullptr) {
		m_sibling2->replaceSibling(this, m_child2);
		m_child2->replaceSibling(nullptr, m_sibling2);
	} else {
		parent->replaceOuterChild(this, m_child2);
	}

	m_child1 = m_child2 = nullptr;
	m_sibling1 = m_sibling2 = nullptr;
	m_parentCNodeId = UNIONFINDINDEX_EMPTY;
	m_parentPNode = nullptr;
	m_childCount = 0;

	OGDF_ASSERT(parent->isValidNode(m_forest));
}

void PCNode::replaceSibling(PCNode* oldS, PCNode* newS) {
	OGDF_ASSERT((newS == nullptr) || (m_forest == newS->m_forest));
	OGDF_ASSERT(newS != this);
	if (oldS == m_sibling1) {
		m_sibling1 = newS;
	} else {
		OGDF_ASSERT(oldS == m_sibling2);
		m_sibling2 = newS;
	}
}

void PCNode::replaceOuterChild(PCNode* oldC, PCNode* newC) {
	OGDF_ASSERT((newC == nullptr) || (m_forest == newC->m_forest));
	OGDF_ASSERT(newC != this);
	if (oldC == m_child1) {
		m_child1 = newC;
	} else {
		OGDF_ASSERT(oldC == m_child2);
		m_child2 = newC;
	}
}

void ogdf::pc_tree::proceedToNextSibling(PCNode*& pred, PCNode*& curr) {
	pred = curr->getNextSibling(pred);
	std::swap(pred, curr);
}

PCNode* PCNode::getNextSibling(const PCNode* pred) const {
	OGDF_ASSERT(pred == nullptr || isSiblingOf(pred));
	if (pred == m_sibling1) {
		return m_sibling2;
	} else {
		OGDF_ASSERT(pred == m_sibling2);
		return m_sibling1;
	}
}

PCNode* PCNode::getOtherOuterChild(const PCNode* child) const {
	OGDF_ASSERT(isParentOf(child));
	if (child == m_child1) {
		return m_child2;
	} else {
		OGDF_ASSERT(child == m_child2);
		return m_child1;
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
			l_next = m_child1;
		} else {
			OGDF_ASSERT(isParentOf(curr));
			l_next = (curr->m_sibling1 != nullptr) ? curr->m_sibling1 : curr->m_sibling2;
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
	if (ofForest && m_forest != ofForest) {
		return false;
	}

#ifdef OGDF_DEBUG
	if (m_parentCNodeId == UNIONFINDINDEX_EMPTY && m_parentPNode == nullptr) {
		OGDF_ASSERT(m_sibling1 == nullptr && m_sibling2 == nullptr);
	} else {
		OGDF_ASSERT(m_parentCNodeId == UNIONFINDINDEX_EMPTY || m_parentPNode == nullptr);
		PCNode* parent = getParent();
		OGDF_ASSERT(parent->m_forest == m_forest);
		int null_sibs = 0;
		if (m_sibling1 == nullptr) {
			null_sibs++;
		} else {
			OGDF_ASSERT(m_sibling1->isSiblingAdjacent(this));
		}
		if (m_sibling2 == nullptr) {
			null_sibs++;
		} else {
			OGDF_ASSERT(m_sibling2->isSiblingAdjacent(this));
		}
		int outsides = 0;
		if (parent->m_child1 == this) {
			outsides++;
		}
		if (parent->m_child2 == this) {
			outsides++;
		}
		OGDF_ASSERT(null_sibs == outsides);
		OGDF_ASSERT((null_sibs > 0) == isOuterChild());
		if (isOuterChild()) {
			OGDF_ASSERT(parent->isChildOuter(this));
		}
	}

	if (m_childCount == 0) {
		OGDF_ASSERT(m_child1 == nullptr);
		OGDF_ASSERT(m_child2 == nullptr);
	} else if (m_childCount == 1) {
		OGDF_ASSERT(m_child1 != nullptr);
		OGDF_ASSERT(m_child2 != nullptr);
		OGDF_ASSERT(m_child1 == m_child2);
	} else {
		OGDF_ASSERT(m_childCount >= 2);
		OGDF_ASSERT(m_child1 != nullptr);
		OGDF_ASSERT(m_child2 != nullptr);
		OGDF_ASSERT(m_child1 != m_child2);
	}
#endif

	if (m_nodeType == PCNodeType::CNode) {
		OGDF_ASSERT(m_forest->m_cNodes.at(m_nodeListIndex) == this);
		return (size_t)m_forest->m_parents.find(m_nodeListIndex) == m_nodeListIndex;
	} else if (m_nodeType == PCNodeType::Leaf) {
		OGDF_ASSERT(getDegree() <= 1);
		// OGDF_ASSERT(forest->leaves.at(nodeListIndex) == this);
		return true;
	} else {
		OGDF_ASSERT(m_nodeType == PCNodeType::PNode);
		return true;
	}
}

PCNode* PCNode::getParent() const {
	if (m_parentPNode != nullptr) {
		OGDF_ASSERT(m_parentCNodeId == UNIONFINDINDEX_EMPTY);
		OGDF_ASSERT(m_parentPNode->m_nodeType == PCNodeType::PNode
				|| m_parentPNode->m_nodeType == PCNodeType::Leaf);
		return m_parentPNode;
	} else if (m_parentCNodeId != UNIONFINDINDEX_EMPTY) {
		m_parentCNodeId = m_forest->m_parents.find(m_parentCNodeId);
		OGDF_ASSERT(m_forest->m_cNodes.at(m_parentCNodeId) != nullptr);
		PCNode* parent = m_forest->m_cNodes[m_parentCNodeId];
		OGDF_ASSERT(parent != this);
		OGDF_ASSERT(parent->m_nodeType == PCNodeType::CNode);
		OGDF_ASSERT(parent->m_nodeListIndex == m_parentCNodeId);
		return parent;
	} else {
		return nullptr;
	}
}

void PCNode::setParent(PCNode* parent) {
	OGDF_ASSERT(isDetached());
	OGDF_ASSERT(parent != nullptr);
	if (parent->m_nodeType == PCNodeType::CNode) {
		m_parentCNodeId = parent->m_nodeListIndex;
	} else {
		m_parentPNode = parent;
	}
}

// only allowed on root children
void PCNode::rotateChildOutside(bool child1) {
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
		if (child1) {
			parent->m_child1 = this;
			parent->m_child2 = (getSibling1() != nullptr) ? getSibling1() : getSibling2();

			OGDF_ASSERT(parent->m_child2 != nullptr);
		} else {
			parent->m_child2 = this;
			parent->m_child1 = (getSibling1() != nullptr) ? getSibling1() : getSibling2();

			OGDF_ASSERT(parent->m_child1 != nullptr);
		}

		// cut connection between node and one of its siblings and promote them to outer children
		parent->m_child1->replaceSibling(parent->m_child2, nullptr);
		parent->m_child2->replaceSibling(parent->m_child1, nullptr);
	} else if (child1 && parent->m_child1 != this) {
		PCNode* otherChild = parent->m_child1;
		parent->m_child1 = this;
		parent->m_child2 = otherChild;
	} else if (!child1 && parent->m_child2 != this) {
		PCNode* otherChild = parent->m_child2;
		parent->m_child2 = this;
		parent->m_child1 = otherChild;
	}

	OGDF_ASSERT((child1 && parent->m_child1 == this) || (!child1 && parent->m_child2 == this));
	OGDF_ASSERT(parent->m_child1->isOuterChild());
	OGDF_ASSERT(parent->m_child2->isOuterChild());
	OGDF_ASSERT(isOuterChild());
}

void PCNode::checkTimestamp() const {
	OGDF_ASSERT(isValidNode());
	if (m_forest->m_timestamp != m_timestamp) {
		OGDF_ASSERT(m_forest->m_timestamp > m_timestamp);
		m_label = NodeLabel::Unknown;
		m_timestamp = m_forest->m_timestamp;
		if (!isLeaf()) {
			m_temp.clear();
		}
	}
}

std::ostream& ogdf::pc_tree::operator<<(std::ostream& os, const ogdf::pc_tree::PCNode& node) {
	return os << &node;
}

std::ostream& ogdf::pc_tree::operator<<(std::ostream& os, const ogdf::pc_tree::PCNode* node) {
	if (node == nullptr) {
		return os << "null-Node";
	}
	os << node->m_nodeType << " " << node->index();
	os << " with " << node->m_childCount << " children";
	os << " [";
	int c = 0;
	for (PCNode *pred = nullptr, *curr = node->m_child1; curr != nullptr;
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
		os << parent->m_nodeType << " " << parent->index();
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
	m_node->proceedToNextNeighbor(m_pred, m_curr);
	return *this;
}

PCNodeIterator PCNodeIterator::operator++(int) {
	PCNodeIterator before = *this;
	m_node->proceedToNextNeighbor(m_pred, m_curr);
	return before;
}

bool PCNodeIterator::isParent() {
	return m_node != nullptr && m_curr != nullptr && m_curr->isParentOf(m_node);
}

PCNodeIterator PCNodeChildrenIterable::begin() const noexcept {
	return PCNodeIterator(m_node, nullptr, m_node->m_child1);
}

PCNodeIterator PCNodeChildrenIterable::end() const noexcept {
	if (m_node->m_child1 != nullptr && !m_node->isDetached()) {
		return PCNodeIterator(m_node, m_node->m_child2, m_node->getParent());
	} else {
		return PCNodeIterator(m_node, m_node->m_child2, m_node->m_child1);
	}
}

unsigned long PCNodeChildrenIterable::count() const { return m_node->m_childCount; }

PCNodeIterator PCNodeNeighborsIterable::begin() const noexcept {
	return PCNodeIterator(m_node, nullptr, m_first);
}

PCNodeIterator PCNodeNeighborsIterable::end() const noexcept {
	PCNode* second = m_node->getNextNeighbor(nullptr, m_first);
	PCNode* last = m_node->getNextNeighbor(second, m_first);
	return PCNodeIterator(m_node, last, m_first);
}

unsigned long PCNodeNeighborsIterable::count() const { return m_node->getDegree(); }
