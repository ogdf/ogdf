/** \file
 * \brief A node in a PC-tree that is either a P-node, C-node or leaf.
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

#include <ogdf/basic/basic.h>
#include <ogdf/basic/memory.h>
#include <ogdf/basic/pctree/PCEnum.h>
#include <ogdf/basic/pctree/PCTreeForest.h>
#include <ogdf/basic/pctree/util/IntrusiveList.h>

#include <array>
#include <cstddef>
#include <iosfwd>
#include <new>
#include <utility>
#include <vector>

namespace ogdf::pc_tree {
class PCTree;
struct PCNodeChildrenIterable;
struct PCNodeNeighborsIterable;

/**
 * A node in a PC-tree that is either a P-node, C-node or leaf.
 * See https://doi.org/10.15475/cpatp.2024 Figure 8.3 for a visualization of the doubly-linked tree structure stored in each node and more details on the changes made and temporary information stored by PCTree::makeConsecutive().
 * Important terminology:
 * - child: direct descendant of this node in the tree
 * - outer child: first or last child
 * - sibling: other node with the same direct parent
 * - adjacent sibling: predecessor or successor in parent's list of children
 * - neighbors: all children and the parent
 */
class OGDF_EXPORT PCNode : public IntrusiveList<PCNode>::node {
	friend OGDF_EXPORT std::ostream&(operator<<)(std::ostream&, const ogdf::pc_tree::PCTree*);
	friend OGDF_EXPORT std::ostream&(operator<<)(std::ostream&, const ogdf::pc_tree::PCNode*);

	friend class PCTree;
	friend class PCTreeForest;
	friend struct PCNodeChildrenIterable;
	friend struct PCNodeNeighborsIterable;

public:
	/**
	 * Temporary information used during each step of the PCTree::makeConsecutive() update operation.
	 */
	struct TempInfo {
		PCNode *predPartial = nullptr, *nextPartial = nullptr;
		PCNode* tpPred = nullptr;
		PCNode* tpPartialPred = nullptr;
		size_t tpPartialHeight = 0;
		PCNode* tpSucc = nullptr;
		std::vector<PCNode*> fullNeighbors;
		PCNode *ebEnd1 = nullptr, *fbEnd1 = nullptr, *fbEnd2 = nullptr, *ebEnd2 = nullptr;

		void replaceNeighbor(PCNode* oldNeigh, PCNode* newNeigh) {
			if (tpPred == oldNeigh) {
				tpPred = newNeigh;
			}
			if (tpPartialPred == oldNeigh) {
				tpPartialPred = newNeigh;
			}
			if (tpSucc == oldNeigh) {
				tpSucc = newNeigh;
			}
			if (ebEnd1 == oldNeigh) {
				ebEnd1 = newNeigh;
			}
			if (ebEnd2 == oldNeigh) {
				ebEnd2 = newNeigh;
			}
			if (fbEnd1 == oldNeigh) {
				fbEnd1 = newNeigh;
			}
			if (fbEnd2 == oldNeigh) {
				fbEnd2 = newNeigh;
			}
		}

		void clear() {
			nextPartial = predPartial = nullptr;
			tpPred = tpPartialPred = tpSucc = nullptr;
			ebEnd1 = fbEnd1 = fbEnd2 = ebEnd2 = nullptr;
			tpPartialHeight = 0;
			fullNeighbors.clear();
		}
	};

	using LeafUserData = std::array<void*, sizeof(TempInfo) / sizeof(void*)>;

private:
	// index in registry
	size_t m_id;

	// global
	PCTreeForest* m_forest;

	// private
	UnionFindIndex m_nodeListIndex = UNIONFINDINDEX_EMPTY;
	PCNodeType m_nodeType;
	PCNode* m_parentPNode = nullptr;
	mutable UnionFindIndex m_parentCNodeId = UNIONFINDINDEX_EMPTY;
	PCNode* m_sibling1 = nullptr;
	PCNode* m_sibling2 = nullptr;
	PCNode* m_child1 = nullptr;
	PCNode* m_child2 = nullptr;
	size_t m_childCount = 0;
	mutable NodeLabel m_label = NodeLabel::Unknown;
	mutable size_t m_timestamp = 0;

	// leaves need no temp info, so they can easily store user data
	union {
		mutable TempInfo m_temp;
		LeafUserData m_userData;
	};

	PCNode(PCTreeForest* forest, size_t id, PCNodeType nodeType)
		: IntrusiveList<PCNode>::node(), m_id(id), m_forest(forest), m_nodeType(nodeType) {
		if (nodeType == PCNodeType::Leaf) {
			new (&m_userData) LeafUserData;
		} else {
			new (&m_temp) TempInfo;
		}
	}

	~PCNode() {
		if (m_nodeType == PCNodeType::Leaf) {
			m_userData.~array();
		} else {
			m_temp.~TempInfo();
		}
	}

public:
	/**
	 * @name Tree structure methods
	 * These methods allow modifying the tree structure or embedding, e.g., when manually constructing a PC-tree.
	 */
	//! @{

	/**
	 * Append a (detached) child node to the begin or end of this nodes' children.
	 */
	void appendChild(PCNode* child, bool begin = false);

	/**
	 * Insert a (detached) child node directly between two adjacent children of this node.
	 */
	void insertBetween(PCNode* sib1, PCNode* sib2);

	/**
	 * Detach this node from its parent. Invalidates but does not change the PC-forest of the node, so don't forget to re-attach it somewhere else in a PC-tree of the same forest to make it valid again.
	 */
	void detach();

	/**
	 * Swaps this node inplace with a (detached) other one. Afterwards, this node will be detached.
	 */
	void replaceWith(PCNode* repl);

	/**
	 * Merges this C-node into its C-node parent.
	 */
	void mergeIntoParent();

	/**
	 * Reverse the stored order of children.
	 */
	void flip() { std::swap(m_child1, m_child2); }

private:
	/**
	 * Notify this node that one of its adjacent siblings changed.
	 */
	void replaceSibling(PCNode* oldS, PCNode* newS);

	/**
	 * Make this node an outer child of its parent. Only works for children of the root node.
	 */
	void rotateChildOutside(bool child1 = true);

	/**
	 * Notify this node that one of its outer children was replaced.
	 */
	void replaceOuterChild(PCNode* oldC, PCNode* newC);

	/**
	 * Notify this node that it has a new parent.
	 */
	void setParent(PCNode* parent);

	/**
	 * detach() without performing checks.
	 */
	void forceDetach();

	/**
	 * Overwrite the type of this node without updating any other data structures.
	 */
	void changeType(PCNodeType newType) {
		if (m_nodeType == PCNodeType::Leaf && newType != PCNodeType::Leaf) {
			m_userData.~array();
			new (&m_temp) TempInfo;
		} else if (m_nodeType != PCNodeType::Leaf && newType == PCNodeType::Leaf) {
			m_temp.~TempInfo();
			new (&m_userData) LeafUserData;
		}
		m_nodeType = newType;
	}

	//! @}

public:
	/**
	 * @name Iterator methods
	 * These methods allow easily walking around the PC-tree.
	 */
	//! @{

	/**
	 * Given the left or right sibling \p pred, return the adjacent sibling on the other side.
	 */
	PCNode* getNextSibling(const PCNode* pred) const;

	/**
	 * Given one outer child, return the outer child on the other side.
	 */
	PCNode* getOtherOuterChild(const PCNode* child) const;

	/**
	 * Method to walk the cyclic order of all neighbors, i.e., all children plus the parent, where this node's parent is considered to be adjacent to this node's two outer children.
	 * The returned node is the one of the two nodes adjacent to \p curr in this cyclic order that is not \p pred, or an arbitrary one of the two if \p pred is null.
	 */
	PCNode* getNextNeighbor(const PCNode* pred, const PCNode* curr) const;

	/**
	 * Iteration-convenience version of getNextNeighbor() that updates the variables \p pred to \p curr and \p curr to the value returned by getNextNeighbor(pred, curr).
	 */
	void proceedToNextNeighbor(PCNode*& pred, PCNode*& curr) const;

	/**
	 * @return the parent node of this node, or null if this node is the root or currently detached. If the parent is a C-node, this requires a look-up in the union-find data structure.
	 */
	PCNode* getParent() const;

	/**
	 * @return iterable for all children
	 */
	PCNodeChildrenIterable children();

	/**
	 * @return iterable for all children plus the parent, optionally selecting a starting node in this cyclic order
	 */
	PCNodeNeighborsIterable neighbors(PCNode* first = nullptr);

	//! @}

public:
	/**
	 * @name Structure check methods
	 * These methods help with checking the current status of this node w.r.t. its surrounding tree structure.
	 */
	//! @{

	/**
	 * @return \c true if this node has no parent, i.e., it is the root of its PC-tree or needs to be attached to some node first before the tree can become valid again.
	 */
	bool isDetached() const {
		if (m_parentCNodeId == UNIONFINDINDEX_EMPTY && m_parentPNode == nullptr) {
			return true;
		} else {
			OGDF_ASSERT(m_parentCNodeId == UNIONFINDINDEX_EMPTY || m_parentPNode == nullptr);
			return false;
		}
	}

	/**
	 * Perform multiple debug checks and return \c true if getForest == \p ofForest
	 */
	bool isValidNode(const PCTreeForest* ofForest = nullptr) const;

	bool isLeaf() const { return m_nodeType == PCNodeType::Leaf; }

	/**
	 * @return \c true, if other->getParent() == this
	 */
	bool isParentOf(const PCNode* other) const {
		OGDF_ASSERT(other != nullptr);
		OGDF_ASSERT(m_forest == other->m_forest);
		return other->getParent() == this;
	}

	/**
	 * @return \c true, if this->getParent() == other->getParent()
	 */
	bool isSiblingOf(const PCNode* other) const {
		OGDF_ASSERT(other != nullptr);
		OGDF_ASSERT(m_forest == other->m_forest);
		return this->getParent() == other->getParent();
	}

	/**
	 * @return \c true, if this->getSibling1() == sibling or this->getSibling2() == sibling
	 */
	bool isSiblingAdjacent(const PCNode* sibling) const {
		OGDF_ASSERT(isSiblingOf(sibling));
		OGDF_ASSERT(this != sibling);
		return m_sibling1 == sibling || m_sibling2 == sibling;
	}

	/**
	 * @return \c true if neigh1 and neigh2 are children of this node and isSiblingAdjacent(neigh1, neigh2) is true,
	 *   or if one of the passed nodes is the parent of this node and the other one is an outer child of this node
	 */
	bool areNeighborsAdjacent(const PCNode* neigh1, const PCNode* neigh2) const;

	/**
	 * @return \c true, if \p child is an outer child of this node, i.e., if this->getChild1() == child or this->getChild2() == child
	 */
	bool isChildOuter(const PCNode* child) const {
		OGDF_ASSERT(isParentOf(child));
		return m_child1 == child || m_child2 == child;
	}

	/**
	 * @return \c true, if this node is an outer child of its parent, i.e., if this->getSibling1() == nullptr or this->getSibling2() == child
	 */
	bool isOuterChild() const { return m_sibling1 == nullptr || m_sibling2 == nullptr; }

	//! @}

public:
	/**
	 * @name makeConsecutive-related temporary information
	 * These methods provide access to temporary information used during a PCTree::makeConsecutive() call.
	 */
	//! @{

	const TempInfo& constTempInfo() const {
		checkTimestamp();
		OGDF_ASSERT(!isLeaf());
		return m_temp;
	}

	bool isFull() const { return getLabel() == NodeLabel::Full; }

	NodeLabel getLabel() const {
		// this operation does not reset the temp info
		return m_forest->m_timestamp == m_timestamp ? m_label : NodeLabel::Empty;
	}

	void setLabel(NodeLabel l) {
		checkTimestamp();
		m_label = l;
	}

private:
	// these methods are slightly faster if we already called checkTimestamp()
	inline NodeLabel getLabelUnchecked() const {
		OGDF_ASSERT(m_forest->m_timestamp == m_timestamp);
		return m_label;
	}

	inline void setLabelUnchecked(NodeLabel l) {
		OGDF_ASSERT(m_forest->m_timestamp == m_timestamp);
		m_label = l;
	}

public:
	/**
	 * @return the user data that can be stored in leaves
	 */
	LeafUserData& leafUserData() {
		OGDF_ASSERT(isLeaf());
		return m_userData;
	}

	/**
	 * @return the user data that can be stored in leaves
	 */
	const LeafUserData& leafUserData() const {
		OGDF_ASSERT(isLeaf());
		return m_userData;
	}

private:
	void checkTimestamp() const;

	TempInfo& tempInfo() {
		checkTimestamp();
		OGDF_ASSERT(!isLeaf());
		return m_temp;
	}

	size_t addFullNeighbor(PCNode* fullNeigh) {
		checkTimestamp();
		OGDF_ASSERT(!isLeaf());
		OGDF_ASSERT(fullNeigh->isFull());
		m_temp.fullNeighbors.push_back(fullNeigh);
		return m_temp.fullNeighbors.size();
	}

	PCNode*& getFullNeighInsertionPoint(PCNode* nonFullNeigh) {
		checkTimestamp();
		OGDF_ASSERT(!isLeaf());
		OGDF_ASSERT(nonFullNeigh != nullptr);
		if (nonFullNeigh == m_temp.ebEnd1) {
			OGDF_ASSERT(areNeighborsAdjacent(m_temp.ebEnd1, m_temp.fbEnd1));
			return m_temp.fbEnd1;
		} else {
			OGDF_ASSERT(nonFullNeigh == m_temp.ebEnd2);
			OGDF_ASSERT(areNeighborsAdjacent(m_temp.ebEnd2, m_temp.fbEnd2));
			return m_temp.fbEnd2;
		}
	}

	//! @}

public:
	/**
	 * @name Getters
	 */
	//! @{
	size_t index() const { return m_id; }

	PCNodeType getNodeType() const { return m_nodeType; }

	size_t getChildCount() const { return m_childCount; }

	size_t getDegree() const { return isDetached() ? m_childCount : m_childCount + 1; }

	PCNode* getChild1() const { return m_child1; }

	PCNode* getChild2() const { return m_child2; }

	/**
	 * Check whether this node has only one child and return it.
	 */
	PCNode* getOnlyChild() const {
		OGDF_ASSERT(m_childCount == 1);
		return m_child1;
	}

	PCNode* getSibling1() const { return m_sibling1; }

	PCNode* getSibling2() const { return m_sibling2; }

	PCTreeForest* getForest() const { return m_forest; }

	//! @}

	OGDF_NEW_DELETE
};

/**
 * Iteration-convenience version of PCNode::getNextSibling() that updates the variables \p pred to \p curr and \p curr to the value returned by PCNode::getNextSibling(pred, curr).
 */
OGDF_EXPORT void proceedToNextSibling(PCNode*& pred, PCNode*& curr);
}
