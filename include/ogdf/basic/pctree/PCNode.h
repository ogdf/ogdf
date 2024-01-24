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

#include <ogdf/basic/basic.h>
#include <ogdf/basic/memory.h>
#include <ogdf/basic/pctree/PCEnum.h>
#include <ogdf/basic/pctree/PCTreeForest.h>
#include <ogdf/basic/pctree/util/IntrusiveList.h>

#include <array>
#include <list>
#include <vector>

namespace pc_tree {
struct OGDF_EXPORT PCNodeChildrenIterable;
struct OGDF_EXPORT PCNodeNeighborsIterable;

class OGDF_EXPORT PCNode : public IntrusiveList<PCNode>::node {
public:
	friend class PCTree;

	friend class PCTreeForest;

	friend struct PCNodeChildrenIterable;
	friend struct PCNodeNeighborsIterable;

	friend std::ostream&(::operator<<)(std::ostream&, const pc_tree::PCTree*);

	friend std::ostream&(::operator<<)(std::ostream&, const pc_tree::PCNode*);

	struct TempInfo {
		PCNode *predPartial = nullptr, *nextPartial = nullptr;
		PCNode* tpPred = nullptr;
		PCNode* tpPartialPred = nullptr;
		int tpPartialHeight = 0;
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
	// used for serializing and debug output?
	int id;

	// global
	PCTreeForest* forest;

	// private
	UnionFindIndex nodeListIndex = UNIONFINDINDEX_EMPTY;
	PCNodeType nodeType;
	PCNode* parentPNode = nullptr;
	mutable UnionFindIndex parentCNodeId = UNIONFINDINDEX_EMPTY;
	PCNode* sibling1 = nullptr;
	PCNode* sibling2 = nullptr;
	PCNode* child1 = nullptr;
	PCNode* child2 = nullptr;
	size_t childCount = 0;
	mutable NodeLabel label = NodeLabel::Unknown;
	mutable int timestamp = 0;

	union {
		mutable TempInfo temp;
		LeafUserData userData;
	};

	PCNode(PCTreeForest* p_forest, int p_id, PCNodeType p_nodeType)
		: IntrusiveList<PCNode>::node(), id(p_id), forest(p_forest), nodeType(p_nodeType) {
		if (p_nodeType == PCNodeType::Leaf) {
			new (&userData) LeafUserData;
		} else {
			new (&temp) TempInfo;
		}
	}

	~PCNode() {
		if (nodeType == PCNodeType::Leaf) {
			userData.~array();
		} else {
			temp.~TempInfo();
		}
	}

public:
	void appendChild(PCNode* p_node, bool begin = false);

	void insertBetween(PCNode* sib1, PCNode* sib2);

	void detach();

	void replaceWith(PCNode* p_node);

	void mergeIntoParent();

	void flip() { std::swap(child1, child2); }

	void replaceSibling(PCNode* oldS, PCNode* newS);

	void rotateChildOutside(bool p_child1 = true);

private:
	void replaceOuterChild(PCNode* oldC, PCNode* newC);

	void setParent(PCNode* parent);

	void forceDetach();

	void changeType(PCNodeType newType) {
		if (nodeType == PCNodeType::Leaf && newType != PCNodeType::Leaf) {
			userData.~array();
			new (&temp) TempInfo;
		} else if (nodeType != PCNodeType::Leaf && newType == PCNodeType::Leaf) {
			temp.~TempInfo();
			new (&userData) LeafUserData;
		}
		nodeType = newType;
	}

public:
	PCNode* getNextSibling(const PCNode* pred) const;

	PCNode* getOtherOuterChild(const PCNode* child) const;

	PCNode* getNextNeighbor(const PCNode* pred, const PCNode* curr) const;

	void proceedToNextNeighbor(PCNode*& pred, PCNode*& curr) const;

	PCNode* getParent() const;

	PCNodeChildrenIterable children();

	PCNodeNeighborsIterable neighbors(PCNode* first = nullptr);

public:
	bool isDetached() const {
		if (parentCNodeId == UNIONFINDINDEX_EMPTY && parentPNode == nullptr) {
			return true;
		} else {
			OGDF_ASSERT(parentCNodeId == UNIONFINDINDEX_EMPTY || parentPNode == nullptr);
			return false;
		}
	}

	bool isValidNode(const PCTreeForest* ofForest = nullptr) const;

	bool isLeaf() const { return nodeType == PCNodeType::Leaf; }

	bool isParentOf(const PCNode* other) const {
		OGDF_ASSERT(other != nullptr);
		OGDF_ASSERT(forest == other->forest);
		return other->getParent() == this;
	}

	bool isSiblingOf(const PCNode* other) const {
		OGDF_ASSERT(other != nullptr);
		OGDF_ASSERT(forest == other->forest);
		return this->getParent() == other->getParent();
	}

	bool isSiblingAdjacent(const PCNode* sibling) const {
		OGDF_ASSERT(isSiblingOf(sibling));
		OGDF_ASSERT(this != sibling);
		return sibling1 == sibling || sibling2 == sibling;
	}

	bool areNeighborsAdjacent(const PCNode* neigh1, const PCNode* neigh2) const;

	bool isChildOuter(const PCNode* child) const {
		OGDF_ASSERT(isParentOf(child));
		return child1 == child || child2 == child;
	}

	bool isOuterChild() const { return sibling1 == nullptr || sibling2 == nullptr; }

public:
	const TempInfo& constTempInfo() const {
		checkTimestamp();
		OGDF_ASSERT(!isLeaf());
		return temp;
	}

	bool isFull() const { return getLabel() == NodeLabel::Full; }

	NodeLabel getLabel() const {
		// this operation does not reset the temp info
		return forest->timestamp == timestamp ? label : NodeLabel::Empty;
	}

	void setLabel(NodeLabel l) {
		checkTimestamp();
		label = l;
	}

	NodeLabel getLabelUnchecked() const {
		OGDF_ASSERT(forest->timestamp == timestamp);
		return label;
	}

	void setLabelUnchecked(NodeLabel l) {
		OGDF_ASSERT(forest->timestamp == timestamp);
		label = l;
	}

	LeafUserData& leafUserData() {
		OGDF_ASSERT(isLeaf());
		return userData;
	}

	const LeafUserData& leafUserData() const {
		OGDF_ASSERT(isLeaf());
		return userData;
	}

	PCNode* getFullNeighInsertionPointConst(PCNode* nonFullNeigh) {
		return getFullNeighInsertionPoint(nonFullNeigh);
	}

private:
	void checkTimestamp() const;

	TempInfo& tempInfo() {
		checkTimestamp();
		OGDF_ASSERT(!isLeaf());
		return temp;
	}

	size_t addFullNeighbor(PCNode* fullNeigh) {
		checkTimestamp();
		OGDF_ASSERT(!isLeaf());
		OGDF_ASSERT(fullNeigh->isFull());
		temp.fullNeighbors.push_back(fullNeigh);
		return temp.fullNeighbors.size();
	}

	PCNode*& getFullNeighInsertionPoint(PCNode* nonFullNeigh) {
		checkTimestamp();
		OGDF_ASSERT(!isLeaf());
		OGDF_ASSERT(nonFullNeigh != nullptr);
		if (nonFullNeigh == temp.ebEnd1) {
			OGDF_ASSERT(areNeighborsAdjacent(temp.ebEnd1, temp.fbEnd1));
			return temp.fbEnd1;
		} else {
			OGDF_ASSERT(nonFullNeigh == temp.ebEnd2);
			OGDF_ASSERT(areNeighborsAdjacent(temp.ebEnd2, temp.fbEnd2));
			return temp.fbEnd2;
		}
	}

public:
	int index() const { return id; }

	PCNodeType getNodeType() const { return nodeType; }

	int getChildCount() const { return childCount; }

	size_t getDegree() const { return isDetached() ? childCount : childCount + 1; }

	PCNode* getChild1() const { return child1; }

	PCNode* getChild2() const { return child2; }

	PCNode* getOnlyChild() const {
		OGDF_ASSERT(childCount == 1);
		return child1;
	}

	PCNode* getSibling1() const { return sibling1; }

	PCNode* getSibling2() const { return sibling2; }

	PCTreeForest* getForest() const { return forest; }

	OGDF_NEW_DELETE
};

OGDF_EXPORT void proceedToNextSibling(PCNode*& pred, PCNode*& curr);
}
