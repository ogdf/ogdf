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

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/pctree/PCEnum.h>
#include <ogdf/basic/pctree/PCNode.h>
#include <ogdf/basic/pctree/PCRegistry.h>
#include <ogdf/basic/pctree/PCTreeForest.h>
#include <ogdf/basic/pctree/PCTreeIterators.h>
#include <ogdf/basic/pctree/util/IntrusiveList.h>

#include <deque>
#include <list>
#include <sstream>
#include <vector>

namespace pc_tree {
bool isTrivialRestriction(int restSize, int leafCount);

int factorial(int n);

namespace uid_utils {
void nodeToID(std::ostream& os, PCNode* n, int pos);

void nodeToPosition(std::ostream& os, PCNode* n, int pos);

void leafToID(std::ostream& os, PCNode* n, int pos);

void leafToPosition(std::ostream& os, PCNode* n, int pos);

bool compareNodesByID(PCNode* a, PCNode* b);
}

class PCTree {
	friend std::ostream&(::operator<<)(std::ostream&, const pc_tree::PCTree*);

	friend std::ostream&(::operator<<)(std::ostream&, const pc_tree::PCNode*);

	friend class PCNode;

	template<class Key>
	friend class PCTreeRegistry;

	friend class PCTreeForest;

public:
	struct Observer;

private:
	// global
	PCTreeForest* forest = nullptr;

	// needs merge
	size_t pNodeCount = 0;
	size_t cNodeCount = 0;
	IntrusiveList<PCNode> leaves;

	// needs special merge
	PCNode* rootNode = nullptr;

	// temp variables
	int partialCount = 0;
	int terminalPathLength = 0;
	PCNode* firstPartial = nullptr;
	PCNode* lastPartial = nullptr;
	PCNode* apexCandidate = nullptr;
	bool apexCandidateIsFix = false;
	PCNode* apexTPPred2 = nullptr;

	// private
	bool externalForest = true;
	std::list<Observer*> observers;

public:
	explicit PCTree() : forest(new PCTreeForest()), externalForest(false) {};

	explicit PCTree(PCTreeForest* p_forest) : forest(p_forest), externalForest(true) {};

	explicit PCTree(int leafNum, std::vector<PCNode*>* added = nullptr);

	explicit PCTree(const PCTree& other, PCTreeNodeArray<PCNode*>& nodeMapping,
			bool keep_ids = false);

	explicit PCTree(const std::string& str, bool keep_ids = false);

	virtual ~PCTree();

public: // Node creation / destruction
	PCNode* newNode(PCNodeType type, PCNode* parent = nullptr, int id = -1);

	void destroyNode(PCNode*& node) {
		destroyNode((PCNode* const&)node);
		node = nullptr;
	}

	void destroyNode(PCNode* const& node);

	void insertLeaves(int count, PCNode* parent, std::vector<PCNode*>* added = nullptr);

	void replaceLeaf(int leafCount, PCNode* leaf, std::vector<PCNode*>* added = nullptr);

	PCNode* mergeLeaves(std::vector<PCNode*>& consecutiveLeaves, bool assumeConsecutive = false) {
		return mergeLeaves(consecutiveLeaves.begin(), consecutiveLeaves.end(), assumeConsecutive);
	}

	template<typename It>
	PCNode* mergeLeaves(It begin, It end, bool assumeConsecutive = false) {
		OGDF_ASSERT(begin != end);

		if (!assumeConsecutive && !makeConsecutive(begin, end)) {
			return nullptr;
		}

		// Remove all consecutive leaves except the first one.
		It back = prev(end);
		for (auto it = begin; it != back; ++it) {
			destroyLeaf(*it);
		}

		OGDF_HEAVY_ASSERT(
				checkValid()); // TODO this breaks if the replacement reduces the leaves to two or less?

		// Return the remaining leaf.
		return *back;
	}

	void destroyLeaf(PCNode* leaf);

	PCNode* setRoot(PCNode* newRoot);

	PCNode* changeRoot(PCNode* newRoot);

	PCNodeType changeNodeType(PCNode* node, PCNodeType newType);

	void insertTree(PCNode* at, PCTree* tree);

private:
	void unregisterNode(PCNode* node);

	void registerNode(PCNode* node);

public: // Restrictions
	bool isTrivialRestriction(int size) const;

	bool makeConsecutive(std::initializer_list<PCNode*> consecutiveLeaves) {
		return makeConsecutive(consecutiveLeaves.begin(), consecutiveLeaves.end());
	}

	bool makeConsecutive(std::vector<PCNode*>& consecutiveLeaves) {
		return makeConsecutive(consecutiveLeaves.begin(), consecutiveLeaves.end());
	}

	template<typename It>
	bool makeConsecutive(It begin, It end) {
		FullLeafIter iter = [&begin, &end]() { return NextFullLeaf<It>(begin, end); };
		for (auto obs : observers) {
			obs->makeConsecutiveCalled(*this, iter);
		}

		OGDF_HEAVY_ASSERT(checkValid());
		resetTempData();

#ifdef OGDF_DEBUG
		for (auto it = begin; it != end; ++it) {
			PCNode* leaf = *it;
			OGDF_ASSERT(leaf);
			OGDF_ASSERT(leaf->isLeaf());
			OGDF_ASSERT(leaf->forest == forest);
		}
#endif
		if (isTrivialRestriction(end - begin)) {
			for (auto obs : observers) {
				obs->makeConsecutiveDone(*this, Observer::Stage::Trivial, true);
			}
			return true;
		}

		// PC_PROFILE_ENTER(1, "label");
		markFull(begin, end);
		// PC_PROFILE_EXIT(1, "label");

		return makeFullNodesConsecutive();
	}

	void resetTempData() {
		forest->timestamp++;
		firstPartial = lastPartial = nullptr;
		partialCount = 0;
		apexCandidate = nullptr;
		apexCandidateIsFix = false;
		terminalPathLength = 0;
		apexTPPred2 = nullptr;
	}

	/**
	 * Only marks leaves full, does not update partial/full info of parents.
	 * Use markFull to also update parents.
	 */
	template<typename It>
	void markLeavesFull(It begin, It end) {
		for (auto it = begin; it != end; ++it) {
			PCNode* leaf = *it;
			OGDF_ASSERT(leaf);
			OGDF_ASSERT(leaf->isLeaf());
			OGDF_ASSERT(leaf->forest == forest);
			leaf->setLabel(NodeLabel::Full);
		}
	}

	/**
	 * Attention: We no longer use a queue to defer processing of partial/full parents to after
	 * all leaves are done, but now directly make parents full if all of their children are full.
	 */
	template<typename It>
	void markFull(It begin, It end, std::vector<PCNode*>* fullNodeOrder = nullptr) {
		if (fullNodeOrder != nullptr) {
			fullNodeOrder->reserve(cNodeCount + pNodeCount);
		}

		for (auto it = begin; it != end; ++it) {
			PCNode* full_parent = markFull(*it, fullNodeOrder);
			while (full_parent != nullptr) {
				full_parent = markFull(full_parent, fullNodeOrder);
			}
		}
	}

	bool makeFullNodesConsecutive();

private:
	PCNode* markFull(PCNode* full_node, std::vector<PCNode*>* fullNodeOrder = nullptr);

	bool findTerminalPath();

	void updateSingletonTerminalPath();

	PCNode* createCentralNode();

	int updateTerminalPath(PCNode* central, PCNode* tpNeigh);

	// labeling / TP finding

	void addPartialNode(PCNode* partial);

	void removePartialNode(PCNode* partial);

	bool checkTPPartialCNode(PCNode* node);

	size_t findEndOfFullBlock(PCNode* node, PCNode* pred, PCNode* curr, PCNode*& fullEnd,
			PCNode*& emptyEnd) const;

	bool setApexCandidate(PCNode* ac, bool fix = false);

	// update

	void replaceTPNeigh(PCNode* central, PCNode* oldTPNeigh, PCNode* newTPNeigh,
			PCNode* newFullNeigh, PCNode* otherEndOfFullBlock);

	PCNode* splitOffFullPNode(PCNode* node, bool skip_parent);

public: // Intersect
	bool intersect(PCTree& other, PCTreeNodeArray<PCNode*>& mapping);

private:
	bool findNodeRestrictions(PCTree& applyTo, PCTreeNodeArray<PCNode*>& mapping,
			PCTreeNodeArray<std::vector<PCNode*>>& blockNodes,
			PCTreeNodeArray<std::vector<PCNode*>>& subtreeNodes,
			PCTreeNodeArray<PCNode*>& leafPartner, PCTreeNodeArray<bool>& isFront);

	void restoreSubtrees(PCTreeNodeArray<std::vector<PCNode*>>& blockNodes,
			PCTreeNodeArray<std::vector<PCNode*>>& subtreeNodes,
			PCTreeNodeArray<PCNode*>& leafPartner, PCTreeNodeArray<bool>& isFront);

public: // Getters
	operator const PCTreeRegistry<PCNode*>&() const { return forest->nodeArrayRegistry; }

	[[nodiscard]] PCTreeForest* getForest() const { return forest; }

	[[nodiscard]] bool isTrivial() const;

	[[nodiscard]] const IntrusiveList<PCNode>& getLeaves() const { return leaves; }

	[[nodiscard]] size_t getLeafCount() const { return leaves.size(); };

	[[nodiscard]] size_t getPNodeCount() const { return pNodeCount; }

	[[nodiscard]] size_t getCNodeCount() const { return cNodeCount; }

	[[nodiscard]] PCNode* getRootNode() const { return rootNode; }

	int getTerminalPathLength() const { return terminalPathLength; }

	FilteringPCTreeDFS allNodes() const { return FilteringPCTreeDFS(*this, rootNode); }

	FilteringPCTreeDFS innerNodes() const {
		return FilteringPCTreeDFS(*this, rootNode, [](PCNode* node) { return !node->isLeaf(); });
	}

	template<typename Container>
	void currentLeafOrder(Container& container) const {
		for (PCNode* leaf : FilteringPCTreeDFS(*this, rootNode)) {
			if (leaf->isLeaf()) {
				container.push_back(leaf);
			}
		}
		OGDF_ASSERT(container.size() == leaves.size());
	}

	std::vector<PCNode*> currentLeafOrder() const {
		std::vector<PCNode*> container;
		currentLeafOrder(container);
		return container;
	}

	bool checkValid(bool allow_small_deg = true) const;

	bool isValidOrder(const std::vector<PCNode*>& order) const;

	void getTree(ogdf::Graph& tree, ogdf::GraphAttributes* g_a,
			PCTreeNodeArray<ogdf::node>& pc_repr, ogdf::NodeArray<PCNode*>* g_repr = nullptr,
			bool mark_full = false, bool show_sibs = false) const;

	void getRestrictions(std::vector<std::vector<PCNode*>>& restrictions,
			PCNode* fixedLeaf = nullptr) const;

	template<typename R>
	R possibleOrders() const;

	std::ostream& uniqueID(std::ostream&,
			const std::function<void(std::ostream& os, PCNode*, int)>& printNode = uid_utils::nodeToID,
			const std::function<bool(PCNode*, PCNode*)>& compareNodes = uid_utils::compareNodesByID);

	std::string uniqueID(
			const std::function<void(std::ostream& os, PCNode*, int)>& printNode = uid_utils::nodeToID,
			[[maybe_unused]] const std::function<bool(PCNode*, PCNode*)>& compareNodes =
					uid_utils::compareNodesByID) {
		std::stringstream sb;
		uniqueID(sb, printNode, compareNodes);
		return sb.str();
	}

public: // Observers
	template<typename It>
	struct NextFullLeaf {
		It it;
		It end;

		NextFullLeaf(It p_it, It an_end) : it(p_it), end(an_end) { }

		PCNode* operator()() {
			if (it == end) {
				return nullptr;
			}
			PCNode* n = *it;
			++it;
			return n;
		}
	};

	using FullLeafIter = std::function<std::function<PCNode*()>()>;

	struct Observer {
		enum class Stage { Trivial, NoPartials, InvalidTP, SingletonTP, Done };

		virtual void onNodeCreate([[maybe_unused]] PCNode* node) {};

		virtual void makeConsecutiveCalled([[maybe_unused]] PCTree& tree,
				[[maybe_unused]] FullLeafIter consecutiveLeaves) {};

		virtual void labelsAssigned([[maybe_unused]] PCTree& tree,
				[[maybe_unused]] PCNode* p_firstPartial, [[maybe_unused]] PCNode* p_lastPartial,
				[[maybe_unused]] int p_partialCount) {};

		virtual void terminalPathFound([[maybe_unused]] PCTree& tree, [[maybe_unused]] PCNode* apex,
				[[maybe_unused]] PCNode* p_apexTPPred2, [[maybe_unused]] int p_terminalPathLength) {};

		virtual void centralCreated([[maybe_unused]] PCTree& tree,
				[[maybe_unused]] PCNode* central) {};

		virtual void beforeMerge([[maybe_unused]] PCTree& tree, [[maybe_unused]] int count,
				[[maybe_unused]] PCNode* tpNeigh) {};

		virtual void afterMerge([[maybe_unused]] PCTree& tree, [[maybe_unused]] PCNode* successor,
				[[maybe_unused]] PCNode* mergedNode) {};

		virtual void whenPNodeMerged([[maybe_unused]] PCTree& tree, [[maybe_unused]] PCNode* tpNeigh,
				[[maybe_unused]] PCNode* tpPred, [[maybe_unused]] PCNode* fullNeigh) {};

		virtual void whenCNodeMerged([[maybe_unused]] PCTree& tree,
				[[maybe_unused]] PCNode* tpNeigh, [[maybe_unused]] bool tpNeighSiblingsFlipped,
				[[maybe_unused]] PCNode* fullNeigh, [[maybe_unused]] PCNode* fullOuterChild) {};

		virtual void fullNodeSplit([[maybe_unused]] PCTree& tree,
				[[maybe_unused]] PCNode* fullNode) {};

		virtual void makeConsecutiveDone([[maybe_unused]] PCTree& tree,
				[[maybe_unused]] Stage stage, [[maybe_unused]] bool success) {};

		virtual void onApexMoved([[maybe_unused]] PCTree& tree,
				[[maybe_unused]] PCNode* p_apexCandidate, [[maybe_unused]] PCNode* central,
				[[maybe_unused]] PCNode* parent) {};

		virtual void nodeDeleted([[maybe_unused]] PCTree& tree,
				[[maybe_unused]] PCNode* toBeDeleted) {};
		virtual void nodeReplaced([[maybe_unused]] PCTree& tree, [[maybe_unused]] PCNode* replaced,
				[[maybe_unused]] PCNode* replacement) {};
	};

	struct LoggingObserver : public Observer {
		void makeConsecutiveCalled(PCTree& tree, FullLeafIter consecutiveLeaves) override;

		void labelsAssigned(PCTree& tree, PCNode* firstPartial, PCNode* lastPartial,
				int partialCount) override;

		void terminalPathFound(PCTree& tree, PCNode* apex, PCNode* apexTPPred2,
				int terminalPathLength) override;

		void centralCreated(PCTree& tree, PCNode* central) override;

		void beforeMerge(PCTree& tree, int count, PCNode* tpNeigh) override;

		void makeConsecutiveDone(PCTree& tree, Stage stage, bool success) override;
	};

	std::list<Observer*>::const_iterator addObserver(Observer* observer) {
		observers.push_back(observer);
		return --observers.end();
	}

	void removeObserver(std::list<Observer*>::const_iterator it) { observers.erase(it); }

	void removeObserver(Observer* observer) { observers.remove(observer); }
};

template<class Key>
int PCTreeRegistry<Key>::calculateArraySize(int add) const {
	return ogdf::calculateTableSize(m_pForest->nextNodeId + add);
}

template<class Key>
int PCTreeRegistry<Key>::maxKeyIndex() const {
	return m_pForest->nextNodeId - 1;
}
}
