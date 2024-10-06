/** \file
 * \brief The main class of the PC-tree.
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
#include <ogdf/basic/basic.h>
#include <ogdf/basic/internal/config_autogen.h>
#include <ogdf/basic/pctree/PCEnum.h>
#include <ogdf/basic/pctree/PCNode.h>
#include <ogdf/basic/pctree/PCRegistry.h>
#include <ogdf/basic/pctree/PCTreeForest.h>
#include <ogdf/basic/pctree/PCTreeIterators.h>
#include <ogdf/basic/pctree/util/IntrusiveList.h>

#include <cstddef>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <list>
#include <sstream>
#include <string>
#include <vector>

namespace ogdf {
class GraphAttributes;
} // namespace ogdf

namespace ogdf::pc_tree {
/**
 * @return \c true if calling PCTree::makeConsecutive() with \p restSize out of \p leafCount total leaves never requires changes to the tree.
 *   This is the case for \p restSize values 0, 1, \p leafCount - 1, and \p leafCount.
 */
OGDF_EXPORT bool isTrivialRestriction(int restSize, int leafCount);

OGDF_EXPORT int factorial(int n);

#ifdef OGDF_DEBUG
/**
 * Allows controlling the frequency of full-tree consistency checks in heavy debug mode.
 * When set to a non-zero value n, only every n'th check will be performed.
 * When set to 1, every check will be executed.
 * When set to 0, entirely skips the checks.
 * Defaults to n=10.
 */
OGDF_EXPORT extern int PCTREE_DEBUG_CHECK_FREQ;
#endif

/**
 * Functions that can be passed to PCTree::uniqueID()
 */
namespace uid_utils {
/**
 * Print the index of a node \p n.
 */
OGDF_EXPORT void nodeToID(std::ostream& os, PCNode* n, int pos);

/**
 * Print the position \p pos of a node \p n.
 */
OGDF_EXPORT void nodeToPosition(std::ostream& os, PCNode* n, int pos);

/**
 * Print the index of a node \p n if it is a leaf.
 */
OGDF_EXPORT void leafToID(std::ostream& os, PCNode* n, int pos);

/**
 * Print the position \p pos of a node \p n if it is a leaf.
 */
OGDF_EXPORT void leafToPosition(std::ostream& os, PCNode* n, int pos);

/**
 * Sort nodes by ascending index.
 */
OGDF_EXPORT bool compareNodesByID(PCNode* a, PCNode* b);
}

/**
 * A PC-tree represents a set of cyclic orders of its leaves by labeling its inner nodes as either P- or C-node
 * and allowing arbitrary permutations of the neighbors of P-nodes while only allowing flips of C-nodes.
 * The operation PCTree::makeConsecutive() updates the tree such that a set of leaves is consecutive in all represented cyclic orders.
 *
 * If you are using this implementation, please cite the following work:
 * \remark Simon D. Fink, Matthias Pfretzschner, and Ignaz Rutter. 2023. Experimental Comparison of PC-Trees and PQ-Trees. ACM J. Exp. Algorithmics 28, Article 1.10 (December 2023). https://doi.org/10.1145/3611653
 *
 * For more details, see also (open access):
 * \remark Simon D. Fink. 2024. Constrained Planarity Algorithms in Theory and Practice. Doctoral Thesis, University of Passau. https://doi.org/10.15475/cpatp.2024
 */
class OGDF_EXPORT PCTree {
	friend OGDF_EXPORT std::ostream&(operator<<)(std::ostream&, const ogdf::pc_tree::PCTree*);
	friend OGDF_EXPORT std::ostream&(operator<<)(std::ostream&, const ogdf::pc_tree::PCNode*);

	friend class PCNode;
	friend class PCTreeRegistry;
	friend class PCTreeForest;

public:
	struct Observer; // pre-declaration

private:
	// global
	PCTreeForest* m_forest = nullptr;

	// needs merge when trees merge
	size_t m_pNodeCount = 0;
	size_t m_cNodeCount = 0;
	IntrusiveList<PCNode> m_leaves;

	// needs special merge
	PCNode* m_rootNode = nullptr;

	// temp makeConsecutive variables
	size_t m_partialCount = 0;
	size_t m_terminalPathLength = 0;
	PCNode* m_firstPartial = nullptr;
	PCNode* m_lastPartial = nullptr;
	PCNode* m_apexCandidate = nullptr;
	bool m_apexCandidateIsFix = false;
	PCNode* m_apexTPPred2 = nullptr;

	// private
	bool m_externalForest = true;
	std::list<Observer*> m_observers;

public:
	/**
	 * Constructs a new PCTree.
	 * The tree will be part of an automatically managed internal forest and thus not allow merging with other trees.
	 */
	explicit PCTree() : PCTree(nullptr) {};

	/**
	 * Constructs a new PCTree that is part of \p forest and can thus be merged with any other PCTree of the same forest.
	 */
	explicit PCTree(PCTreeForest* forest) {
		if (forest) {
			m_forest = forest;
			m_externalForest = true;
		} else {
			m_forest = new PCTreeForest();
			m_externalForest = false;
		}
	};

	/**
	 * Convenience method generating a PCTree consisting of a single P-node with \p leafNum leaves, which are all copied to the optional list \p added.
	 * Automatically creates and manages a forest if \p forest is null.
	 */
	explicit PCTree(int leafNum, std::vector<PCNode*>* added = nullptr,
			PCTreeForest* forest = nullptr);

	/**
	 * Copy a PCTree.
	 * @param other The PCTree to copy.
	 * @param nodeMapping Will be assigned a mapping from the nodes of \p other to the newly created nodes in this tree.
	 * @param keep_ids Set to \c true to use the same node IDs as in \p other, otherwise consecutive IDs will be generated.
	 * @param forest Automatically create and manages a forest if \p forest is null.
	 */
	explicit PCTree(const PCTree& other, PCTreeNodeArray<PCNode*>& nodeMapping,
			bool keep_ids = false, PCTreeForest* forest = nullptr);

	/**
	 * Deserialize a PCTree from a string \p str as generated by operator<<(std::ostream&, const PCTree*) or PCTree::uniqueID().
	 * Automatically creates and manages a forest if \p forest is null.
	 */
	explicit PCTree(const std::string& str, bool keep_ids = false, PCTreeForest* forest = nullptr);

	virtual ~PCTree();

public:
	/**
	 * @name Node creation / destruction
	 */
	//! @{
	/**
	 * Create a new node.
	 *
	 * @param type The type of the node.
	 * @param parent Parent to attach this node to, may only be null if this tree is empty and thus has no root.
	 * @param id ID for the new node, or -1 to automatically use the next free one.
	 * @return the new node
	 */
	PCNode* newNode(PCNodeType type, PCNode* parent = nullptr, int id = -1);

	/**
	 * @copydoc destroyNode(PCNode* const&)
	 * @param node will be set to null afterwards
	 */
	void destroyNode(PCNode*& node) {
		destroyNode((PCNode* const&)node);
		node = nullptr;
	}

	/**
	 * Destroy a node.
	 * The node must be detached and may not be the root of this tree.
	 */
	void destroyNode(PCNode* const& node);

	/**
	 * Attach \p count leaves to P- or C-node \p parent and optionally store the new leaves in a vector \p added.
	 */
	void insertLeaves(int count, PCNode* parent, std::vector<PCNode*>* added = nullptr);

	/**
	 * Convert \p leaf into a P-node and attach \p leafCount new leaves to it.
	 */
	void replaceLeaf(int leafCount, PCNode* leaf, std::vector<PCNode*>* added = nullptr);

	/**
	 * Merge multiple leaves into a single one and return it.
	 *
	 * @param consecutiveLeaves The leaves that shall be merged.
	 * @param assumeConsecutive Set to \c true if you already called makeConsecutive() on the leaves.
	 * @return The entry of \p consecutiveLeaves into which all other leaves got merged.
	 */
	PCNode* mergeLeaves(const std::vector<PCNode*>& consecutiveLeaves,
			bool assumeConsecutive = false) {
		return mergeLeaves(consecutiveLeaves.begin(), consecutiveLeaves.end(), assumeConsecutive);
	}

	/**
	 * Merge multiple leaves into a single one and return it.
	 *
	 * @param begin, end Iterator range spanning the leaves that shall be merged.
	 * @param assumeConsecutive Set to \c true if you already called makeConsecutive() on the leaves.
	 * @return The entry of \p consecutiveLeaves into which all other leaves got merged.
	 */
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

		OGDF_HEAVY_ASSERT(checkValid());

		// Return the remaining leaf.
		return *back;
	}

	/**
	 * Remove a leaf and also any newly-introduced inner degree-2 or -1 nodes.
	 * Unlike destroyNode(), the leaf \p leaf must still be attached (so that its parent's degree
	 * can be checked)
	 * @sa destroyNode()
	 */
	void destroyLeaf(PCNode* leaf);

	/**
	 * Overwrite the stored root for this PC-tree.
	 *
	 * Note that the passed \p newRoot needs to be valid as root for this tree.
	 * @return The old node stored as root.
	 * @sa changeRoot()
	 */
	PCNode* setRoot(PCNode* newRoot);

	/**
	 * Change the orientation of edges such that \p newRoot becomes the root of the tree.
	 *
	 * @return The previous root.
	 * @sa setRoot()
	 */
	PCNode* changeRoot(PCNode* newRoot);

	/**
	 * Change the type of a node and update all its registrations.
	 * @return the previous type of the node.
	 */
	PCNodeType changeNodeType(PCNode* node, PCNodeType newType);

	/**
	 * Insert tree \p tree into this tree at node \p at.
	 *
	 * Both trees need to be part of the same forest. All observers of \p tree will be moved to this tree.
	 * If \p at is a leaf, it will be replaced by \p tree, otherwise \p tree will be appended as child of \p at.
	 */
	void insertTree(PCNode* at, PCTree* tree);

private:
	void unregisterNode(PCNode* node);

	void registerNode(PCNode* node);

	//! @}

public:
	/**
	 * @name Restrictions
	 */
	//! @{
	/**
	 * @return \c true if calling makeConsecutive() with \p size leaves never requires changes to the tree.
	 *   This is the case for \p size values 0, 1, getLeafCount() - 1, and getLeafCount().
	 * @sa ogdf::pc_tree::isTrivialRestriction()
	 */
	bool isTrivialRestriction(int size) const;

	bool makeConsecutive(std::initializer_list<PCNode*> consecutiveLeaves) {
		return makeConsecutive(consecutiveLeaves.begin(), consecutiveLeaves.end());
	}

	bool makeConsecutive(const std::vector<PCNode*>& consecutiveLeaves) {
		return makeConsecutive(consecutiveLeaves.begin(), consecutiveLeaves.end());
	}

	/**
	 * Make the leaves contained in the range denoted by iterators \p begin (inclusive) to
	 * \p end (exclusive) consecutive in all represented orders.
	 * This is equivalent to calling resetTempData() and markFull(It, It, std::vector<PCNode*>*) followd by makeFullNodesConsecutive().
	 * @return \c true if the update was successful, \c false if the leaves cannot be made
	 *    consecutive and the tree was left unchanged.
	 */
	template<typename It>
	bool makeConsecutive(It begin, It end) {
		FullLeafIter iter = [&begin, &end]() { return NextFullLeaf<It>(begin, end); };
		for (auto obs : m_observers) {
			obs->makeConsecutiveCalled(*this, iter);
		}

		OGDF_HEAVY_ASSERT(checkValid());
		resetTempData();

#ifdef OGDF_DEBUG
		for (auto it = begin; it != end; ++it) {
			PCNode* leaf = *it;
			OGDF_ASSERT(leaf);
			OGDF_ASSERT(leaf->isLeaf());
			OGDF_ASSERT(leaf->m_forest == m_forest);
		}
#endif
		if (isTrivialRestriction(end - begin)) {
			for (auto obs : m_observers) {
				obs->makeConsecutiveDone(*this, Observer::Stage::Trivial, true);
			}
			return true;
		}

		// PC_PROFILE_ENTER(1, "label");
		markFull(begin, end);
		// PC_PROFILE_EXIT(1, "label");

		return makeFullNodesConsecutive();
	}

	/**
	 * Reset all makeConsecutive()-related temporary information, especially which leaves are full (should be made consecutive).
	 */
	void resetTempData() {
		m_forest->m_timestamp++;
		m_firstPartial = m_lastPartial = nullptr;
		m_partialCount = 0;
		m_apexCandidate = nullptr;
		m_apexCandidateIsFix = false;
		m_terminalPathLength = 0;
		m_apexTPPred2 = nullptr;
	}

	/**
	 * Only marks leaves full, does not update partial/full info of parents.
	 * Use markFull(It, It, std::vector<PCNode*>*) to also update parents for use with makeFullNodesConsecutive().
	 */
	template<typename It>
	void markLeavesFull(It begin, It end) {
		for (auto it = begin; it != end; ++it) {
			PCNode* leaf = *it;
			OGDF_ASSERT(leaf);
			OGDF_ASSERT(leaf->isLeaf());
			OGDF_ASSERT(leaf->m_forest == m_forest);
			leaf->setLabel(NodeLabel::Full);
		}
	}

	/**
	 * Marks the leaves contained in the range denoted by iterators \p begin (inclusive) to
	 * \p end (exclusive) full, that is marked for being made consecutive in all represented orders.
	 * Also propagates the markings to parents, which is required for makeFullNodesConsecutive().
	 */
	template<typename It>
	void markFull(It begin, It end, std::vector<PCNode*>* fullNodeOrder = nullptr) {
		// Attention: This method no longer uses a queue to defer processing of partial/full parents to after
		// all leaves are done, but now directly make parents full if all of their children are full.
		if (fullNodeOrder != nullptr) {
			fullNodeOrder->reserve(m_cNodeCount + m_pNodeCount);
		}

		for (auto it = begin; it != end; ++it) {
			PCNode* full_parent = markFull(*it, fullNodeOrder);
			while (full_parent != nullptr) {
				full_parent = markFull(full_parent, fullNodeOrder);
			}
		}
	}

	/**
	 * Updates the tree to make all leaves marked as full consecutive in all represented orders.
	 * Requires labels of parents to be correctly set by markFull(It, It, std::vector<PCNode*>*).
	 * @return \c true if the update was successful, \c false if the leaves cannot be made
	 *    consecutive and the tree was left unchanged.
	 */
	bool makeFullNodesConsecutive();

private:
	/* see the paper for more info on how the update works with the following methods */

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

	//! @}

public:
	/**
	 * @name Intersect
	 */
	//! @{
	/**
	 * Intersect the restrictions represented by this tree with those represented by \p other,
	 * given a bijection \p mapping from the leaves of \p other to this trees' leaves.
	 * @return \c true if the intersection is non-empty and this now represented by this tree.
	 *   Otherwise, the intersection is empty and the state of this tree is undefined.
	 */
	bool intersect(PCTree& other, PCTreeNodeArray<PCNode*>& mapping);

private:
	bool findNodeRestrictions(PCTree& applyTo, PCTreeNodeArray<PCNode*>& mapping,
			PCTreeNodeArray<std::vector<PCNode*>>& blockNodes,
			PCTreeNodeArray<std::vector<PCNode*>>& subtreeNodes,
			PCTreeNodeArray<PCNode*>& leafPartner, PCTreeNodeArray<bool>& isFront);

	void restoreSubtrees(PCTreeNodeArray<std::vector<PCNode*>>& blockNodes,
			PCTreeNodeArray<std::vector<PCNode*>>& subtreeNodes,
			PCTreeNodeArray<PCNode*>& leafPartner, PCTreeNodeArray<bool>& isFront);

	//! @}

public:
	/**
	 * @name Getters
	 */
	//! @{
	operator const PCTreeRegistry&() const { return m_forest->m_nodeArrayRegistry; }

	//! The PCTreeForest this PCTree belongs to, or nullptr.
	[[nodiscard]] PCTreeForest* getForest() const { return m_forest; }

	//! Whether this PCTree allows all orders (consists of a single P-node).
	[[nodiscard]] bool isTrivial() const;

	[[nodiscard]] const IntrusiveList<PCNode>& getLeaves() const { return m_leaves; }

	[[nodiscard]] size_t getLeafCount() const { return m_leaves.size(); };

	[[nodiscard]] size_t getPNodeCount() const { return m_pNodeCount; }

	[[nodiscard]] size_t getCNodeCount() const { return m_cNodeCount; }

	[[nodiscard]] PCNode* getRootNode() const { return m_rootNode; }

	/**
	 * @return the length of the terminal path in the last update operation
	 */
	int getTerminalPathLength() const { return m_terminalPathLength; }

	//! An iterable through all nodes of this PCTree.
	FilteringPCTreeDFS allNodes() const { return FilteringPCTreeDFS(*this, m_rootNode); }

	//! An iterable through all inner (non-leaf) nodes of this PCTree.
	FilteringPCTreeDFS innerNodes() const {
		return FilteringPCTreeDFS(*this, m_rootNode, [](PCNode* node) { return !node->isLeaf(); });
	}

	//! Store the order of leaves currently represented by this tree in \p container.
	template<typename Container>
	void currentLeafOrder(Container& container) const {
		for (PCNode* leaf : FilteringPCTreeDFS(*this, m_rootNode)) {
			if (leaf->isLeaf()) {
				container.push_back(leaf);
			}
		}
		OGDF_ASSERT(container.size() == m_leaves.size());
	}

	std::vector<PCNode*> currentLeafOrder() const {
		std::vector<PCNode*> container;
		currentLeafOrder(container);
		return container;
	}

	//! Validity check for debugging assertions.
	bool checkValid(bool allow_small_deg = true) const;

	//! Check whether the order \p order is represented by this tree.
	bool isValidOrder(const std::vector<PCNode*>& order) const;

	//! Get a graphical representation of this tree as Graph.
	void getTree(ogdf::Graph& tree, ogdf::GraphAttributes* g_a,
			PCTreeNodeArray<ogdf::node>& pc_repr, ogdf::NodeArray<PCNode*>* g_repr = nullptr,
			bool mark_full = false, bool show_sibs = false) const;

	/**
	 * Get a list of all cyclic restrictions used to generate this tree.
	 * If a \p fixedLeaf was given, the restrictions will linear with none of them containing \p fixedLeaf.
	 */
	void getRestrictions(std::vector<std::vector<PCNode*>>& restrictions,
			PCNode* fixedLeaf = nullptr) const;

	//! Calculate the number of cyclic orders represented by this tree.
	template<typename R>
	R possibleOrders() const {
		R orders(1);
		for (PCNode* node : innerNodes()) {
			if (node->getNodeType() == PCNodeType::CNode) {
				orders *= 2;
			} else {
				R children(node->getChildCount());
				if (node == m_rootNode) {
					children -= 1; // don't count circular shifts
				}
				orders *= factorial(children);
			}
		}
		return orders;
	}

	/**
	 * Print a deterministic and unique representation of this PCTree to \p os.
	 * Unique node IDs and a deterministic order of nodes' children is generated using \p printNode and \p compareNodes, respectively.
	 * @sa ogdf::pc_tree::uid_utils
	 */
	std::ostream& uniqueID(std::ostream& os,
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

	//! @}

public:
	/**
	 * @name  Observers
	 */
	//! @{
	template<typename It>
	struct NextFullLeaf {
		It m_it;
		It m_end;

		NextFullLeaf(It it, It an_end) : m_it(it), m_end(an_end) { }

		PCNode* operator()() {
			if (m_it == m_end) {
				return nullptr;
			}
			PCNode* n = *m_it;
			++m_it;
			return n;
		}
	};

	using FullLeafIter = std::function<std::function<PCNode*()>()>;

	//! Interface for Observers that can be notified of all changes made to the tree during an update.
	struct Observer {
		enum class Stage { Trivial, NoPartials, InvalidTP, SingletonTP, Done };

		virtual void onNodeCreate([[maybe_unused]] PCNode* node) {};

		virtual void makeConsecutiveCalled([[maybe_unused]] PCTree& tree,
				[[maybe_unused]] FullLeafIter consecutiveLeaves) {};

		virtual void labelsAssigned([[maybe_unused]] PCTree& tree,
				[[maybe_unused]] PCNode* firstPartial, [[maybe_unused]] PCNode* lastPartial,
				[[maybe_unused]] int partialCount) {};

		virtual void terminalPathFound([[maybe_unused]] PCTree& tree, [[maybe_unused]] PCNode* apex,
				[[maybe_unused]] PCNode* apexTPPred2, [[maybe_unused]] int terminalPathLength) {};

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
				[[maybe_unused]] PCNode* apexCandidate, [[maybe_unused]] PCNode* central,
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
		m_observers.push_back(observer);
		return --m_observers.end();
	}

	void removeObserver(std::list<Observer*>::const_iterator it) { m_observers.erase(it); }

	void removeObserver(Observer* observer) { m_observers.remove(observer); }

	//! @}
};

int PCTreeRegistry::keyToIndex(PCNode* key) { return key->index(); }

}
