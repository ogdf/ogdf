/** \file
 * \brief Implementation for ogdf::pc_tree::PCTree update/restriction methods
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

#include <ogdf/basic/basic.h>
#include <ogdf/basic/pctree/PCEnum.h>
#include <ogdf/basic/pctree/PCNode.h>
#include <ogdf/basic/pctree/PCTree.h>
#include <ogdf/basic/pctree/PCTreeIterators.h>
#include <ogdf/basic/pctree/util/IntrusiveList.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <list>
#include <string>
#include <utility>
#include <vector>

#ifdef LIKWID_PERFMON

#	include <likwid.h>

#	define PC_PROFILE_ENTER(level, msg) PC_PROFILE_ENTER_##level(msg)
#	define PC_PROFILE_EXIT(level, msg) PC_PROFILE_EXIT_##level(msg)

#	ifdef PC_PROFILE_LEVEL_1
#		define PC_PROFILE_ENTER_1(msg) likwid_markerStartRegion(msg)
#		define PC_PROFILE_EXIT_1(msg) likwid_markerStopRegion(msg)
#	else
#		define PC_PROFILE_ENTER_1(msg)
#		define PC_PROFILE_EXIT_1(msg)
#	endif

#	ifdef PC_PROFILE_LEVEL_2
#		define PC_PROFILE_ENTER_2(msg) likwid_markerStartRegion(msg)
#		define PC_PROFILE_EXIT_2(msg) likwid_markerStopRegion(msg)
#	else
#		define PC_PROFILE_ENTER_2(msg)
#		define PC_PROFILE_EXIT_2(msg)
#	endif

#else
#	define PC_PROFILE_ENTER(level, msg)
#	define PC_PROFILE_EXIT(level, msg)
#endif

using namespace ogdf::pc_tree;

#ifdef UFPC_DEBUG
#	define log   \
		if (true) \
		std::cout

#	include <ogdf/fileformats/GraphIO.h>
#	include <ogdf/misclayout/CircularLayout.h>
#	include <ogdf/tree/TreeLayout.h>

void dump(PCTree& T, const std::string& name) {
	Graph G;
	GraphAttributes GA(G, GraphAttributes::all);
	PCTreeNodeArray<ogdf::node> pc_repr(T, nullptr);
	T.getTree(G, &GA, pc_repr, nullptr, true);
	CircularLayout cl;
	cl.call(GA);
	GraphIO::write(GA, name + "-cl.svg");
	TreeLayout tl;
	G.reverseAllEdges();
	tl.call(GA);
	GraphIO::write(GA, name + "-tl.svg");
}

#else
#	define log    \
		if (false) \
		std::cout

void dump(PCTree& T, const std::string& name) { }

#endif


namespace ogdf::pc_tree {

bool isTrivialRestriction(int restSize, int leafCount) {
	return restSize <= 1 || restSize >= leafCount - 1;
}

bool PCTree::isTrivialRestriction(int size) const {
	return ::isTrivialRestriction(size, getLeafCount());
}

int factorial(int n) { return (int)std::tgamma(n + 1); }

void PCTree::LoggingObserver::makeConsecutiveCalled(PCTree& tree, FullLeafIter consecutiveLeaves) {
	log << "Tree " << tree << " with consecutive leaves [";
	auto it = consecutiveLeaves();
	for (auto l = it(); l != nullptr; l = it()) {
		log << l->index() << ", ";
	}
	log << "]" << std::endl;
	dump(tree, "before");
}

void PCTree::LoggingObserver::labelsAssigned(PCTree& tree, PCNode* firstPartial,
		PCNode* lastPartial, int partialCount) {
	log << "Found " << partialCount << " partial nodes: [";
#ifdef OGDF_DEBUG
	int found = 0;
#endif
	for (PCNode* node = firstPartial; node != nullptr; node = node->tempInfo().nextPartial) {
		log << node << ", ";
#ifdef OGDF_DEBUG
		found++;
#endif
	}
	log << "]" << std::endl;
	OGDF_ASSERT(partialCount == found);
	dump(tree, "labelled");
}

void PCTree::LoggingObserver::terminalPathFound(PCTree& tree, PCNode* apex, PCNode* apexTPPred2,
		int terminalPathLength) {
	log << "Apex is " << apex->getLabel() << " " << apex << ", TP length is " << terminalPathLength
		<< ", first path is [";
	int length = 1;
	for (PCNode* cur = apex->constTempInfo().tpPred; cur != nullptr;
			cur = cur->constTempInfo().tpPred) {
		log << cur << ", ";
		length++;
	}
	log << "], second path is [";
	for (PCNode* cur = apexTPPred2; cur != nullptr; cur = cur->constTempInfo().tpPred) {
		log << cur << ", ";
		length++;
	}
	log << "] (effective length " << length << ")" << std::endl;
	OGDF_ASSERT(terminalPathLength == length);
}

void PCTree::LoggingObserver::centralCreated(PCTree& tree, PCNode* central) {
	log << "Tree before merge into central " << central->index() << ": " << tree << std::endl;
	dump(tree, "before-merge");
}

void PCTree::LoggingObserver::beforeMerge(PCTree& tree, int count, PCNode* tpNeigh) {
	PCNode* central = tpNeigh->getParent();
	PCNode* fullNeigh = central->getFullNeighInsertionPoint(tpNeigh);
	log << "Merging child #" << count << " " << tpNeigh->index() << " next to "
		<< fullNeigh->index() << " into parent " << central->index() << std::endl
		<< "\t" << tree << std::endl;
	dump(tree, "merge-" + std::to_string(count) + "-" + std::to_string(tpNeigh->index()));
}

void PCTree::LoggingObserver::makeConsecutiveDone(PCTree& tree, Stage stage, bool success) {
	log << "Resulting tree after stage ";
	switch (stage) {
	case Stage::Trivial:
		log << "Trivial";
		break;
	case Stage::NoPartials:
		log << "NoPartials";
		break;
	case Stage::InvalidTP:
		log << "InvalidTP";
		break;
	case Stage::SingletonTP:
		log << "SingletonTP";
		break;
	case Stage::Done:
		log << "Done";
		break;
	default:
		log << "???";
		break;
	}
	log << ", restriction " << (success ? "successful" : "invalid") << ": " << tree << std::endl;
}

extern int PCTREE_DEBUG_CHECK_CNT;

bool PCTree::makeFullNodesConsecutive() {
	if (m_firstPartial == nullptr) {
		OGDF_ASSERT(m_lastPartial == nullptr);
		OGDF_ASSERT(m_partialCount == 0);
		for (auto obs : m_observers) {
			obs->makeConsecutiveDone(*this, Observer::Stage::NoPartials, true);
		}
		return true;
	}
	OGDF_ASSERT(m_lastPartial != nullptr);
	OGDF_ASSERT(m_partialCount > 0);
	for (auto obs : m_observers) {
		obs->labelsAssigned(*this, m_firstPartial, m_lastPartial, m_partialCount);
	}

	PC_PROFILE_ENTER(1, "find_tp");
	bool find_tp = findTerminalPath();
	PC_PROFILE_EXIT(1, "find_tp");
	if (!find_tp) {
		for (auto obs : m_observers) {
			obs->makeConsecutiveDone(*this, Observer::Stage::InvalidTP, false);
		}
		return false;
	}
	OGDF_ASSERT(m_apexCandidate != nullptr);
	OGDF_ASSERT(m_apexCandidateIsFix == true);
	for (auto obs : m_observers) {
		obs->terminalPathFound(*this, m_apexCandidate, m_apexTPPred2, m_terminalPathLength);
	}

	PC_PROFILE_ENTER(1, "update_tp");
	if (m_terminalPathLength == 1) {
		OGDF_ASSERT(m_apexCandidate->tempInfo().tpPred == nullptr);
		updateSingletonTerminalPath();
		PC_PROFILE_EXIT(1, "update_tp");
		for (auto obs : m_observers) {
			obs->makeConsecutiveDone(*this, Observer::Stage::SingletonTP, true);
		}
		return true;
	}
	OGDF_ASSERT(m_apexCandidate->tempInfo().tpPred != nullptr);
	PC_PROFILE_ENTER(2, "update_tp_central");
	PCNode* central = createCentralNode();
	PC_PROFILE_EXIT(2, "update_tp_central");
	for (auto obs : m_observers) {
		obs->centralCreated(*this, central);
	}

	PCNode::TempInfo& ctinfo = central->tempInfo();
#ifdef OGDF_DEBUG
	size_t merged =
#endif
			updateTerminalPath(central, ctinfo.tpPred);
	if (m_apexTPPred2 != nullptr) {
#ifdef OGDF_DEBUG
		merged +=
#endif
				updateTerminalPath(central, m_apexTPPred2);
	}
	OGDF_ASSERT(merged == m_terminalPathLength - 1);
	PC_PROFILE_EXIT(1, "update_tp");

	for (auto obs : m_observers) {
		obs->makeConsecutiveDone(*this, Observer::Stage::Done, true);
	}
#ifdef OGDF_HEAVY_DEBUG
	OGDF_HEAVY_ASSERT(checkValid());
	if (PCTREE_DEBUG_CHECK_FREQ != 0 && PCTREE_DEBUG_CHECK_CNT % PCTREE_DEBUG_CHECK_FREQ == 0) {
		std::list<PCNode*> order;
		currentLeafOrder(order);
		while (!order.front()->isFull()) {
			order.push_back(order.front());
			order.pop_front();
		}
		while (order.back()->isFull()) {
			order.push_front(order.back());
			order.pop_back();
		}
		OGDF_ASSERT(order.front()->isFull());
		OGDF_ASSERT(order.size() == m_leaves.size());
		while (order.front()->isFull()) {
			order.pop_front();
		}
		while (!order.empty()) {
			OGDF_ASSERT(!order.front()->isFull());
			order.pop_front();
		}
	}
#endif
	return true;
}

void PCTree::addPartialNode(PCNode* partial) {
	OGDF_ASSERT(partial->tempInfo().predPartial == nullptr);
	OGDF_ASSERT(partial->tempInfo().nextPartial == nullptr);
	partial->tempInfo().predPartial = m_lastPartial;
	if (m_firstPartial == nullptr) {
		m_firstPartial = partial;
	} else {
		m_lastPartial->tempInfo().nextPartial = partial;
	}
	m_lastPartial = partial;
	m_partialCount++;
}

void PCTree::removePartialNode(PCNode* partial) {
	PCNode::TempInfo& tinfo = partial->tempInfo();
	OGDF_ASSERT((tinfo.predPartial == nullptr) == (m_firstPartial == partial));
	if (tinfo.predPartial == nullptr) {
		m_firstPartial = tinfo.nextPartial;
	} else {
		tinfo.predPartial->tempInfo().nextPartial = tinfo.nextPartial;
	}
	OGDF_ASSERT((tinfo.nextPartial == nullptr) == (m_lastPartial == partial));
	if (tinfo.nextPartial == nullptr) {
		m_lastPartial = tinfo.predPartial;
	} else {
		tinfo.nextPartial->tempInfo().predPartial = tinfo.predPartial;
	}
	tinfo.predPartial = tinfo.nextPartial = nullptr;
	m_partialCount--;
}

PCNode* PCTree::markFull(PCNode* full_node, std::vector<PCNode*>* fullNodeOrder) {
	OGDF_HEAVY_ASSERT(full_node->isValidNode(m_forest));
	if (!full_node->isLeaf()) {
		OGDF_ASSERT(full_node->tempInfo().fullNeighbors.size() == full_node->getDegree() - 1);
		OGDF_ASSERT(full_node->getLabel() == NodeLabel::Partial);
		removePartialNode(full_node);
	}
	full_node->setLabel(NodeLabel::Full);

	PCNode* partial_neigh = full_node->getParent();
	if (partial_neigh == nullptr || partial_neigh->isFull()) {
		// if we are the root or our parent node got full before us, we need to find our one non-full neighbor
		PC_PROFILE_ENTER(2, "label_process_neigh");
		PCNode* pred = nullptr;
		partial_neigh = full_node->m_child1;
		while (partial_neigh != nullptr && partial_neigh->isFull()) {
			proceedToNextSibling(pred, partial_neigh);
		}
		PC_PROFILE_EXIT(2, "label_process_neigh");
	}
	OGDF_ASSERT(partial_neigh != nullptr);
	OGDF_ASSERT(!partial_neigh->isFull());
	if (partial_neigh->isLeaf()) {
		// when a leaf becomes partial, all other leaves need to be full, which only happens during getRestriction
		OGDF_ASSERT(fullNodeOrder != nullptr);
		OGDF_ASSERT(fullNodeOrder->size() == getPNodeCount() + getCNodeCount());
		return nullptr;
	}

	size_t fullNeighborCounter = partial_neigh->addFullNeighbor(full_node);
	OGDF_ASSERT(fullNeighborCounter >= 1);
	OGDF_ASSERT(fullNeighborCounter + 1 <= partial_neigh->getDegree());
	if (fullNeighborCounter == 1) {
		OGDF_ASSERT(partial_neigh->getLabel() == NodeLabel::Unknown);
		partial_neigh->setLabel(NodeLabel::Partial);
		addPartialNode(partial_neigh);
	} else {
		OGDF_ASSERT(partial_neigh->getLabel() == NodeLabel::Partial);
	}
	if (fullNeighborCounter == partial_neigh->getDegree() - 1) {
		if (fullNodeOrder != nullptr) {
			fullNodeOrder->push_back(partial_neigh);
		}
		return partial_neigh;
	} else {
		return nullptr;
	}
}

bool PCTree::checkTPPartialCNode(PCNode* node) {
	// check that C node's full neighbors are consecutive
	PCNode::TempInfo& tinfo = node->tempInfo();
	if (tinfo.ebEnd1 == nullptr) {
		PC_PROFILE_ENTER(2, "find_tp_cnode");
		PCNode* fullChild = tinfo.fullNeighbors.front();
		PCNode* sib1 = node->getNextNeighbor(nullptr, fullChild);
		PCNode* sib2 = node->getNextNeighbor(sib1, fullChild);
		size_t count = 1;
		count += findEndOfFullBlock(node, fullChild, sib1, tinfo.fbEnd1, tinfo.ebEnd1);
		count += findEndOfFullBlock(node, fullChild, sib2, tinfo.fbEnd2, tinfo.ebEnd2);
		PC_PROFILE_EXIT(2, "find_tp_cnode");
		if (count != tinfo.fullNeighbors.size()) {
			log << "C-node's full-block isn't consecutive, abort!" << std::endl;
			return false;
		}
	}
	OGDF_ASSERT(tinfo.ebEnd1 != nullptr);
	OGDF_ASSERT(tinfo.ebEnd2 != nullptr);
	OGDF_ASSERT(tinfo.fbEnd1 != nullptr);
	OGDF_ASSERT(tinfo.fbEnd2 != nullptr);
	if (tinfo.tpPred != nullptr) {
		if (tinfo.tpPred != tinfo.ebEnd1 && tinfo.tpPred != tinfo.ebEnd2) {
			log << "C-node's TP pred is not adjacent to its empty block, abort!" << std::endl;
			return false;
		}
	}
	if (node == m_apexCandidate && m_apexTPPred2 != nullptr) {
		if (m_apexTPPred2 != tinfo.ebEnd1 && m_apexTPPred2 != tinfo.ebEnd2) {
			log << "C-node's TP second pred is not adjacent to its empty block, abort!" << std::endl;
			return false;
		}
	}
	return true;
}

bool PCTree::findTerminalPath() {
	while (m_firstPartial != nullptr) {
		OGDF_ASSERT(m_lastPartial != nullptr);
		OGDF_ASSERT(m_partialCount > 0);

		PCNode* node = m_firstPartial;
		removePartialNode(node);
		PCNode* parent = node->getParent();
		PCNode::TempInfo& tinfo = node->tempInfo();
		OGDF_ASSERT(!node->isFull());
		log << "Processing " << node->getLabel() << " " << node;
		if (parent != nullptr) {
			log << " (" << parent->getLabel() << ")";
		}
		log << ", current TP length is " << m_terminalPathLength << ": ";

		if (node->m_nodeType == PCNodeType::CNode && node->getLabelUnchecked() == NodeLabel::Partial
				&& !checkTPPartialCNode(node)) {
			return false;
		}
		if (node->getLabelUnchecked() == NodeLabel::Partial) {
			tinfo.tpPartialPred = node;
			tinfo.tpPartialHeight = 0;
		}

		OGDF_ASSERT((parent == nullptr) == (node == m_rootNode));
		if (tinfo.tpSucc != nullptr) {
			// we never process a node twice, proceeding to the next entry if we detect this case
			log << "dupe!" << std::endl;
			continue;
		}

		if (m_firstPartial == nullptr && m_apexCandidate == nullptr) {
			// we can stop early if the queue size reached 1 right before we removed the current node,
			// but we have not found an apex, marking the node the remaining arc points to as apex candidate
			log << "early stop I-shaped!" << std::endl;
#ifdef OGDF_DEBUG
			bool s =
#endif
					setApexCandidate(node, false);
			OGDF_ASSERT(s);
		} else if (m_firstPartial == nullptr && m_apexCandidate != nullptr && m_apexTPPred2 != nullptr
				&& tinfo.tpPartialPred == m_apexCandidate->tempInfo().tpPartialPred) {
			log << "early stop A-shaped!" << std::endl;
			OGDF_ASSERT(m_apexCandidate->tempInfo().tpPartialHeight <= tinfo.tpPartialHeight);
			OGDF_ASSERT(m_apexCandidateIsFix);
#ifdef OGDF_DEBUG
			// setApexCandidate here should make no changes, but just re-do the checks from this branch and keep the old apex
			PCNode* oldApex = m_apexCandidate;
			bool s = setApexCandidate(node, false);
			OGDF_ASSERT(s);
			OGDF_ASSERT(m_apexCandidate == oldApex);
#endif
		} else if (node == m_rootNode || parent->getLabel() == NodeLabel::Full) {
			// we can't ascend from the root node or if our parent is full
			log << "can't ascend from root / node with full parent!" << std::endl;
			if (!setApexCandidate(node, false)) {
				return false;
			}
		} else {
			// check that we can ascend through a C-Node (P-Nodes always work)
			if (node->m_nodeType == PCNodeType::CNode) {
				if (node->getLabelUnchecked() == NodeLabel::Empty) {
					if (!node->isChildOuter(tinfo.tpPred)) {
						log << "can't ascend from empty C-Node where TP pred is not adjacent to parent!"
							<< std::endl;
						if (!setApexCandidate(node, false)) {
							return false;
						}
						continue;
					} else {
						tinfo.ebEnd1 = tinfo.fbEnd2 = tinfo.tpPred;
						tinfo.ebEnd2 = tinfo.fbEnd1 = parent;
					}
				} else {
					OGDF_ASSERT(node->getLabelUnchecked() == NodeLabel::Partial);
					if (!node->isChildOuter(tinfo.fbEnd1) && !node->isChildOuter(tinfo.fbEnd2)) {
						log << "can't ascend from partial C-Node where full block is not adjacent to parent!"
							<< std::endl;
						if (!setApexCandidate(node, false)) {
							return false;
						}
						continue;
					}
				}
			}

			// ascend to parent
			m_terminalPathLength++;
			tinfo.tpSucc = parent;
			PCNode::TempInfo& parent_tinfo = parent->tempInfo();
			if (parent_tinfo.tpPred == nullptr) {
				// attach current path to parent
				parent_tinfo.tpPred = node;
				if (parent->getLabelUnchecked() != NodeLabel::Partial) { // parent might be empty or full
					parent_tinfo.tpPartialPred = tinfo.tpPartialPred;
					parent_tinfo.tpPartialHeight = tinfo.tpPartialHeight + 1;
					OGDF_ASSERT(parent_tinfo.tpPartialPred != nullptr);
					addPartialNode(parent);
					log << "proceed to non-partial parent (whose partial height is "
						<< parent_tinfo.tpPartialHeight << ")" << std::endl;
				} else {
					OGDF_ASSERT(parent->getLabelUnchecked() == NodeLabel::Partial);
					log << "partial parent is already queued" << std::endl;
				}
			} else if (parent_tinfo.tpPred != node) {
				log << "parent is A-shaped apex!" << std::endl;
				if (!setApexCandidate(parent, true)) {
					return false;
				}
				if (m_apexTPPred2 != nullptr && m_apexTPPred2 != node) {
					log << "Conflicting apexTPPred2!" << std::endl;
					return false;
				}
				m_apexTPPred2 = node;
				if (parent->m_nodeType == PCNodeType::CNode
						&& parent->getLabelUnchecked() == NodeLabel::Empty) {
					if (!parent->areNeighborsAdjacent(parent_tinfo.tpPred, m_apexTPPred2)) {
						log << "Apex is empty C-Node, but partial predecessors aren't adjacent!"
							<< std::endl;
						return false;
					} else {
						parent_tinfo.ebEnd1 = parent_tinfo.fbEnd2 = parent_tinfo.tpPred;
						parent_tinfo.ebEnd2 = parent_tinfo.fbEnd1 = m_apexTPPred2;
					}
				}
			}

			// if the parent is a partial C-node, it might have been processed before we were added as its tpPred(2), so re-run the checks
			if (parent->m_nodeType == PCNodeType::CNode
					&& parent->getLabelUnchecked() == NodeLabel::Partial) {
				if (!checkTPPartialCNode(parent)) {
					return false;
				}
			}
		}
	}
	OGDF_ASSERT(m_lastPartial == nullptr);
	OGDF_ASSERT(m_partialCount == 0);
	OGDF_ASSERT(m_apexCandidate != nullptr);
	if (!m_apexCandidateIsFix) {
		if (m_apexCandidate->getLabel() != NodeLabel::Partial) {
			log << "Backtracking from " << m_apexCandidate << " to "
				<< m_apexCandidate->tempInfo().tpPartialPred << ", TP length "
				<< m_terminalPathLength << "-" << m_apexCandidate->tempInfo().tpPartialHeight << "=";
			m_terminalPathLength -= m_apexCandidate->tempInfo().tpPartialHeight;
			m_apexCandidate = m_apexCandidate->tempInfo().tpPartialPred;
			log << m_terminalPathLength << std::endl;
			OGDF_ASSERT(m_apexCandidate->tempInfo().tpPartialHeight
					<= m_terminalPathLength); // make sure we didn't ascent too far
		}
		m_apexCandidateIsFix = true;
	}
	if (m_apexCandidate->m_nodeType == PCNodeType::CNode
			&& m_apexCandidate->getLabel() == NodeLabel::Empty) {
		OGDF_ASSERT(m_apexCandidate->tempInfo().tpPred != nullptr);
		OGDF_ASSERT(m_apexTPPred2 != nullptr);
	}
	m_apexCandidate->tempInfo().tpSucc = nullptr;
	return true;
}

PCNode* PCTree::createCentralNode() {
	PCNode::TempInfo& atinfo = m_apexCandidate->tempInfo();
	PCNode::TempInfo* ctinfo;
	PCNode* parent = m_apexCandidate->getParent();
	PCNode* central;
	PCNode* tpStubApex1 = atinfo.tpPred;
	PCNode* tpStubApex2 = m_apexTPPred2;
	OGDF_ASSERT(tpStubApex1 != nullptr);
	OGDF_ASSERT(m_apexCandidate->getLabelUnchecked() != NodeLabel::Unknown || tpStubApex2 != nullptr);

	if (m_apexCandidate->m_nodeType == PCNodeType::PNode) {
		// calculate the sizes of all sets
		bool isParentFull = false;
		if (parent != nullptr) {
			OGDF_ASSERT(parent->getLabel() != NodeLabel::Partial);
			isParentFull = parent->isFull();
		}
		size_t fullNeighbors = atinfo.fullNeighbors.size();
		size_t partialNeighbors = 1; // tpStubApex1
		if (tpStubApex2 != nullptr) {
			partialNeighbors = 2;
		}
		OGDF_ASSERT(m_apexCandidate->getDegree() >= fullNeighbors + partialNeighbors);
		size_t emptyNeighbors = m_apexCandidate->getDegree() - fullNeighbors - partialNeighbors;
		log << "Apex is " << m_apexCandidate->getLabel() << " " << m_apexCandidate
			<< " with neighbors: full=" << fullNeighbors << ", partial=" << partialNeighbors
			<< ", empty=" << emptyNeighbors << std::endl;
		if (parent != nullptr) {
			log << "Parent is " << parent->getLabel() << " " << parent << std::endl;
		}

		// isolate the apex candidate
		tpStubApex1->detach();
		if (tpStubApex2 != nullptr) {
			tpStubApex2->detach();
		}

		// create the new central node
		central = newNode(PCNodeType::CNode);
		ctinfo = &central->tempInfo();
		log << "New Central has index " << central->index();

		// first step of assembly: tpStubApex1 == apexCandidate.tempInfo().tpPred
		central->appendChild(tpStubApex1);
		tpStubApex1->tempInfo().replaceNeighbor(m_apexCandidate, central);

		// second step of assembly: fullNode
		PCNode* fullNode = nullptr;
		if (fullNeighbors == 1 && isParentFull) {
			m_apexCandidate->replaceWith(central);
			for (auto obs : m_observers) {
				obs->nodeReplaced(*this, m_apexCandidate, central);
			}
			fullNode = parent;
			log << ", full parent is node " << parent->index();
			OGDF_ASSERT(atinfo.fullNeighbors.front() == parent);
		} else if (fullNeighbors > 0) {
			fullNode = splitOffFullPNode(m_apexCandidate, isParentFull);
			if (isParentFull) {
				OGDF_ASSERT(!m_apexCandidate->isDetached());
				m_apexCandidate->replaceWith(fullNode);
				for (auto obs : m_observers) {
					obs->nodeReplaced(*this, m_apexCandidate, fullNode);
				}
				fullNode->appendChild(central);
			} else {
				central->appendChild(fullNode);
			}
			log << ", full " << (fullNeighbors == 1 ? "neigh" : "node") << " is " << fullNode;
			OGDF_ASSERT(fullNeighbors == 1 || fullNode->getDegree() == fullNeighbors + 1);
		}
		if (fullNode) {
			central->setLabelUnchecked(NodeLabel::Partial);
			ctinfo->fullNeighbors = {fullNode};
		} else {
			central->setLabelUnchecked(NodeLabel::Empty);
		}

		// third step of assembly: tpStubApex2 == apexTPPred2
		if (tpStubApex2 != nullptr) {
			central->appendChild(tpStubApex2, isParentFull);
			tpStubApex2->tempInfo().replaceNeighbor(m_apexCandidate, central);
		}

		// fourth step of assembly: apexCandidate (the empty node)
		PCNode* emptyNode = nullptr;
		if (emptyNeighbors == 1) {
			if (isParentFull || parent == nullptr) {
				emptyNode = m_apexCandidate->m_child1;
				emptyNode->detach();
				if (fullNode && tpStubApex2) {
					emptyNode->insertBetween(tpStubApex2, tpStubApex1);
				} else {
					central->appendChild(emptyNode, isParentFull);
				}
				log << ", empty neigh is " << emptyNode;
			} else {
				emptyNode = parent;
				m_apexCandidate->replaceWith(central);
				for (auto obs : m_observers) {
					obs->nodeReplaced(*this, m_apexCandidate, central);
				}
				log << ", empty parent is node " << parent->index();
			}
		} else if (emptyNeighbors > 1) {
			emptyNode = m_apexCandidate;
			if (isParentFull) {
				OGDF_ASSERT(fullNode);
				if (tpStubApex2) {
					m_apexCandidate->insertBetween(tpStubApex2, tpStubApex1);
				} else {
					central->appendChild(m_apexCandidate, true);
				}
			} else {
				m_apexCandidate->appendChild(central);
			}
			log << ", empty node is " << m_apexCandidate;
			OGDF_ASSERT(m_apexCandidate->getDegree() == emptyNeighbors + 1);
		}

		for (auto obs : m_observers) {
			obs->onApexMoved(*this, m_apexCandidate, central, parent);
		}

		if (emptyNeighbors <= 1) {
			if (m_apexCandidate == m_rootNode) {
				log << ", replaced root by central";
				m_rootNode = central;
			}
			destroyNode(m_apexCandidate);
		} else {
			OGDF_ASSERT(m_apexCandidate->getDegree() > 2);
			OGDF_ASSERT(m_apexCandidate->isDetached() == (m_rootNode == m_apexCandidate));
		}
		log << std::endl;
		OGDF_ASSERT(central->isDetached() == (m_rootNode == central));

		// assembly done, verify
#ifdef OGDF_DEBUG
		log << "Central " << central << " and neighbors [";
		std::vector<PCNode*> cn {tpStubApex1};
		log << "TP1 " << tpStubApex1->index() << ", ";
		if (fullNode) {
			log << "F " << fullNode->index() << ", ";
			cn.push_back(fullNode);
		}
		if (tpStubApex2) {
			log << "TP2 " << tpStubApex2->index() << ", ";
			cn.push_back(tpStubApex2);
		}
		if (emptyNode) {
			log << "E " << emptyNode->index();
			cn.push_back(emptyNode);
		}
		log << "]" << std::endl;
		OGDF_ASSERT(central->getDegree() == cn.size());
		OGDF_ASSERT(central->getDegree()
				== partialNeighbors + (fullNeighbors > 0 ? 1 : 0) + (emptyNeighbors > 0 ? 1 : 0));
		OGDF_ASSERT(central->getDegree() >= 3);
		std::list<PCNode*> actual_neighbors {
				central->neighbors().begin(), central->neighbors().end()};
		OGDF_ASSERT(actual_neighbors.size() == cn.size());
		for (size_t i = 0; actual_neighbors.front() != cn.front(); i++) {
			actual_neighbors.push_back(actual_neighbors.front());
			actual_neighbors.pop_front();
			OGDF_ASSERT(i < cn.size());
		}
		if (!std::equal(actual_neighbors.begin(), actual_neighbors.end(), cn.begin(), cn.end())) {
			actual_neighbors.push_back(actual_neighbors.front());
			actual_neighbors.pop_front();
			OGDF_ASSERT(std::equal(actual_neighbors.rbegin(), actual_neighbors.rend(), cn.begin(),
					cn.end()));
		}
#endif

		// update temporary pointers around central
		ctinfo->tpPred = tpStubApex1;
		if (fullNode) {
			ctinfo->ebEnd1 = tpStubApex1;
			ctinfo->fbEnd1 = ctinfo->fbEnd2 = fullNode;
			if (tpStubApex2) {
				ctinfo->ebEnd2 = tpStubApex2;
			} else {
				ctinfo->ebEnd2 = emptyNode;
			}
			OGDF_ASSERT(ctinfo->ebEnd2 != nullptr);
		} else {
			OGDF_ASSERT(tpStubApex2 != nullptr);
			ctinfo->ebEnd1 = ctinfo->fbEnd2 = tpStubApex1;
			ctinfo->ebEnd2 = ctinfo->fbEnd1 = tpStubApex2;
		}
	} else {
		central = m_apexCandidate;
		ctinfo = &atinfo;
		if (m_apexCandidate->getLabelUnchecked() == NodeLabel::Partial) {
			if (atinfo.ebEnd2 == tpStubApex1) {
				OGDF_ASSERT(atinfo.ebEnd1 == tpStubApex2 || tpStubApex2 == nullptr);
				std::swap(atinfo.ebEnd1, atinfo.ebEnd2);
				std::swap(atinfo.fbEnd1, atinfo.fbEnd2);
			}
		} else {
			OGDF_ASSERT(m_apexCandidate->getLabelUnchecked() == NodeLabel::Empty);
			OGDF_ASSERT(tpStubApex2 != nullptr);
			OGDF_ASSERT(central->areNeighborsAdjacent(tpStubApex1, tpStubApex2));
			atinfo.ebEnd1 = tpStubApex1;
			atinfo.fbEnd1 = tpStubApex2;
			atinfo.fbEnd2 = tpStubApex1;
			atinfo.ebEnd2 = tpStubApex2;
		}
	}

	OGDF_ASSERT(ctinfo->tpPred == tpStubApex1);
	OGDF_ASSERT(ctinfo->ebEnd1 == tpStubApex1);
	OGDF_ASSERT(central->getFullNeighInsertionPoint(tpStubApex1) == ctinfo->fbEnd1);
	if (m_apexTPPred2 != nullptr) {
		OGDF_ASSERT(ctinfo->ebEnd2 == tpStubApex2);
		OGDF_ASSERT(central->getFullNeighInsertionPoint(tpStubApex2) == ctinfo->fbEnd2);
	}
	return central;
}

int PCTree::updateTerminalPath(PCNode* central, PCNode* tpNeigh) {
	int count = 0;
	PCNode*& fullNeigh = central->getFullNeighInsertionPoint(tpNeigh);
	while (tpNeigh != nullptr) {
		PCNode::TempInfo& tinfo = tpNeigh->tempInfo();
		OGDF_ASSERT(central->areNeighborsAdjacent(tpNeigh, fullNeigh));
		OGDF_ASSERT(tinfo.tpSucc != nullptr);
		OGDF_ASSERT(tpNeigh->getLabelUnchecked() != NodeLabel::Full);
		OGDF_ASSERT(tpNeigh->getLabelUnchecked() != NodeLabel::Unknown || tinfo.tpPred != nullptr);
		for (auto obs : m_observers) {
			obs->beforeMerge(*this, count, tpNeigh);
		}
		PCNode* nextTPNeigh = tinfo.tpPred;
		PCNode* otherEndOfFullBlock;
		if (tpNeigh->m_nodeType == PCNodeType::PNode) {
			PC_PROFILE_ENTER(2, "update_tp_pnode");
			if (tpNeigh->getLabelUnchecked() == NodeLabel::Partial) {
				PCNode* fullNode = splitOffFullPNode(tpNeigh, false);
				fullNode->insertBetween(tpNeigh, fullNeigh);
				otherEndOfFullBlock = fullNeigh = fullNode;
				log << "\tFull child is " << fullNode << std::endl;
			} else {
				otherEndOfFullBlock = nullptr;
			}
			if (tinfo.tpPred != nullptr) {
				tinfo.tpPred->detach();
				log << "\ttpPred child is " << tinfo.tpPred << std::endl;
				tinfo.tpPred->insertBetween(fullNeigh, tpNeigh);
			}

			log << "\tThere are " << tpNeigh->m_childCount << " children left" << std::endl;

			for (auto obs : m_observers) {
				obs->whenPNodeMerged(*this, tpNeigh, tinfo.tpPred, otherEndOfFullBlock);
			}

			if (tpNeigh->m_childCount == 0) {
				tpNeigh->detach();
				for (auto obs : m_observers) {
					obs->nodeDeleted(*this, tpNeigh);
				}
				destroyNode(std::as_const(tpNeigh));
			} else if (tpNeigh->m_childCount == 1) {
				PCNode* child = tpNeigh->m_child1;
				child->detach();
				tpNeigh->replaceWith(child);
				for (auto obs : m_observers) {
					obs->nodeReplaced(*this, tpNeigh, child);
				}
				destroyNode(std::as_const(tpNeigh));
			}

			PC_PROFILE_EXIT(2, "update_tp_pnode");
		} else {
			PC_PROFILE_ENTER(2, "update_tp_cnode");

			auto print_tp_neigh = [](PCNode* node, int n) {
				//                cout << n << "   node: " << node << endl << "children: " << endl;
				//                for (auto child : node->children()) {
				//                    cout << child;
				//                    if (child->isLeaf()) {
				//                        cout << "    edge: "
				//                             << reinterpret_cast<adjEntry&>(child->leafUserData()[0]);
				//                    }
				//                    cout << endl;
				//                }
			};
			print_tp_neigh(tpNeigh, 1);

			OGDF_ASSERT(tpNeigh->m_nodeType == PCNodeType::CNode);
			PCNode* otherNeigh = central->getNextNeighbor(fullNeigh, tpNeigh);
			bool tpNeighSiblingsFlipped = false;
			if (tpNeigh->m_sibling1 == fullNeigh || tpNeigh->m_sibling2 == otherNeigh) {
				tpNeighSiblingsFlipped = true;
				std::swap(tpNeigh->m_sibling1, tpNeigh->m_sibling2);
			}
			OGDF_ASSERT(tpNeigh->m_sibling1 == otherNeigh
					|| (tpNeigh->m_sibling1 == nullptr && central->getParent() == otherNeigh)
					|| (tpNeigh->m_sibling1 == nullptr && central->getParent() == nullptr
							&& central->getOtherOuterChild(tpNeigh) == otherNeigh));
			OGDF_ASSERT(tpNeigh->m_sibling2 == fullNeigh
					|| (tpNeigh->m_sibling2 == nullptr && central->getParent() == fullNeigh)
					|| (tpNeigh->m_sibling2 == nullptr && central->getParent() == nullptr
							&& central->getOtherOuterChild(tpNeigh) == fullNeigh));

			print_tp_neigh(tpNeigh, 2);

			if (tpNeigh->m_child1 == tinfo.fbEnd1 || tpNeigh->m_child1 == tinfo.fbEnd2) {
				tpNeigh->flip();
			}
#ifdef OGDF_DEBUG
			PCNode* emptyOuterChild = tpNeigh->m_child1;
#endif
			PCNode* fullOuterChild = tpNeigh->m_child2;

			print_tp_neigh(tpNeigh, 3);

			if (fullOuterChild == tinfo.fbEnd2) {
				std::swap(tinfo.ebEnd1, tinfo.ebEnd2);
				std::swap(tinfo.fbEnd1, tinfo.fbEnd2);
			}
			OGDF_ASSERT(fullOuterChild == tinfo.fbEnd1);
			OGDF_ASSERT(tpNeigh->getParent() == tinfo.ebEnd1);
			OGDF_ASSERT(tinfo.tpPred == nullptr || tinfo.tpPred == tinfo.ebEnd2);

			print_tp_neigh(tpNeigh, 4);

			for (auto obs : m_observers) {
				obs->whenCNodeMerged(*this, tpNeigh, tpNeighSiblingsFlipped, fullNeigh,
						fullOuterChild);
			}

			auto parent = tpNeigh->getParent();

			tpNeigh->mergeIntoParent();

			print_tp_neigh(parent, 6);

			OGDF_ASSERT(central->areNeighborsAdjacent(fullNeigh, fullOuterChild));
			OGDF_ASSERT(central->areNeighborsAdjacent(otherNeigh, emptyOuterChild));

			if (tpNeigh->getLabelUnchecked() == NodeLabel::Partial) {
				OGDF_ASSERT(fullOuterChild->isFull());
				fullNeigh = tinfo.fbEnd2;
				otherEndOfFullBlock = tinfo.fbEnd1;
			} else {
				OGDF_ASSERT(tpNeigh->getLabelUnchecked() == NodeLabel::Empty);
				OGDF_ASSERT(fullOuterChild == tinfo.tpPred);
				otherEndOfFullBlock = nullptr;
			}
			OGDF_ASSERT(tinfo.tpPred == tinfo.ebEnd2 || tinfo.tpPred == nullptr);
			destroyNode(std::as_const(tpNeigh));
			PC_PROFILE_EXIT(2, "update_tp_cnode");
		}

		replaceTPNeigh(central, tpNeigh, nextTPNeigh, fullNeigh, otherEndOfFullBlock);
		if (nextTPNeigh != nullptr) {
			nextTPNeigh->tempInfo().replaceNeighbor(tpNeigh, central);
		}
		auto mergedNode = tpNeigh;
		tpNeigh = nextTPNeigh;
		count++;
		for (auto obs : m_observers) {
			obs->afterMerge(*this, tpNeigh, mergedNode);
		}
	}
	return count;
}

void PCTree::replaceTPNeigh(PCNode* central, PCNode* oldTPNeigh, PCNode* newTPNeigh,
		PCNode* newFullNeigh, PCNode* otherEndOfFullBlock) {
	PCNode::TempInfo& cinfo = central->tempInfo();
	OGDF_ASSERT(oldTPNeigh != nullptr);
	OGDF_ASSERT(cinfo.tpPred == cinfo.ebEnd1);
	OGDF_ASSERT(m_apexTPPred2 == nullptr || m_apexTPPred2 == cinfo.ebEnd2);
	if (oldTPNeigh == cinfo.tpPred) {
		cinfo.tpPred = cinfo.ebEnd1 = newTPNeigh;
		cinfo.fbEnd1 = newFullNeigh;
		if (cinfo.fbEnd2 == oldTPNeigh) {
			if (otherEndOfFullBlock != nullptr) {
				cinfo.fbEnd2 = otherEndOfFullBlock;
			} else {
				cinfo.fbEnd2 = newTPNeigh;
			}
		}
	} else {
		OGDF_ASSERT(oldTPNeigh == m_apexTPPred2);
		m_apexTPPred2 = cinfo.ebEnd2 = newTPNeigh;
		cinfo.fbEnd2 = newFullNeigh;
		if (cinfo.fbEnd1 == oldTPNeigh) {
			if (otherEndOfFullBlock != nullptr) {
				cinfo.fbEnd1 = otherEndOfFullBlock;
			} else {
				cinfo.fbEnd1 = newFullNeigh;
			}
		}
	}
	OGDF_ASSERT(m_apexTPPred2 != oldTPNeigh);
	OGDF_ASSERT(cinfo.ebEnd1 != oldTPNeigh);
	OGDF_ASSERT(cinfo.ebEnd2 != oldTPNeigh);
	OGDF_ASSERT(cinfo.fbEnd1 != oldTPNeigh);
	OGDF_ASSERT(cinfo.fbEnd2 != oldTPNeigh);
	if (cinfo.ebEnd1 != nullptr) {
		OGDF_ASSERT(central->areNeighborsAdjacent(cinfo.ebEnd1, cinfo.fbEnd1));
	}
	if (cinfo.ebEnd2 != nullptr) {
		OGDF_ASSERT(m_apexTPPred2 == nullptr || m_apexTPPred2 == cinfo.ebEnd2);
		OGDF_ASSERT(central->areNeighborsAdjacent(cinfo.ebEnd2, cinfo.fbEnd2));
	}
}

///////////////////////////////////////////////////////////////////////////////

// findTerminalPath utils

size_t PCTree::findEndOfFullBlock(PCNode* node, PCNode* pred, PCNode* curr, PCNode*& fullEnd,
		PCNode*& emptyEnd) const {
#ifdef OGDF_DEBUG
	PCNode* start = pred;
#endif
	fullEnd = pred;
	size_t count = 0;
	while (curr->getLabel() == NodeLabel::Full) {
		PCNode* next = node->getNextNeighbor(pred, curr);
		OGDF_ASSERT(next != start);
		fullEnd = pred = curr;
		curr = next;
		count++;
	}
	emptyEnd = curr;
	return count;
}

bool PCTree::setApexCandidate(PCNode* ac, bool fix) {
	if (ac->tempInfo().tpSucc == nullptr) {
		m_terminalPathLength++;
	} else if (ac->tempInfo().tpSucc != ac) {
		log << "  Note: already ascended from new " << (m_apexCandidateIsFix ? "" : "non-")
			<< "fix apex (" << ac << ") to (" << ac->tempInfo().tpSucc << ")" << std::endl;
	}
	ac->tempInfo().tpSucc = ac; // mark ac as processed
	if (m_apexCandidate == nullptr) {
		m_apexCandidate = ac;
		m_apexCandidateIsFix = fix;
		return true;
	} else if (m_apexCandidate == ac) {
		if (!m_apexCandidateIsFix && fix) {
			m_apexCandidateIsFix = true;
		}
		return true;
	} else {
		// we reached a node from which we can't ascend (nonFix), but also found the actual apex
		// (fix) check whether we ascended too far, i.e. the nonFix node is a parent of the fix apex
		if ((!fix && m_apexCandidateIsFix) || (fix && !m_apexCandidateIsFix)) {
			PCNode *fixAC, *nonFixAC;
			if (!fix && m_apexCandidateIsFix) {
				fixAC = m_apexCandidate;
				nonFixAC = ac;
			} else {
				OGDF_ASSERT(fix && !m_apexCandidateIsFix);
				fixAC = ac;
				nonFixAC = m_apexCandidate;
			}
			PCNode::TempInfo& fixTI = fixAC->tempInfo();
			PCNode::TempInfo& nonFixTI = nonFixAC->tempInfo();
			if (nonFixAC->getLabelUnchecked() == NodeLabel::Empty) {
				if (fixAC->getLabelUnchecked() == NodeLabel::Partial) {
					// if the fix apex is partial, it must be the tpPartialPred of the ascended-too-far empty/nonFix node
					if (nonFixTI.tpPartialPred == fixAC) {
						m_terminalPathLength -= nonFixTI.tpPartialHeight;
						m_apexCandidate = fixAC;
						m_apexCandidateIsFix = true;
						return true;
					}
				} else {
					OGDF_ASSERT(fixAC->getLabelUnchecked() == NodeLabel::Empty);
					// otherwise the fix apex is also empty, but they must have the same tpPartialPred
					if (nonFixTI.tpPartialPred == fixTI.tpPartialPred) {
						m_terminalPathLength -= (nonFixTI.tpPartialHeight - fixTI.tpPartialHeight);
						m_apexCandidate = fixAC;
						m_apexCandidateIsFix = true;
						return true;
					}
				}
			}
		}
		log << "Conflicting " << (m_apexCandidateIsFix ? "" : "non-") << "fix apex candidates ("
			<< m_apexCandidate << ") and (" << ac << ")" << std::endl;
		return false;
	}
}

///////////////////////////////////////////////////////////////////////////////

// updateTerminalPath utils

void PCTree::updateSingletonTerminalPath() {
	PCNode::TempInfo& atinfo = m_apexCandidate->tempInfo();
	OGDF_ASSERT(atinfo.tpPred == nullptr);
	OGDF_ASSERT(m_apexTPPred2 == nullptr);
	int fullNeighbors = atinfo.fullNeighbors.size();
	int emptyNeighbors = m_apexCandidate->getDegree() - fullNeighbors;
	if (m_apexCandidate->m_nodeType == PCNodeType::PNode && fullNeighbors > 1 && emptyNeighbors > 1) {
		// parent is handled and degree is always greater than 1
		PCNode* fullNode = splitOffFullPNode(m_apexCandidate, true);
		PCNode* parent = m_apexCandidate->getParent();
		if (parent != nullptr && parent->getLabel() == NodeLabel::Full) {
			m_apexCandidate->replaceWith(fullNode);
			for (auto obs : m_observers) {
				obs->nodeReplaced(*this, m_apexCandidate, fullNode);
			}
			fullNode->appendChild(m_apexCandidate);
			for (auto obs : m_observers) {
				obs->onApexMoved(*this, m_apexCandidate, fullNode, parent);
			}
		} else {
			m_apexCandidate->appendChild(fullNode);
		}
	}
}

PCNode* PCTree::splitOffFullPNode(PCNode* node, bool skip_parent) {
	auto& tinfo = node->tempInfo();
	auto* parent = node->getParent();
	if (tinfo.fullNeighbors.size() == 1) {
		PCNode* fullNode = tinfo.fullNeighbors.front();
		OGDF_ASSERT(fullNode != parent);
		OGDF_ASSERT(node->isParentOf(fullNode));
		fullNode->detach();
		for (auto obs : m_observers) {
			obs->fullNodeSplit(*this, fullNode);
		}
		return fullNode;
	}
	PCNode* fullNode = newNode(PCNodeType::PNode);
	fullNode->setLabel(NodeLabel::Full);
	for (PCNode* fullChild : tinfo.fullNeighbors) {
		if (skip_parent) {
			if (fullChild == parent) {
				continue;
			}
		} else {
			OGDF_ASSERT(fullChild != parent);
		}
		OGDF_ASSERT(node->isParentOf(fullChild));
		fullChild->detach();
		fullNode->appendChild(fullChild);
		fullNode->tempInfo().fullNeighbors.push_back(fullChild);
	}
	if (skip_parent && parent != nullptr && parent->getLabel() == NodeLabel::Full) {
		OGDF_ASSERT(fullNode->getDegree() >= 1);
	} else {
		OGDF_ASSERT(fullNode->getDegree() >= 2);
	}
	for (auto obs : m_observers) {
		obs->fullNodeSplit(*this, fullNode);
	}
	return fullNode;
}

}
