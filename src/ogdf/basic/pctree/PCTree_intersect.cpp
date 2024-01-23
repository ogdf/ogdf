#include <queue>

#include "PCTree.h"

bool pc_tree::PCTree::intersect(PCTree& other, PCTreeNodeArray<PCNode*>& mapping) {
	OGDF_HEAVY_ASSERT(checkValid() && other.checkValid());
	OGDF_ASSERT(leaves.size() == other.leaves.size());
	OGDF_ASSERT(mapping.registeredAt() == &other.forest->nodeArrayRegistry);
	if (other.isTrivial()) {
		return true;
	}

	size_t oldLeaves = leaves.size();
	PCTreeNodeArray<std::vector<PCNode*>> blockNodes(*this);
	PCTreeNodeArray<std::vector<PCNode*>> subtreeNodes(*this);
	PCTreeNodeArray<PCNode*> leafPartner(*this, nullptr);
	PCTreeNodeArray<bool> isFront(*this, false);
	bool possible = other.findNodeRestrictions(*this, mapping, blockNodes, subtreeNodes,
			leafPartner, isFront);
	restoreSubtrees(blockNodes, subtreeNodes, leafPartner, isFront);
	OGDF_ASSERT(oldLeaves == leaves.size());

	return possible;
}

bool pc_tree::PCTree::findNodeRestrictions(PCTree& applyTo, PCTreeNodeArray<PCNode*>& mapping,
		PCTreeNodeArray<std::vector<PCNode*>>& blockNodes,
		PCTreeNodeArray<std::vector<PCNode*>>& subtreeNodes, PCTreeNodeArray<PCNode*>& leafPartner,
		PCTreeNodeArray<bool>& isFront) {
	std::vector<PCNode*> fullNodeOrder;
	resetTempData();
	markFull(leaves.begin(), IntrusiveList<PCNode>::iterator(leaves.back()), &fullNodeOrder);
	changeRoot(leaves.back() == rootNode ? leaves.back()->child1 : leaves.back()->getParent());
	applyTo.changeRoot(mapping[leaves.back()] == applyTo.rootNode
					? mapping[leaves.back()]->child1
					: mapping[leaves.back()]->getParent());

	for (PCNode* node : fullNodeOrder) {
		auto nonLeafNeighborIt = std::find_if(node->neighbors().begin(), node->neighbors().end(),
				[](PCNode* n) { return !n->isLeaf(); });

		PCNode* startWith = nonLeafNeighborIt == node->neighbors().end() ? node->neighbors().first
																		 : *nonLeafNeighborIt;

		std::vector<PCNode*> consecutiveOriginal;
		std::vector<PCNode*> consecutiveOther;
		auto neighbors = node->neighbors(startWith);
		for (auto it = std::next(neighbors.begin()); it != neighbors.end(); ++it) {
			PCNode* current = *it;

			OGDF_ASSERT(current->isLeaf());
			consecutiveOther.push_back(current);
			consecutiveOriginal.push_back(mapping[current]);
			if (leafPartner[mapping[current]] != nullptr) {
				consecutiveOriginal.push_back(leafPartner[mapping[current]]);
			}

			if (node->getNodeType() == PCNodeType::CNode && consecutiveOther.size() >= 2) {
				std::vector<PCNode*> pair;
				PCNode* n1 = consecutiveOther[consecutiveOther.size() - 2];
				PCNode* n2 = consecutiveOther.back();
				pair.push_back(mapping[n1]);
				if (leafPartner[mapping[n1]] != nullptr) {
					pair.push_back(leafPartner[mapping[n1]]);
				}
				pair.push_back(mapping[n2]);
				if (leafPartner[mapping[n2]] != nullptr) {
					pair.push_back(leafPartner[mapping[n2]]);
				}
				if (!applyTo.makeConsecutive(pair)) {
					return false;
				}
			}
		}

		if (node != fullNodeOrder.back()) {
			PCNode* merged = mergeLeaves(consecutiveOther);
			if (!applyTo.makeConsecutive(consecutiveOriginal)) {
				return false;
			}
			std::vector<PCNode*> nodeOrder;
			applyTo.resetTempData();
			applyTo.markFull(consecutiveOriginal.begin(), consecutiveOriginal.end(), &nodeOrder);
			PCNode* partialNode = applyTo.firstPartial;
			OGDF_ASSERT(partialNode);
			OGDF_ASSERT(partialNode == applyTo.lastPartial);
			auto& fullNeighbors = partialNode->tempInfo().fullNeighbors;
			PCNode* ebEnd1;
			auto fbEnd1It = std::find_if(fullNeighbors.begin(), fullNeighbors.end(), [&](PCNode* n) {
				PCNode* sib1 = partialNode->getNextNeighbor(nullptr, n);
				if (!sib1->isFull()) {
					ebEnd1 = sib1;
					return true;
				}
				PCNode* sib2 = partialNode->getNextNeighbor(sib1, n);
				if (!sib2->isFull()) {
					ebEnd1 = sib2;
					return true;
				}
				return false;
			});
			OGDF_ASSERT(fbEnd1It != fullNeighbors.end());
			PCNode* fbEnd1 = *fbEnd1It;
			PCNode *fbEnd2, *ebEnd2;
			size_t count = applyTo.findEndOfFullBlock(partialNode, ebEnd1, fbEnd1, fbEnd2, ebEnd2);
			OGDF_ASSERT(partialNode->areNeighborsAdjacent(fbEnd1, ebEnd1));
			OGDF_ASSERT(partialNode->areNeighborsAdjacent(fbEnd2, ebEnd2));
			PCNode* newLeaf = applyTo.newNode(PCNodeType::Leaf);
			PCNode* newLeaf2 = nullptr;
			if (count >= 2) {
				/*
                 * If the block connecting the subtree to the partial node has length >= 2, it is important that the
                 * orientation of the block remains the same when reinserted in restoreSubtrees(). Because a single leaf
                 * could be rotated during makeConsecutive() calls, we insert two leaves instead to ensure the
                 * orientation remains consistent.
                 */
				OGDF_ASSERT(partialNode->getNodeType() == PCNodeType::CNode);
				newLeaf2 = applyTo.newNode(PCNodeType::Leaf);
				leafPartner[newLeaf] = newLeaf2;
				isFront[newLeaf] = true;
				leafPartner[newLeaf2] = newLeaf;
			}

			subtreeNodes[newLeaf].assign(nodeOrder.begin(), nodeOrder.end());
			subtreeNodes[newLeaf].insert(subtreeNodes[newLeaf].end(), consecutiveOriginal.begin(),
					consecutiveOriginal.end());

			PCNode* current = fbEnd1;
			while (current != ebEnd2) {
				OGDF_ASSERT(current->isFull());

				blockNodes[newLeaf].push_back(current);
				PCNode* next = partialNode->getNextNeighbor(ebEnd1, current);
				current->detach();
				current = next;
			}
			OGDF_ASSERT(blockNodes[newLeaf].size() == count);

			newLeaf->insertBetween(ebEnd1, ebEnd2);
			if (newLeaf2 != nullptr) {
				newLeaf2->insertBetween(newLeaf, ebEnd2);
				subtreeNodes[newLeaf2] = subtreeNodes[newLeaf];
				blockNodes[newLeaf2] = blockNodes[newLeaf];
			}

			for (PCNode* n : subtreeNodes[newLeaf]) {
				int oldId = n->nodeListIndex;
				applyTo.unregisterNode(n);
				if (n->getNodeType() == PCNodeType::CNode) {
					n->nodeListIndex = oldId;
				}
			}

			mapping[merged] = newLeaf;
			OGDF_HEAVY_ASSERT(applyTo.checkValid());
			OGDF_HEAVY_ASSERT(checkValid());
		}
	}
	return true;
}

void pc_tree::PCTree::restoreSubtrees(PCTreeNodeArray<std::vector<PCNode*>>& blockNodes,
		PCTreeNodeArray<std::vector<PCNode*>>& subtreeNodes, PCTreeNodeArray<PCNode*>& leafPartner,
		PCTreeNodeArray<bool>& isFront) {
	PCTreeNodeArray<bool> visited(*this, false);
	std::queue<PCNode*> queue;
	queue.push(leaves.front()->getParent());

	while (!queue.empty()) {
		PCNode* node = queue.front();
		queue.pop();

		PCNode* previous = node->neighbors().first;
		PCNode* current = node->getNextNeighbor(nullptr, previous);
		if (leafPartner[current] != nullptr && leafPartner[current] == previous) {
			previous = node->getNextNeighbor(previous, current);
		}
		int degree = node->neighbors().count();
		for (int i = 0; i < degree; i++) {
			OGDF_ASSERT(current->isValidNode(forest));

			if (!blockNodes[current].empty()) {
				// Replace the merged leaf (or two leaves) with the subtree it represents.
				visited[current] = true;
				PCNode* neighbor1 = previous;
				OGDF_ASSERT(node->areNeighborsAdjacent(current, previous));
				PCNode* neighbor2 = leafPartner[current] == nullptr
						? node->getNextNeighbor(previous, current)
						: node->getNextNeighbor(current, leafPartner[current]);
				OGDF_ASSERT(neighbor1 != neighbor2);
				if (!isFront[current]) {
					std::reverse(blockNodes[current].begin(), blockNodes[current].end());
				}

				for (PCNode* n : subtreeNodes[current]) {
					if (n->getNodeType() == PCNodeType::CNode) {
						cNodeCount++;
						OGDF_ASSERT(getForest()->cNodes.at(n->nodeListIndex) == nullptr);
						getForest()->cNodes[n->nodeListIndex] = n;
					} else {
						registerNode(n);
					}
				}

				current->detach();
				if (leafPartner[current] != nullptr) {
					leafPartner[current]->detach();
				}

				for (PCNode* n : blockNodes[current]) {
					n->insertBetween(neighbor1, neighbor2);
					neighbor1 = n;
				}

				degree += blockNodes[current].size();
				PCNode* tmp = current;
				current = blockNodes[current].front();
				blockNodes[tmp].clear();
				if (leafPartner[tmp] != nullptr) {
					destroyNode(leafPartner[tmp]);
				}
				destroyNode(tmp);
			} else {
				if (!current->isLeaf() && !visited[current]) {
					queue.push(current);
				}
				visited[current] = true;
				PCNode* tmp = current;
				current = node->getNextNeighbor(previous, current);
				previous = tmp;
				OGDF_ASSERT(previous != current);
			}
		}
	}

	OGDF_ASSERT(checkValid());
}
