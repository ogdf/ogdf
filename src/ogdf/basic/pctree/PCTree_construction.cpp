#include <regex>

#include "PCNode.h"
#include "PCTree.h"

using namespace pc_tree;

PCTree::PCTree(int leafNum, std::vector<PCNode*>* added) : PCTree() {
	OGDF_ASSERT(leafNum > 2);
	rootNode = newNode(PCNodeType::PNode);
	insertLeaves(leafNum, rootNode, added);
}

PCTree::PCTree(const std::string& str, bool keep_ids) : PCTree() {
	std::string s = std::regex_replace(str, std::regex("\\s+"), ""); //remove whitespaces

	std::stringstream ss(s);
	std::stack<PCNode*> stack;
	int nextIndex = 0;
	bool indexUsed = true;
	char previousChar = ' ';

	while (!ss.eof()) {
		char nextChar = ss.peek();

		if (isdigit(nextChar) || nextChar == '-') {
			ss >> nextIndex;
			if (keep_ids) {
				forest->nextNodeId = std::max(nextIndex + 1, forest->nextNodeId);
			} else {
				nextIndex = forest->nextNodeId++;
			}
			indexUsed = false;
		} else {
			ss.ignore();

			PCNode* parent = stack.empty() ? nullptr : stack.top();
			PCNode* created = nullptr;

			switch (nextChar) {
			case ',':
				if (previousChar != ']' && previousChar != ')') {
					if (stack.empty()) {
						throw std::invalid_argument("Invalid PC-Tree");
					}

					created = newNode(PCNodeType::Leaf, parent, nextIndex);
				}
				break;
			case '{':
				if (stack.empty() && getLeafCount() > 0) {
					throw std::invalid_argument("Invalid PC-Tree");
				}

				OGDF_ASSERT(parent == nullptr);
				created = newNode(PCNodeType::Leaf, parent, nextIndex);
				stack.push(created);
				break;
			case '[':
				if (stack.empty() && getLeafCount() > 0) {
					throw std::invalid_argument("Invalid PC-Tree");
				}

				created = newNode(PCNodeType::CNode, parent, nextIndex);
				stack.push(created);
				break;
			case '(':
				if (stack.empty() && getLeafCount() > 0) {
					throw std::invalid_argument("Invalid PC-Tree");
				}

				created = newNode(PCNodeType::PNode, parent, nextIndex);
				stack.push(created);
				break;
			case ']':
				if (stack.empty() || stack.top()->nodeType != PCNodeType::CNode) {
					throw std::invalid_argument("Invalid PC-Tree");
				}

				if (previousChar != ']' && previousChar != ')') {
					created = newNode(PCNodeType::Leaf, parent, nextIndex);
				}
				stack.pop();
				break;
			case ')':
				if (stack.empty() || stack.top()->nodeType != PCNodeType::PNode) {
					throw std::invalid_argument("Invalid PC-Tree");
				}

				if (previousChar != ']' && previousChar != ')') {
					created = newNode(PCNodeType::Leaf, parent, nextIndex);
				}
				stack.pop();
				break;
			case '}':
				if (stack.empty() || stack.top()->nodeType != PCNodeType::Leaf) {
					throw std::invalid_argument("Invalid PC-Tree");
				}

				if (previousChar != ']' && previousChar != ')') {
					created = newNode(PCNodeType::Leaf, parent, nextIndex);
				}
				stack.pop();
				OGDF_ASSERT(stack.empty());
				break;
			default:
				break;
			}

			if (created) {
				if (indexUsed) {
					throw std::invalid_argument("Invalid PC-Tree");
				}
				indexUsed = true;
			}

			previousChar = nextChar;
		}
	}
	if (!stack.empty()) {
		throw std::invalid_argument("Invalid PC-Tree");
	}
}

PCTree::PCTree::PCTree(const PCTree& other, PCTreeNodeArray<PCNode*>& nodeMapping, bool keep_ids)
	: PCTree() {
	nodeMapping.init(other);
	for (PCNode* other_node : other.allNodes()) {
		PCNode* parent = other_node->getParent();
		int id = -1;
		if (keep_ids) {
			id = other_node->id;
			forest->nextNodeId = std::max(id + 1, forest->nextNodeId);
		}
		OGDF_ASSERT((parent == nullptr) == (other.rootNode == other_node));
		if (parent == nullptr) {
			rootNode = nodeMapping[other_node] = newNode(other_node->getNodeType(), nullptr, id);
		} else {
			nodeMapping[other_node] = newNode(other_node->getNodeType(), nodeMapping[parent], id);
		}
	}
	OGDF_ASSERT(nodeMapping[other.rootNode] == rootNode);
	OGDF_ASSERT(other.getLeafCount() == getLeafCount());
	OGDF_ASSERT(other.getPNodeCount() == getPNodeCount());
	OGDF_ASSERT(other.getCNodeCount() == getCNodeCount());
}

PCTree::~PCTree() {
	OGDF_ASSERT(checkValid());
	// get rid of any degree <= 2 root, including a root leaf and possible degree 2 descendants
	while (rootNode != nullptr && rootNode->childCount <= 2) {
		PCNode* old_root = rootNode;
		if (rootNode->childCount == 2) {
			rootNode = old_root->getChild1();
			rootNode->detach();
			PCNode* child = old_root->getChild2();
			child->detach();
			rootNode->appendChild(child);
		} else if (rootNode->childCount == 1) {
			rootNode = old_root->getOnlyChild();
			rootNode->detach();
		} else {
			rootNode = nullptr;
		}
		destroyNode(old_root);
	}
	OGDF_ASSERT(checkValid());

	while (!leaves.empty()) {
		PCNode* node = leaves.back();
		PCNode* parent = node->getParent();
		PCNodeType type = node->getNodeType();
		bool is_root = node == rootNode;
		OGDF_ASSERT((parent == nullptr) == is_root);
		if (is_root) {
			OGDF_ASSERT(leaves.size() == 1);
			rootNode = nullptr;
		}
		node->detach();
		destroyNode(node);
		if (type != PCNodeType::Leaf) {
			leaves.pop_back(); // destroyNode(leaf) automatically removes it from the list
		}
		if (parent != nullptr && parent->childCount == 0) {
			leaves.push_back(parent);
		}
		if (is_root) {
			OGDF_ASSERT(leaves.empty());
			OGDF_ASSERT(rootNode == nullptr);
		} else {
			OGDF_ASSERT(!leaves.empty());
			OGDF_ASSERT(rootNode != nullptr);
		}
	}

	if (!externalForest) {
		delete forest;
	}

	OGDF_ASSERT(rootNode == nullptr);
	OGDF_ASSERT(pNodeCount == 0);
	OGDF_ASSERT(cNodeCount == 0);
}

void PCTree::registerNode(PCNode* node) {
	if (node->nodeType == PCNodeType::Leaf) {
		leaves.push_back(node);
	} else if (node->nodeType == PCNodeType::PNode) {
		pNodeCount++;
	} else {
		OGDF_ASSERT(node->nodeType == PCNodeType::CNode);
		node->nodeListIndex = forest->parents.makeSet();
		OGDF_ASSERT(forest->cNodes.size() == node->nodeListIndex);
		forest->cNodes.push_back(node);
		cNodeCount++;
	}
}

void PCTree::unregisterNode(PCNode* node) {
	if (node->nodeType == PCNodeType::Leaf) {
		leaves.erase(node);
	} else if (node->nodeType == PCNodeType::PNode) {
		pNodeCount--;
	} else {
		OGDF_ASSERT(node->nodeType == PCNodeType::CNode);
		OGDF_ASSERT(forest->cNodes.at(node->nodeListIndex) == node);
		forest->cNodes[node->nodeListIndex] = nullptr;
		cNodeCount--;
		node->nodeListIndex = UNIONFINDINDEX_EMPTY;
	}
}

PCNode* PCTree::newNode(PCNodeType type, PCNode* parent, int id) {
	int oldTableSize = forest->nodeArrayRegistry.keyArrayTableSize();
	PCNode* node;
#ifdef PCTREE_REUSE_NODES
	if (forest->reusableNodes) {
		node = forest->reusableNodes;
		forest->reusableNodes = forest->reusableNodes->parentPNode;
		node->parentPNode = nullptr;
		node->timestamp = 0;
		if (id >= 0) {
			node->id = id;
			forest->nextNodeId = std::max(forest->nextNodeId, id + 1);
		} // else we can't re-use the old ID after clear was called
		else {
			node->id = forest->nextNodeId++;
		}
		node->changeType(type);
	} else
#endif
	{
		if (id < 0) {
			id = forest->nextNodeId++;
		} else {
			forest->nextNodeId = std::max(forest->nextNodeId, id + 1);
		}
		node = new PCNode(forest, id, type);
	}
	registerNode(node);
	if (parent != nullptr) {
		parent->appendChild(node);
	} else if (rootNode == nullptr) {
		rootNode = node;
	}
	if (oldTableSize != forest->nodeArrayRegistry.keyArrayTableSize()) {
		forest->nodeArrayRegistry.enlargeArrayTables();
	}

	for (auto obs : observers) {
		obs->onNodeCreate(node);
	}

	return node;
}

void PCTree::destroyNode(PCNode* const& node) {
	OGDF_ASSERT(node->forest == forest);
	OGDF_ASSERT(node->isDetached());
	OGDF_ASSERT(node->childCount == 0);
	OGDF_ASSERT(node->child1 == nullptr);
	OGDF_ASSERT(node->child2 == nullptr);
	OGDF_ASSERT(node != rootNode);
	unregisterNode(node);
#ifdef PCTREE_REUSE_NODES
	node->parentPNode = forest->reusableNodes;
	forest->reusableNodes = node;
#else
	delete node;
#endif
}

PCNodeType PCTree::changeNodeType(PCNode* node, PCNodeType newType) {
	PCNodeType oldType = node->nodeType;
	if (oldType == newType) {
		return oldType;
	}

	UnionFindIndex oldIndex = node->nodeListIndex;
	unregisterNode(node);
	node->changeType(newType);
	registerNode(node);

	if (oldType == PCNodeType::CNode || newType == PCNodeType::CNode) {
		PCNode* pred = nullptr;
		PCNode* curr = node->child1;
		size_t children = 0;
		while (curr != nullptr) {
			if (oldType == PCNodeType::CNode) {
				OGDF_ASSERT(curr->parentPNode == nullptr);
				OGDF_ASSERT(forest->parents.find(curr->parentCNodeId) == oldIndex);
			} else {
				OGDF_ASSERT(curr->parentPNode == node);
				OGDF_ASSERT(curr->parentCNodeId == UNIONFINDINDEX_EMPTY);
			}
			if (newType == PCNodeType::CNode) {
				curr->parentPNode = nullptr;
				curr->parentCNodeId = node->nodeListIndex;
			} else {
				curr->parentPNode = node;
				curr->parentCNodeId = UNIONFINDINDEX_EMPTY;
			}
			children++;
			proceedToNextSibling(pred, curr);
		}
		OGDF_ASSERT(children == node->childCount);
		OGDF_ASSERT(pred == node->child2);
	}

	return oldType;
}

void PCTree::insertLeaves(int count, PCNode* parent, std::vector<PCNode*>* added) {
	OGDF_ASSERT(parent != nullptr);
	OGDF_ASSERT(parent->forest == forest);
	if (added) {
		added->reserve(added->size() + count);
	}
	for (int i = 0; i < count; i++) {
		PCNode* leaf = newNode(PCNodeType::Leaf);
		parent->appendChild(leaf);
		if (added) {
			added->push_back(leaf);
		}
	}
}

void PCTree::replaceLeaf(int leafCount, PCNode* leaf, std::vector<PCNode*>* added) {
	OGDF_ASSERT(leaf && leaf->forest == forest);
	OGDF_ASSERT(leaf->isLeaf());
	OGDF_ASSERT(leafCount > 1);
	if (getLeafCount() <= 2) {
		changeNodeType(leaf->getParent(), PCNodeType::PNode);
		insertLeaves(leafCount, leaf->getParent(), added);
		leaf->detach();
		destroyNode(leaf);
	} else {
		changeNodeType(leaf, PCNodeType::PNode);
		insertLeaves(leafCount, leaf, added);
	}
}

void PCTree::destroyLeaf(PCNode* leaf) {
	OGDF_ASSERT(leaf->getNodeType() == PCNodeType::Leaf);
	OGDF_ASSERT(leaf != rootNode);

	PCNode* parent = leaf->getParent();
	leaf->detach();
	destroyNode(leaf);

	/* assume the PC-tree is valid, so a childCount of 0 is impossible, since every inner node must
     * have at least 2 children (except for child of root node) */

	if (parent->childCount == 1) {
		if (rootNode->getNodeType() == PCNodeType::Leaf) {
			if (parent->getChild1()->getNodeType() != PCNodeType::Leaf
					|| rootNode->getChild1() != parent) {
				PCNode* child = parent->getChild1();
				child->detach();
				parent->replaceWith(child);
				destroyNode(parent);
			}
		} else {
			PCNode* child = parent->getChild1();
			if (parent != rootNode) {
				child->detach();
				parent->replaceWith(child);
				destroyNode(parent);
			} else if (child->getNodeType() != PCNodeType::Leaf) {
				PCNode* root = rootNode;
				root->detach();
				root->childCount = 0;
				root->child1 = root->child2 = nullptr;
				child->parentCNodeId = UNIONFINDINDEX_EMPTY;
				child->parentPNode = nullptr;
				setRoot(child);
				destroyNode(root);
			}
		}
	}
}

// identified `inserted`s root leaf with `base`s base_leaf
void PCTree::insertTree(PCNode* at, PCTree* inserted) {
	//    OGDF_ASSERT(a != nullptr && b != nullptr);
	//    OGDF_ASSERT(a != b);
	//    OGDF_ASSERT(a->forest == b->forest);
	//    OGDF_ASSERT(a->externalForest && b->externalForest);
	//    OGDF_ASSERT(a->getRootNode()->getNodeType() == PCNodeType::Leaf
	//                && b->getRootNode()->getNodeType() == PCNodeType::Leaf);
	OGDF_ASSERT(checkValid());
	OGDF_ASSERT(inserted->checkValid());

	observers.splice(observers.end(), inserted->observers);
	leaves.splice(leaves.end(), inserted->leaves);
	pNodeCount += inserted->pNodeCount;
	cNodeCount += inserted->cNodeCount;

	PCNode* root = inserted->rootNode;
	inserted->rootNode = nullptr;
	inserted->pNodeCount = inserted->cNodeCount = 0;
	delete inserted;
	inserted = nullptr;

	while (root->getChildCount() == 1) {
		PCNode* new_root = root->getOnlyChild();
		new_root->detach();
		destroyNode(root);
		root = new_root;
	}

	if (at->isLeaf()) {
		OGDF_ASSERT(!at->isDetached());
		PCNode* parent = at->getParent();
		if (!root->isLeaf() && parent->getDegree() == 2) {
			// remove degree 2 node in a tree with two leaves if it got further leaves
			at->detach();
			if (parent->isDetached()) {
				// use the inserted root as root
				rootNode = root;
				PCNode* child = parent->getOnlyChild();
				child->detach();
				rootNode->appendChild(child);
			} else {
				// we can use the other leaf as root
				parent->replaceWith(root);
			}
			destroyNode(parent);
		} else {
			at->replaceWith(root);
		}
		destroyNode(at);
	} else {
		at->appendChild(root);
	}
	OGDF_ASSERT(checkValid());
}
