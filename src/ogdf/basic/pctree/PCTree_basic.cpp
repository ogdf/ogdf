#include <ogdf/basic/pctree/PCNode.h>
#include <ogdf/basic/pctree/PCTree.h>

#include <queue>
#include <stack>
#include <variant>

using namespace pc_tree;
using namespace ogdf;

bool PCTree::isTrivial() const {
	if (leaves.empty()) {
		OGDF_ASSERT(rootNode == nullptr);
		return true;
	}
	return rootNode->getNodeType() == PCNodeType::PNode && rootNode->childCount == leaves.size();
}

void PCTree::getTree(Graph& tree, GraphAttributes* g_a, PCTreeNodeArray<ogdf::node>& pc_repr,
		ogdf::NodeArray<PCNode*>* g_repr, bool mark_full, bool show_sibs) const {
	tree.clear();

	if (leaves.empty()) {
		return;
	}

	bool nodeGraphics = false, nodeLabel = false;
	if (g_a != nullptr) {
		nodeGraphics = g_a->has(GraphAttributes::nodeGraphics);
		nodeLabel = g_a->has(GraphAttributes::nodeLabel);
	}

	for (PCNode* pc_node : allNodes()) {
		ogdf::node g_node = pc_repr[pc_node] = tree.newNode();
		if (g_repr != nullptr) {
			(*g_repr)[g_node] = pc_node;
		}

		if (nodeGraphics) {
			if (pc_node->nodeType == PCNodeType::CNode) {
				g_a->shape(g_node) = Shape::Rhomb;
			} else if (pc_node->nodeType == PCNodeType::PNode) {
				g_a->shape(g_node) = Shape::Ellipse;
			} else {
				OGDF_ASSERT(pc_node->isLeaf());
				g_a->shape(g_node) = Shape::Triangle;
			}
			if (mark_full) {
				NodeLabel label = pc_node->getLabel();
				if (label == NodeLabel::Full) {
					g_a->fillColor(g_node) = (Color(Color::Name::Darkblue));
				} else if (label == NodeLabel::Partial) {
					g_a->fillColor(g_node) = (Color(Color::Name::Lightblue));
				}
			}
		}
		if (nodeLabel) {
			g_a->label(g_node) = std::to_string(pc_node->index());
		}
		PCNode* parent = pc_node->getParent();
		if (parent != nullptr) {
			ogdf::edge e = tree.newEdge(pc_repr[parent], g_node);
			if (pc_node->parentPNode != nullptr) {
				g_a->strokeType(e) = ogdf::StrokeType::Solid;
			} else if (pc_node->parentCNodeId != UNIONFINDINDEX_EMPTY) {
				g_a->strokeType(e) = ogdf::StrokeType::Dot;
			} else {
				g_a->strokeType(e) = ogdf::StrokeType::None;
			}

			if (parent->getChild1() == pc_node) {
				g_a->strokeColor(e) = (Color(Color::Name::Red));
				g_a->label(e) = "c1";
			}
			if (parent->getChild2() == pc_node) {
				g_a->strokeColor(e) = (Color(Color::Name::Green));
				g_a->label(e) = "c2";
			}
		}
	}
	if (!show_sibs) {
		return;
	}
	for (PCNode* pc_node : allNodes()) {
		ogdf::node g_node = pc_repr[pc_node];
		//        PCNode* parent = pc_node->getParent();
		if (pc_node->getSibling1() != nullptr) {
			ogdf::edge e = tree.newEdge(g_node, pc_repr[pc_node->getSibling1()]);
			g_a->strokeType(e) = ogdf::StrokeType::Dash;
			g_a->strokeColor(e) = Color(Color::Name::Red);
			g_a->label(e) = "s1";
		}
		if (pc_node->getSibling2() != nullptr) {
			ogdf::edge e = tree.newEdge(g_node, pc_repr[pc_node->getSibling2()]);
			g_a->strokeType(e) = ogdf::StrokeType::Dash;
			g_a->strokeColor(e) = Color(Color::Name::Green);
			g_a->label(e) = "s2";
		}
	}
}

std::ostream& operator<<(std::ostream& os, const PCTree& tree) { return os << &tree; }

std::ostream& operator<<(std::ostream& os, const PCTree* tree) {
	std::stack<std::variant<PCNode*, std::string>> stack;
	stack.push(tree->rootNode);
	if (tree->rootNode == nullptr) {
		return os << "empty";
	}

	while (!stack.empty()) {
		auto next = stack.top();
		stack.pop();

		// Next stack element is either a string we need to append or an arc pointing to a
		// subtree we still need to process.
		if (std::holds_alternative<std::string>(next)) {
			os << std::get<std::string>(next);
			continue;
		}

		PCNode* base = std::get<PCNode*>(next);
		if (base->nodeType == PCNodeType::CNode) {
			os << base->id << ":[";
			stack.push("]");
		} else if (base->nodeType == PCNodeType::PNode) {
			os << base->id << ":(";
			stack.push(")");
		} else {
			OGDF_ASSERT(base->nodeType == PCNodeType::Leaf);
			if (base != tree->rootNode) {
				os << base->id;
				continue;
			} else {
				os << base->id << ":{";
				stack.push("}");
			}
		}

		bool space = false;
		for (PCNode* node : base->children()) {
			if (space) {
				stack.push(", ");
			}
			OGDF_ASSERT(node != nullptr);
			stack.emplace(node);
			space = true;
		}
	}

	return os;
}

void pc_tree::uid_utils::nodeToID(std::ostream& os, PCNode* n, int pos) {
	os << n->index();
	if (!n->isLeaf()) {
		os << ":";
	}
}

void pc_tree::uid_utils::nodeToPosition(std::ostream& os, PCNode* n, int pos) {
	os << pos;
	if (!n->isLeaf()) {
		os << ":";
	}
}

void pc_tree::uid_utils::leafToID(std::ostream& os, PCNode* n, int pos) {
	if (n->isLeaf()) {
		os << n->index();
	}
}

void pc_tree::uid_utils::leafToPosition(std::ostream& os, PCNode* n, int pos) {
	if (n->isLeaf()) {
		os << pos;
	}
}

bool pc_tree::uid_utils::compareNodesByID(PCNode* a, PCNode* b) { return a->index() < b->index(); }

std::ostream& PCTree::uniqueID(std::ostream& os,
		const std::function<void(std::ostream& os, PCNode*, int)>& printNode,
		const std::function<bool(PCNode*, PCNode*)>& compareNodes) {
	if (rootNode == nullptr) {
		return os << "empty";
	}
	std::vector<PCNode*> sortedLeaves(leaves.begin(), leaves.end());
	std::sort(sortedLeaves.begin(), sortedLeaves.end(), compareNodes);

	PCTreeNodeArray<int> order(*this, -1);
	int i = 0;
	for (PCNode* leaf : sortedLeaves) {
		order[leaf] = i++;
	}

	PCNode* lastLeaf = sortedLeaves.back();
	sortedLeaves.resize(leaves.size() - 1);
	std::vector<PCNode*> fullNodeOrder;
	{
		resetTempData();
		markFull(sortedLeaves.begin(), sortedLeaves.end(), &fullNodeOrder);
	}
	for (PCNode* node : fullNodeOrder) {
		order[node] = i++;
	}

	std::stack<std::variant<PCNode*, std::string>> stack;
	stack.push(fullNodeOrder.back());
	while (!stack.empty()) {
		auto next = stack.top();
		stack.pop();

		if (std::holds_alternative<std::string>(next)) {
			os << std::get<std::string>(next);
			continue;
		}

		PCNode* node = std::get<PCNode*>(next);
		std::list<PCNode*> children;
		if (node->nodeType == PCNodeType::CNode) {
			printNode(os, node, order[node]);
			os << "[";
			stack.push("]");
			if (node == fullNodeOrder.back()) {
				children.assign(node->neighbors().begin(), node->neighbors().end());
				PCNode* minChild = *std::min_element(children.begin(), children.end(),
						[&order](PCNode* elem, PCNode* min) { return order[elem] < order[min]; });
				while (children.front() != minChild) {
					children.push_back(children.front());
					children.pop_front();
				}
				PCNode* second = *(++children.begin());
				if (order[second] > order[children.back()]) {
					children.push_back(children.front());
					children.pop_front();
					children.reverse();
				}
				second = *(++children.begin());
				OGDF_ASSERT(children.front() == minChild);
				OGDF_ASSERT(order[second] < order[children.back()]);
			} else {
				PCNode* informedNeighbor = nullptr;
				for (PCNode* neigh : node->neighbors()) {
					if (order[neigh] > order[node]) {
						OGDF_ASSERT(informedNeighbor == nullptr);
						informedNeighbor = neigh;
					}
				}
				OGDF_ASSERT(informedNeighbor != nullptr);
				PCNode* neigh1 = node->getNextNeighbor(nullptr, informedNeighbor);
				PCNode* neigh2 = node->getNextNeighbor(neigh1, informedNeighbor);
				if (order[neigh2] < order[neigh1]) {
					std::swap(neigh1, neigh2);
				}
				for (PCNode *pred = informedNeighbor, *curr = neigh1; curr != informedNeighbor;
						node->proceedToNextNeighbor(pred, curr)) {
					children.push_back(curr);
				}
			}
		} else if (node->nodeType == PCNodeType::PNode) {
			printNode(os, node, order[node]);
			if (node->getDegree() <= 3) {
				os << "[";
				stack.push("]");
			} else {
				os << "(";
				stack.push(")");
			}
			std::vector<PCNode*>& fullNeighbors = node->tempInfo().fullNeighbors;
			children.assign(fullNeighbors.begin(), fullNeighbors.end());
			if (node == fullNodeOrder.back()) {
				children.push_back(lastLeaf);
			}
			children.sort([&order](PCNode* a, PCNode* b) { return order[a] < order[b]; });
		} else {
			OGDF_ASSERT(node->nodeType == PCNodeType::Leaf);
			printNode(os, node, order[node]);
			continue;
		}
		if (node == fullNodeOrder.back()) {
			OGDF_ASSERT(children.size() == node->tempInfo().fullNeighbors.size() + 1);
		} else {
			OGDF_ASSERT(children.size() == node->tempInfo().fullNeighbors.size());
		}
		OGDF_ASSERT(order[children.front()] < order[children.back()]);

		bool space = false;
		for (PCNode* child : children) {
			if (space) {
				stack.push(", ");
			}
			OGDF_ASSERT(child != nullptr);
			stack.push(child);
			space = true;
		}
	}

	return os;
}

std::ostream& operator<<(std::ostream& os, const pc_tree::PCNodeType t) {
	switch (t) {
	case pc_tree::PCNodeType::Leaf:
		return os << "Leaf";
	case pc_tree::PCNodeType::PNode:
		return os << "PNode";
	case pc_tree::PCNodeType::CNode:
		return os << "CNode";
	default:
		OGDF_ASSERT(false);
		return os << "PCNodeType???";
	}
}

std::ostream& operator<<(std::ostream& os, const pc_tree::NodeLabel l) {
	switch (l) {
	case pc_tree::NodeLabel::Unknown:
		return os << "Empty/Unknown";
	case pc_tree::NodeLabel::Partial:
		return os << "Partial";
	case pc_tree::NodeLabel::Full:
		return os << "Full";
	default:
		OGDF_ASSERT(false);
		return os << "NodeLabel???";
	}
}

template class pc_tree::PCTreeRegistry<PCNode*>;

bool PCTree::isValidOrder(const std::vector<PCNode*>& order) const {
	OGDF_ASSERT(order.size() == leaves.size());
	PCTreeNodeArray<PCNode*> leafMapping(*this);
	PCTree copy(*this, leafMapping);
	PCNode* previous = nullptr;
	for (PCNode* node : order) {
		OGDF_ASSERT(node->forest == forest);
		if (previous == nullptr) {
			previous = node;
			continue;
		}
		if (!copy.makeConsecutive({leafMapping[previous], leafMapping[node]})) {
			return false;
		}
		previous = node;
	}
#ifdef OGDF_DEBUG
	OGDF_ASSERT(copy.possibleOrders<int>() == 2);
	OGDF_ASSERT(copy.makeConsecutive({leafMapping[order.front()], leafMapping[order.back()]}));
	std::list<PCNode*> res_order;
	copy.currentLeafOrder(res_order);
	OGDF_ASSERT(res_order.size() == order.size());
	for (size_t i = 0; res_order.front() != leafMapping[order.front()]; i++) {
		res_order.push_back(res_order.front());
		res_order.pop_front();
		OGDF_ASSERT(i < order.size());
	}
	if (res_order.back() != leafMapping[order.back()]) {
		std::reverse(res_order.begin(), res_order.end());
		res_order.push_front(res_order.back());
		res_order.pop_back();
	}
	OGDF_ASSERT(res_order.size() == order.size());
	int i = 0;
	for (PCNode* n : res_order) {
		OGDF_ASSERT(n == leafMapping[order.at(i)]);
		++i;
	}
#endif

	return true;
}

bool PCTree::checkValid(bool allow_small_deg) const {
	if (rootNode == nullptr) {
		OGDF_ASSERT(leaves.size() == 0);
		OGDF_ASSERT(pNodeCount == 0);
		OGDF_ASSERT(cNodeCount == 0);
		return true;
	}

	OGDF_ASSERT(rootNode->getParent() == nullptr);
	for (PCNode* leaf : leaves) {
		OGDF_ASSERT(leaf->isLeaf());
		if (leaf == rootNode) {
			OGDF_ASSERT(leaf->getChildCount() <= 1);
		} else {
			OGDF_ASSERT(leaf->getChildCount() == 0);
		}
	}

	if (leaves.size() == 0) {
		OGDF_ASSERT(rootNode->getDegree() == 0);
		OGDF_ASSERT(!rootNode->isLeaf());
		OGDF_ASSERT(pNodeCount + cNodeCount == 1);
		return true;
	} else if (leaves.size() == 1) {
		if (rootNode->isLeaf()) {
			OGDF_ASSERT(rootNode == leaves.front());
			if (rootNode->getDegree() == 1) {
				OGDF_ASSERT(pNodeCount + cNodeCount == 1);
				OGDF_ASSERT(rootNode->getOnlyChild()->isValidNode(forest));
				OGDF_ASSERT(rootNode->getOnlyChild()->getChildCount() == 0);
				OGDF_ASSERT(!rootNode->getOnlyChild()->isLeaf());
			} else {
				OGDF_ASSERT(pNodeCount + cNodeCount == 0);
			}
		} else {
			OGDF_ASSERT(rootNode->getOnlyChild()->isValidNode(forest));
			OGDF_ASSERT(rootNode->getOnlyChild() == leaves.front());
			OGDF_ASSERT(rootNode->getOnlyChild()->isLeaf());
		}
		return true;
	} else if (leaves.size() == 2) {
		OGDF_ASSERT(pNodeCount + cNodeCount == 1);
		if (rootNode->isLeaf()) {
			OGDF_ASSERT(!rootNode->getOnlyChild()->isLeaf());
			if (rootNode == leaves.front()) {
				OGDF_ASSERT(rootNode->getOnlyChild()->getOnlyChild() == leaves.back());
			} else {
				OGDF_ASSERT(rootNode == leaves.back());
				OGDF_ASSERT(rootNode->getOnlyChild()->getOnlyChild() == leaves.front());
			}
		} else {
			if (rootNode->getChild1() == leaves.front()) {
				OGDF_ASSERT(rootNode->getChild2() == leaves.back());
			} else {
				OGDF_ASSERT(rootNode->getChild2() == leaves.front());
				OGDF_ASSERT(rootNode->getChild1() == leaves.back());
			}
		}
		return true;
	}
	OGDF_ASSERT(leaves.size() > 2);

	// bottom-up
	std::queue<PCNode*> todo;
	for (PCNode* leaf : leaves) {
		if (leaf != rootNode) {
			todo.push(leaf);
		}
	}
	if (rootNode->isLeaf()) {
		OGDF_ASSERT(todo.size() == leaves.size() - 1);
	} else {
		OGDF_ASSERT(todo.size() == leaves.size());
	}
	PCNode* lastLeaf = todo.back();
	bool leaves_done = false, root_found = false;
	size_t leaves_found = 0, p_nodes_found = 0, c_nodes_found = 0;
	std::vector<PCNode*> id_seen(forest->nextNodeId, nullptr);
	while (!todo.empty()) {
		PCNode* node = todo.front();
		todo.pop();

		if (id_seen.at(node->id) == node) {
			OGDF_ASSERT(leaves_done);
			continue;
		}
		OGDF_ASSERT(node->isValidNode(forest));
		OGDF_ASSERT(id_seen[node->id] == nullptr);
		id_seen[node->id] = node;
		PCNode* parent = node->getParent();
		OGDF_ASSERT((node == rootNode) == (parent == nullptr));
		if (node == rootNode) {
			OGDF_ASSERT(leaves_done);
			root_found = true;
		} else {
			OGDF_ASSERT(leaves_done == (node->nodeType != PCNodeType::Leaf));
			todo.push(parent);
		}

		if (node->nodeType == PCNodeType::PNode) {
			p_nodes_found++;
		} else if (node->nodeType == PCNodeType::CNode) {
			c_nodes_found++;
		} else {
			OGDF_ASSERT(node->isLeaf());
			leaves_found++;
			if (node == rootNode) {
				OGDF_ASSERT(node->getChildCount() == 1);
			} else {
				OGDF_ASSERT(node->getChildCount() == 0);
			}
		}

		if (node == lastLeaf) {
			leaves_done = true;
		}
	}
	OGDF_ASSERT(leaves_done);
	OGDF_ASSERT(root_found);
	OGDF_ASSERT(leaves_found == leaves.size());
	OGDF_ASSERT(p_nodes_found == pNodeCount);
	OGDF_ASSERT(c_nodes_found == cNodeCount);

	// top-down
	leaves_found = p_nodes_found = c_nodes_found = 0;
	OGDF_ASSERT(todo.empty());
	todo.push(rootNode);
	while (!todo.empty()) {
		PCNode* node = todo.front();
		todo.pop();
		OGDF_ASSERT(id_seen[node->id] == node);

		if (node->nodeType == PCNodeType::PNode) {
			p_nodes_found++;
		} else if (node->nodeType == PCNodeType::CNode) {
			c_nodes_found++;
		} else {
			OGDF_ASSERT(node->isLeaf());
			leaves_found++;
		}

		OGDF_ASSERT(node->getDegree() >= 1);
		if (node->nodeType != PCNodeType::Leaf && !allow_small_deg) {
			OGDF_ASSERT(node->getDegree() >= 3);
		}

		// also check that all my children know me and that my degree is right
		PCNode* pred = nullptr;
		PCNode* curr = node->child1;
		size_t children = 0;
		while (curr != nullptr) {
			todo.push(curr);
			OGDF_ASSERT(curr->getParent() == node);
			if (node->getNodeType() == PCNodeType::CNode) {
				OGDF_ASSERT(curr->parentPNode == nullptr);
				OGDF_ASSERT(curr->parentCNodeId
						== node->nodeListIndex); // getParent() updates parentCNodeId
			} else {
				OGDF_ASSERT(curr->parentPNode == node);
				OGDF_ASSERT(curr->parentCNodeId == UNIONFINDINDEX_EMPTY);
			}
			children++;
			proceedToNextSibling(pred, curr);
		}
		OGDF_ASSERT(children == node->childCount);
		OGDF_ASSERT(pred == node->child2);
	}
	OGDF_ASSERT(leaves_found == leaves.size());
	OGDF_ASSERT(p_nodes_found == pNodeCount);
	OGDF_ASSERT(c_nodes_found == cNodeCount);

	// list c-nodes
	OGDF_ASSERT(forest->cNodes.size() >= cNodeCount);
	/*c_nodes_found = 0;
    for (int cid = 0; cid < cNodes.size(); cid++) {
        PCNode *node = cNodes.at(cid);
        OGDF_ASSERT(node == nullptr || parents.getRepresentative(cid) == cid);
        if (node == nullptr) continue;
        c_nodes_found++;
        OGDF_ASSERT(id_seen.at(node->id) == node);
    }
    OGDF_ASSERT(c_nodes_found == cNodeCount);*/

	return true;
}

void PCTree::getRestrictions(std::vector<std::vector<PCNode*>>& restrictions,
		PCNode* fixedLeaf) const {
	PCTreeNodeArray<size_t> readyChildren(*this, 0);
	PCTreeNodeArray<std::list<PCNode*>> subtreeLeaves(*this);
	std::queue<PCNode*> todo;
	for (PCNode* leaf : leaves) {
		if (leaf == fixedLeaf) {
			continue;
		}
		subtreeLeaves[leaf].push_back(leaf);
		PCNode* next = leaf == rootNode ? leaf->child1 : leaf->getParent();
		if ((readyChildren[next] += 1) == next->getDegree() - 1) {
			todo.push(next);
		}
	}
	PCNode* central = nullptr;
	while (!todo.empty()) {
		PCNode* node = todo.front();
		todo.pop();
		OGDF_ASSERT(node != fixedLeaf);

		PCNode* next = nullptr;
		PCNode* parent = node->getParent();
		if (parent != nullptr && subtreeLeaves[parent].empty()
				&& !(parent == rootNode && parent->isLeaf())) {
			next = parent;
		}
#ifndef OGDF_DEBUG
		else
#endif
		{
			for (PCNode* neigh : node->neighbors()) {
				if (subtreeLeaves[neigh].empty()) {
					OGDF_ASSERT(next == nullptr || (next == parent && neigh == parent));
					next = neigh;
#ifndef OGDF_DEBUG
					break;
#endif
				}
			}
		}
		if (next == nullptr) {
			OGDF_ASSERT(fixedLeaf == nullptr);
			OGDF_ASSERT(central == nullptr);
			OGDF_ASSERT(todo.empty());
			central = node;
		}

		PCNode* pred = nullptr;
		for (PCNode* curr : node->neighbors(next)) {
			if (curr == next) {
				continue;
			}
			OGDF_ASSERT(!subtreeLeaves[curr].empty());
			if (node->nodeType == PCNodeType::CNode && pred != nullptr) {
				unsigned long size = subtreeLeaves[pred].size() + subtreeLeaves[curr].size();
				if (!isTrivialRestriction(size)) {
					std::vector<PCNode*>& back = restrictions.emplace_back();
					back.reserve(size);
					back.insert(back.end(), subtreeLeaves[pred].begin(), subtreeLeaves[pred].end());
					back.insert(back.end(), subtreeLeaves[curr].begin(), subtreeLeaves[curr].end());
				}
			}
			if (pred != nullptr) {
				subtreeLeaves[node].splice(subtreeLeaves[node].end(), subtreeLeaves[pred]);
			}
			pred = curr;
		}
		if (pred != next) {
			subtreeLeaves[node].splice(subtreeLeaves[node].end(), subtreeLeaves[pred]);
		}

		if (node->nodeType == PCNodeType::PNode && !isTrivialRestriction(subtreeLeaves[node].size())) {
			restrictions.emplace_back(subtreeLeaves[node].begin(), subtreeLeaves[node].end());
		}

		if (next != nullptr && (readyChildren[next] += 1) == next->getDegree() - 1
				&& next != fixedLeaf) {
			todo.push(next);
		}
	}
	if (fixedLeaf != nullptr) {
		OGDF_ASSERT(central == nullptr);
		central = fixedLeaf == rootNode ? fixedLeaf->child1 : fixedLeaf->getParent();
		OGDF_ASSERT(readyChildren[central] == central->getDegree() - 1);
		OGDF_ASSERT(subtreeLeaves[central].size() == getLeafCount() - 1);
	} else {
		OGDF_ASSERT(central != nullptr);
		OGDF_ASSERT(readyChildren[central] == central->getDegree());
		OGDF_ASSERT(subtreeLeaves[central].size() == getLeafCount());
	}
}

template<typename R>
R PCTree::possibleOrders() const {
	R orders(1);
	for (PCNode* node : innerNodes()) {
		if (node->getNodeType() == PCNodeType::CNode) {
			orders *= 2;
		} else {
			R children(node->getChildCount());
			if (node == rootNode) {
				children -= 1; // don't count circular shifts
			}
			orders *= factorial(children);
		}
	}
	return orders;
}

PCNode* PCTree::setRoot(PCNode* newRoot) {
	OGDF_ASSERT(newRoot->isDetached());
	PCNode* oldRoot = rootNode;
	rootNode = newRoot;
	OGDF_ASSERT(checkValid());
	return oldRoot;
}

PCNode* PCTree::changeRoot(PCNode* newRoot) {
	OGDF_ASSERT(checkValid());
	std::stack<PCNode*> path;
	for (PCNode* node = newRoot; node != nullptr; node = node->getParent()) {
		OGDF_ASSERT(node != nullptr);
		OGDF_ASSERT(node->forest == forest);
		path.push(node);
	}
	while (path.size() > 1) {
		PCNode* old_parent = path.top();
		path.pop();
		PCNode* new_parent = path.top();

		if (!old_parent->isLeaf()) {
			PCNode* sib1 = old_parent->getNextNeighbor(nullptr, new_parent);
			PCNode* sib2 = old_parent->getNextNeighbor(sib1, new_parent);
			new_parent->detach();
			old_parent->child1->replaceSibling(nullptr, old_parent->child2);
			old_parent->child2->replaceSibling(nullptr, old_parent->child1);
			old_parent->child1 = sib1;
			old_parent->child2 = sib2;
			sib1->replaceSibling(sib2, nullptr);
			sib2->replaceSibling(sib1, nullptr);
		} else {
			new_parent->detach();
		}
		new_parent->appendChild(old_parent);
	}
	OGDF_ASSERT(path.size() == 1);
	OGDF_ASSERT(path.top() == newRoot);

	return setRoot(newRoot);
}
