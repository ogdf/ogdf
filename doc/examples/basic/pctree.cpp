#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/pctree/PCEnum.h>
#include <ogdf/basic/pctree/PCNode.h>
#include <ogdf/basic/pctree/PCRegistry.h>
#include <ogdf/basic/pctree/PCTree.h>
#include <ogdf/basic/pctree/PCTreeIterators.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/tree/RadialTreeLayout.h>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace ogdf;
using namespace pc_tree;

int main() {
	// create a tree with a given number of leaves
	std::vector<PCNode*> leaves;
	PCTree tree(7, &leaves);

	// the tree allows any order of its leaves
	OGDF_ASSERT(tree.isTrivial());
	OGDF_ASSERT(tree.getLeafCount() == 7);
	OGDF_ASSERT(tree.getPNodeCount() == 1);
	OGDF_ASSERT(tree.getCNodeCount() == 0);
	OGDF_ASSERT(tree.possibleOrders<int>() == factorial(6));
	OGDF_ASSERT(tree.uniqueID(uid_utils::nodeToPosition) == "7:(6, 5, 4, 3, 2, 1, 0)");

	// now we can force some leaves to be consecutive in all represented orders
	tree.makeConsecutive({leaves.at(3), leaves.at(4), leaves.at(5)});
	OGDF_ASSERT(tree.uniqueID(uid_utils::nodeToPosition) == "8:(7:(5, 4, 3), 6, 2, 1, 0)");
	OGDF_ASSERT(tree.getPNodeCount() == 2);
	OGDF_ASSERT(tree.possibleOrders<int>() == factorial(3) * factorial(4));

	// further updates further change the tree
	tree.makeConsecutive({leaves.at(4), leaves.at(5)});
	OGDF_ASSERT(tree.uniqueID(uid_utils::nodeToPosition) == "9:(8:[7:[5, 4], 3], 6, 2, 1, 0)");

	tree.makeConsecutive({leaves.at(0), leaves.at(1)});
	OGDF_ASSERT(tree.uniqueID(uid_utils::nodeToPosition) == "10:(9:[8:[5, 4], 3], 7:[1, 0], 6, 2)");

	tree.makeConsecutive({leaves.at(1), leaves.at(2)});
	OGDF_ASSERT(tree.uniqueID(uid_utils::nodeToPosition) == "10:[9:[8:[5, 4], 3], 7:[2, 1, 0], 6]");

	tree.makeConsecutive({leaves.at(2), leaves.at(3), leaves.at(4), leaves.at(5)});
	OGDF_ASSERT(tree.uniqueID(uid_utils::nodeToPosition) == "9:[6, 8:[7:[5, 4], 3], 2, 1, 0]");

	// if some leaves cannot be made consecutive, the update returns false and leaves the tree unchanged
	OGDF_ASSERT(tree.makeConsecutive({leaves.at(2), leaves.at(0)}) == false);
	OGDF_ASSERT(tree.uniqueID(uid_utils::nodeToPosition) == "9:[6, 8:[7:[5, 4], 3], 2, 1, 0]");
	OGDF_ASSERT(tree.possibleOrders<int>() == 8);

	// we can query the (arbitrary) current order of leaves and check whether any order is represented by the tree
	std::vector<PCNode*> currentLeafOrder =
			tree.currentLeafOrder(); // same entries as `leaves`, different order
	OGDF_ASSERT(tree.isValidOrder(currentLeafOrder));
	// use tree.getLeaves() to get the leaves in creation order

	// iterator for all inner nodes
	auto innerNode_it = tree.innerNodes();
	std::vector<PCNode*> innerNodes {innerNode_it.begin(), innerNode_it.end()};
	OGDF_ASSERT(innerNodes.size() == tree.getPNodeCount() + tree.getCNodeCount());
	OGDF_ASSERT(innerNodes.front() == tree.getRootNode());

	// we can also manually walk the tree
	PCNode* root = tree.getRootNode();
	OGDF_ASSERT(root->getChildCount() == 5); // 6, 8:[...], 2, 1, 0
	OGDF_ASSERT(root->getChild1()->isLeaf());
	for (PCNode* n : root->children()) {
		std::cout << n->index() << " ";
	}
	std::cout << std::endl;

	// we can visualize the PCTree
	ogdf::Graph G;
	ogdf::GraphAttributes GA(G);
	PCTreeNodeArray<ogdf::node> pc_repr(tree);
	tree.getTree(G, &GA, pc_repr);
	RadialTreeLayout l;
	l.rootSelection(RadialTreeLayout::RootSelectionType::Source);
	l.call(GA);
	GraphIO::write(GA, "pctree.svg");

	// PCTrees can also be reconstructed from strings
	PCTree from_string {"9:[6, 8:[7:[5, 4], 3], 2, 1, 0]", true};
	OGDF_ASSERT(from_string.getLeafCount() == tree.getLeafCount());
	OGDF_ASSERT(from_string.possibleOrders<int>() == tree.possibleOrders<int>());
	OGDF_ASSERT(from_string.uniqueID(uid_utils::nodeToPosition) == "9:[6, 8:[7:[5, 4], 3], 2, 1, 0]");
	// tree.getRestrictions() yields a list of updated via which the tree can also be reconstructed

	// alternatively, they can also be constructed manually
	PCTree manual;
	auto nr = manual.newNode(pc_tree::PCNodeType::CNode);
	auto n1 = manual.newNode(pc_tree::PCNodeType::PNode, nr);
	manual.insertLeaves(3, nr);
	auto n2 = manual.newNode(pc_tree::PCNodeType::PNode, nr);
	manual.insertLeaves(3, n2);
	manual.insertLeaves(3, nr);
	manual.insertLeaves(3, n1);
	std::stringstream ss;
	ss << manual;
	OGDF_ASSERT(ss.str() == "0:[11, 10, 9, 5:(8, 7, 6), 4, 3, 2, 1:(14, 13, 12)]");

	return 0;
}
