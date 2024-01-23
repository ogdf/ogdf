#include <bandit/bandit.h>

#include "PCNode.h"
#include "PCTree.h"

using namespace pc_tree;
using namespace snowhouse;
using namespace bandit;

std::string treeToString(const PCTree& tree) {
    std::stringstream ss;
    ss << tree;
    return ss.str();
}

template<class Iterable>
void checkOrder(std::vector<PCNode*>& nodes, Iterable iter, std::vector<int> ids) {
    int index = 0;
    for (PCNode* node : iter) {
        AssertThat(node, Equals(nodes.at(ids.at(index))));
        index++;
    }
    AssertThat(index, Equals(ids.size()));
}

go_bandit([]() {
    describe("PCTree construction", []() {
        it("allows creating a trivial instance", []() {
            std::vector<PCNode*> leaves;
            PCTree tree(5, &leaves);
            AssertThat(tree.checkValid(), IsTrue());
            AssertThat(tree.isTrivial(), IsTrue());
            AssertThat(tree.getLeafCount(), Equals(5));
            AssertThat(tree.getPNodeCount(), Equals(1));
            AssertThat(tree.getCNodeCount(), Equals(0));
            AssertThat(std::equal(tree.getLeaves().begin(), tree.getLeaves().end(), leaves.begin(),
                               leaves.end()),
                    IsTrue());
            AssertThat(tree.possibleOrders(), Equals(Dodecahedron::Bigint(24)));
            AssertThat(treeToString(tree), Equals("0:(5, 4, 3, 2, 1)"));
            AssertThat(tree.uniqueID(uid_utils::nodeToPosition), Equals("5:(4, 3, 2, 1, 0)"));
            PCNode* root = tree.getRootNode();
            AssertThat(root->getNodeType(), Equals(PCNodeType::PNode));
            AssertThat(root->getChildCount(), Equals(5));

            tree.makeConsecutive({leaves.at(1), leaves.at(2)});
            AssertThat(treeToString(tree), Equals("0:(6:(3, 2), 5, 4, 1)"));
            AssertThat(tree.uniqueID(uid_utils::nodeToPosition), Equals("6:(5:[2, 1], 4, 3, 0)"));
            AssertThat(tree.checkValid(), IsTrue());
            AssertThat(tree.isTrivial(), IsFalse());
        });

        it("correctly handles the first case where JPPZanetti fails", []() {
            std::vector<PCNode*> leaves;
            PCTree tree(7, &leaves);
            AssertThat(tree.checkValid(), IsTrue());
            AssertThat(tree.uniqueID(uid_utils::nodeToPosition), Equals("7:(6, 5, 4, 3, 2, 1, 0)"));
            tree.makeConsecutive({leaves.at(4), leaves.at(5)});
            AssertThat(tree.checkValid(), IsTrue());
            AssertThat(tree.uniqueID(uid_utils::nodeToPosition),
                    Equals("8:(7:[5, 4], 6, 3, 2, 1, 0)"));
            tree.makeConsecutive({leaves.at(3), leaves.at(4), leaves.at(5)});
            AssertThat(tree.checkValid(), IsTrue());
            AssertThat(tree.uniqueID(uid_utils::nodeToPosition),
                    Equals("9:(8:[7:[5, 4], 3], 6, 2, 1, 0)"));
            tree.makeConsecutive({leaves.at(0), leaves.at(1)});
            AssertThat(tree.checkValid(), IsTrue());
            AssertThat(tree.uniqueID(uid_utils::nodeToPosition),
                    Equals("10:(9:[8:[5, 4], 3], 7:[1, 0], 6, 2)"));
            tree.makeConsecutive({leaves.at(1), leaves.at(2)});
            AssertThat(tree.checkValid(), IsTrue());
            AssertThat(tree.uniqueID(uid_utils::nodeToPosition),
                    Equals("10:[9:[8:[5, 4], 3], 7:[2, 1, 0], 6]"));
            tree.makeConsecutive({leaves.at(2), leaves.at(3), leaves.at(4), leaves.at(5)});
            AssertThat(tree.checkValid(), IsTrue());
            AssertThat(tree.uniqueID(uid_utils::nodeToPosition),
                    Equals("9:[6, 8:[7:[5, 4], 3], 2, 1, 0]"));
            tree.makeConsecutive({leaves.at(3), leaves.at(4), leaves.at(5), leaves.at(6)});
            AssertThat(tree.checkValid(), IsTrue());
            AssertThat(tree.uniqueID(uid_utils::nodeToPosition),
                    Equals("9:[6, 8:[7:[5, 4], 3], 2, 1, 0]"));
            tree.makeConsecutive({leaves.at(3), leaves.at(4)});
            AssertThat(tree.checkValid(), IsTrue());
            AssertThat(tree.uniqueID(uid_utils::nodeToPosition),
                    Equals("8:[6, 7:[5, 4, 3], 2, 1, 0]"));
            tree.makeConsecutive({leaves.at(2), leaves.at(3), leaves.at(4), leaves.at(5)});
            AssertThat(tree.checkValid(), IsTrue());
            AssertThat(tree.uniqueID(uid_utils::nodeToPosition),
                    Equals("8:[6, 7:[5, 4, 3], 2, 1, 0]"));
        });

        it("is iterated correctly", []() {
            std::vector<PCNode*> leaves;
            PCTree tree(7, &leaves);
            AssertThat(tree.checkValid(), IsTrue());
            AssertThat(treeToString(tree), Equals("0:(7, 6, 5, 4, 3, 2, 1)"));
            tree.makeConsecutive({leaves.at(0), leaves.at(1), leaves.at(2), leaves.at(3)});
            AssertThat(tree.checkValid(), IsTrue());
            AssertThat(treeToString(tree), Equals("0:(8:(4, 3, 2, 1), 7, 6, 5)"));

            FilteringPCTreeDFS walk = tree.allNodes();
            std::vector<PCNode*> nodes {walk.begin(), walk.end()};
            std::sort(nodes.begin(), nodes.end(),
                    [](PCNode* a, PCNode* b) { return a->index() < b->index(); });
            AssertThat(nodes.size(), Equals(9));
            AssertThat(nodes.front()->index(), Equals(0));
            AssertThat(nodes.back()->index(), Equals(8));

            checkOrder(nodes, tree.getLeaves(), {1, 2, 3, 4, 5, 6, 7});
            checkOrder(nodes, tree.allNodes(), {0, 8, 4, 3, 2, 1, 7, 6, 5});
            checkOrder(nodes, tree.innerNodes(), {0, 8});

            FilteringPCTreeBFS bfs(tree, tree.getRootNode());
            checkOrder(nodes, bfs, {0, 5, 6, 7, 8, 1, 2, 3, 4});

            PCNode* node = tree.getRootNode()->getChild2();
            AssertThat(node->index(), Equals(8));
            checkOrder(nodes, node->children(), {1, 2, 3, 4});
            checkOrder(nodes, node->neighbors(), {1, 2, 3, 4, 0});
        });
    });
});

int main(int argc, char* argv[]) {
    return bandit::run(argc, argv);
}