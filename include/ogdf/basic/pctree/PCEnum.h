#pragma once

#include <ostream>

#include "utils/RegisteredArray.h"
#include "utils/RegisteredElementSet.h"

namespace pc_tree {
    enum class NodeLabel { Unknown, Partial, Full, Empty = Unknown };

    enum class PCNodeType { PNode, CNode, Leaf };

    class PCTree;

    class PCTreeForest;

    template<class Key>
    class PCTreeRegistry;

    class PCNode;

    template<typename Value>
    using PCTreeNodeArray = ogdf::RegisteredArray<PCTreeRegistry<PCNode*>, PCNode*, Value>;

    using PCTreeNodeSet = ogdf::RegisteredElementSet<PCNode*, PCTreeRegistry<PCNode*>>;
}

std::ostream& operator<<(std::ostream&, pc_tree::NodeLabel);

std::ostream& operator<<(std::ostream&, pc_tree::PCNodeType);

std::ostream& operator<<(std::ostream&, const pc_tree::PCTree*);

std::ostream& operator<<(std::ostream&, const pc_tree::PCNode*);

std::ostream& operator<<(std::ostream&, const pc_tree::PCTree&);

std::ostream& operator<<(std::ostream&, const pc_tree::PCNode&);
