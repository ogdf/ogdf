#pragma once

#include <ogdf/basic/RegisteredArray.h>
#include <ogdf/basic/RegisteredSet.h>

#include <ostream>

namespace pc_tree {
enum class NodeLabel { Unknown, Partial, Full, Empty = Unknown };

enum class PCNodeType { PNode, CNode, Leaf };

class PCTree;

class PCTreeForest;

template<class Key>
class PCTreeRegistry;

class PCNode;

#define OGDF_DECL_REG_ARRAY_TYPE(v, c) ogdf::RegisteredArray<PCTreeRegistry<PCNode*>, v, c>
OGDF_DECL_REG_ARRAY(PCTreeNodeArray)
#undef OGDF_DECL_REG_ARRAY_TYPE

template<bool SupportFastSizeQuery = true>
using PCTreeNodeSet = ogdf::RegisteredSet<PCTreeRegistry<PCNode*>, SupportFastSizeQuery>;
}

std::ostream& operator<<(std::ostream&, pc_tree::NodeLabel);

std::ostream& operator<<(std::ostream&, pc_tree::PCNodeType);

std::ostream& operator<<(std::ostream&, const pc_tree::PCTree*);

std::ostream& operator<<(std::ostream&, const pc_tree::PCNode*);

std::ostream& operator<<(std::ostream&, const pc_tree::PCTree&);

std::ostream& operator<<(std::ostream&, const pc_tree::PCNode&);
