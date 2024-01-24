#pragma once

#include <ogdf/basic/DisjointSets.h>
#include <ogdf/basic/pctree/PCEnum.h>
#include <ogdf/basic/pctree/PCRegistry.h>

#include <cstdint>
#include <vector>

#define PCTREE_REUSE_NODES

#define UNIONFINDINDEX_EMPTY (-1)

namespace pc_tree {
using UnionFindIndex = int;

class PCTreeForest {
	friend class PCNode;
	friend class PCTree;

	template<class Key>
	friend class PCTreeRegistry;

private:
	std::vector<PCTree*> trees;
	std::vector<PCNode*> cNodes;
	ogdf::DisjointSets<> parents {1 << 8};
	int nextNodeId = 0;
	int timestamp = 0;
	PCTreeRegistry<PCNode*> nodeArrayRegistry;
	bool autodelete;

#ifdef PCTREE_REUSE_NODES
	// TODO: also reuse PCTrees?
	PCNode* reusableNodes = nullptr;
#endif

public:
	PCTreeForest(bool p_autodelete = true) : nodeArrayRegistry(this), autodelete(p_autodelete) {};

	virtual ~PCTreeForest();

	PCTree* makeTree(void);

	bool merge(PCTree* a, PCTree* b);

	void clear(void);

	operator const PCTreeRegistry<PCNode*>&() const { return nodeArrayRegistry; }
};
}
