#pragma once

#include <ogdf/basic/Graph.h>
#include <ogdf/cluster/sync_plan/utils/Bijection.h>

#include <ostream>

using namespace ogdf;

struct EncapsulatedBlock {
	node bicon = nullptr;
	node bicon_rep = nullptr;
	node star_rep = nullptr;
	PipeBij bij;

	explicit EncapsulatedBlock(node _bicon) : bicon(_bicon) { }

	friend std::ostream& operator<<(std::ostream& os, const EncapsulatedBlock& block);
};

std::ostream& operator<<(std::ostream& os, const EncapsulatedBlock& block);
