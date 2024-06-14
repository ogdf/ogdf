#pragma once

#include <ogdf/basic/Graph.h>

#include "Iterators.h"

using namespace ogdf;

using PipeBijIterator = ZipIterator<adjEntry, internal::GraphObjectContainer<AdjElement>::iterator,
		adjEntry, internal::GraphObjectContainer<AdjElement>::reverse_iterator>;
using PipeBijPair = std::pair<adjEntry, adjEntry>;
using FrozenPipeBijPair = std::pair<int, int>;
using PipeBij = List<PipeBijPair>; // TODO convert to vector?
using FrozenPipeBij = List<FrozenPipeBijPair>;

OGDF_DECLARE_COMPARER(PipeBijCmp, PipeBijPair, int, x.first->theEdge()->index());

OGDF_DECLARE_COMPARER(FrozenPipeBijCmp, FrozenPipeBijPair, int, x.first);

PipeBijIterator getPipeBijection(node u, node v);

void getPipeBijection(node u, node v, PipeBij& out);

void getPipeBijection(node u, node v, AdjEntryArray<adjEntry>& out);

void getPipeBijection(node u, node v, EdgeArray<edge>& out);

void getFrozenPipeBijection(node u, node v, FrozenPipeBij& out);

void freezePipeBijection(const PipeBij& in, FrozenPipeBij& out);
