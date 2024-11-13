/** \file
 * \brief Implements graph generator for hierarchical graphs.
 *
 * \author Carsten Gutwenger, Christoph Buchheim, Simon D. Fink
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.md in the OGDF root directory for details.
 *
 * \par
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * Version 2 or 3 as published by the Free Software Foundation;
 * see the file LICENSE.txt included in the packaging of this file
 * for details.
 *
 * \par
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * \par
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, see
 * http://www.gnu.org/copyleft/gpl.html
 */


#include <ogdf/basic/Array.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/comparer.h>
#include <ogdf/basic/graph_generators/deterministic.h>
#include <ogdf/basic/graph_generators/randomized.h>
#include <ogdf/basic/memory.h>

#include <algorithm>
#include <random>
#include <vector>

using std::minstd_rand;
using std::mt19937;
using std::uniform_int_distribution;
using std::uniform_real_distribution;
using std::vector;

namespace ogdf {


class BEdge {
public:
	int head, tail, id, pos;
	BEdge* next;

	BEdge(int t, int h, int c) : head(h), tail(t), id(c), pos(-1), next(nullptr) { }

	OGDF_NEW_DELETE
};

using bEdge = BEdge*;


OGDF_DECLARE_COMPARER(CmpTail, bEdge, int, x->tail);
OGDF_DECLARE_COMPARER(CmpHead, bEdge, int, x->head);

void randomHierarchy(Graph& G, int numberOfNodes, int numberOfEdges, bool planar, bool singleSource,
		bool longEdges) {
	Array<node> nnr(3 * numberOfNodes);
	Array<int> vrt(3 * numberOfNodes);
	Array<int> fst(numberOfNodes + 1);

	/** Place nodes **/

	emptyGraph(G, numberOfNodes);

	minstd_rand rng(randomSeed());
	uniform_real_distribution<> dist_0_1(0.0, 1.0);

	int numberOfLayers = 0, totNumber = 0, realCount = 0;
	fst[0] = 0;
	for (node v : G.nodes) {
		if (longEdges && numberOfLayers) {
			vrt[totNumber++] = 1;
		}

		nnr[totNumber] = v;
		vrt[totNumber++] = 0;
		realCount++;
		double r = dist_0_1(rng);
		if ((totNumber == 1 && singleSource) || realCount == numberOfNodes
				|| r * r * numberOfNodes < 1) {
			if (longEdges && numberOfLayers) {
				vrt[totNumber++] = 1;
			}
			fst[++numberOfLayers] = totNumber;
		}
	}

	/** Determine allowed neighbours **/

	Array<int> leftN(totNumber);
	Array<int> rightN(totNumber);
	for (int layer = 1; layer < numberOfLayers; layer++) {
		if (planar) {
			int n1 = fst[layer - 1];
			int n2 = fst[layer];
			leftN[n2] = n1;
			while (n1 < fst[layer] && n2 < fst[layer + 1]) {
				double r = dist_0_1(rng);
				if (n1 != fst[layer] - 1
						&& (n2 == fst[layer + 1] - 1
								|| r < (double)(fst[layer] - fst[layer - 1])
												/ (double)(fst[layer + 1] - fst[layer - 1]))) {
					n1++;
				} else {
					rightN[n2] = n1;
					if (++n2 < fst[layer + 1]) {
						leftN[n2] = n1;
					}
				}
			}
		} else {
			for (int n2 = fst[layer]; n2 < fst[layer + 1]; n2++) {
				leftN[n2] = fst[layer - 1];
				rightN[n2] = fst[layer] - 1;
			}
		}
	}

	/** Insert edges **/

	List<bEdge> startEdges;
	Array<SList<bEdge>> edgeIn(totNumber);
	Array<SList<bEdge>> edgeOut(totNumber);

	if (numberOfLayers) {
		double x1 = numberOfEdges;
		double x2 = 0;
		for (int n2 = fst[1]; n2 < totNumber; n2++) {
			if (!vrt[n2]) {
				x2 += rightN[n2] - leftN[n2] + 1;
			}
		}

		int idc = 0;
		for (int n2 = fst[1]; n2 < totNumber; n2++) {
			if (!vrt[n2]) {
				bool connected = !singleSource;
				for (int n1 = leftN[n2]; n1 <= rightN[n2] || !connected; n1++) {
					double r = dist_0_1(rng);
					if (r < x1 / x2 || n1 > rightN[n2]) {
						int next = (n1 <= rightN[n2]
										? n1
										: uniform_int_distribution<>(leftN[n2], rightN[n2])(rng));
						int act = n2;
						bEdge nextEdge = new BEdge(next, act, idc++);
						while (vrt[next]) {
							act = next;
							next = uniform_int_distribution<>(leftN[act], rightN[act])(rng);
							edgeOut[act].pushBack(nextEdge);
							nextEdge = new BEdge(next, act, idc++);
							edgeIn[act].pushBack(nextEdge);
						}
						startEdges.pushBack(nextEdge);
						connected = true;
						x1 -= 1;
					}
					if (n1 <= rightN[n2]) {
						x2 -= 1;
					}
				}
			}
		}
	}

	if (planar) {
		for (int act = 0; act < totNumber; act++) {
			CmpTail cmpTail;
			edgeIn[act].quicksort(cmpTail);
			CmpHead cmpHead;
			edgeOut[act].quicksort(cmpHead);
		}
	}

	for (int act = 0; act < totNumber; act++) {
		for (bEdge nextEdge : edgeIn[act]) {
			nextEdge->next = edgeOut[act].popFrontRet();
		}
	}

	for (bEdge actEdge : startEdges) {
		bEdge nextEdge = actEdge;
		while (vrt[nextEdge->head]) {
			nextEdge = nextEdge->next;
		}
		G.newEdge(nnr[actEdge->tail], nnr[nextEdge->head]);
	}

	/** Clean up **/
	for (bEdge nextEdge : startEdges) {
		bEdge toDelete = nextEdge;
		while (vrt[nextEdge->head]) {
			nextEdge = nextEdge->next;
			delete toDelete;
			toDelete = nextEdge;
		}
		delete toDelete;
	}
}

void randomProperMaximalLevelPlaneGraph(Graph& G, std::vector<std::vector<node>>& emb, int N, int K,
		bool radial) {
	// checks
	OGDF_ASSERT(N > 0);
	OGDF_ASSERT(K > 0);

	// init arrays
	G.clear();
	emb.clear();
	emb.resize(K);

	// create nodes, at least one on each level
	for (int i = 0; i < K; ++i) {
		emb[i].push_back(G.newNode());
	}
	for (int i = K; i < N; ++i) {
		emb[randomNumber(0, K - 1)].push_back(G.newNode());
	}

	// create all edges on each level
	for (int l = 0; l < K - 1; ++l) {
		int lp = 0;
		int up = 0;
		bool wrap_forw = false; // whether there is an edge from emb[l].back() to emb[l+1].front()
		bool wrap_back = false; // whether there is an edge from emb[l].front() to emb[l+1].back()
		while (true) {
			G.newEdge(emb[l][lp], emb[l + 1][up]);

			wrap_forw = wrap_forw || (lp + 1 == emb[l].size() && up == 0);
			wrap_back = wrap_back || (lp == 0 && up + 1 == emb[l + 1].size());

			// select the next node either on the upper or the lower level
			if (lp + 1 < emb[l].size() && up + 1 < emb[l + 1].size()) {
				int r = randomNumber(0, 1);
				lp += r;
				up += 1 - r;
			} else if (lp + 1 < emb[l].size()) {
				lp++;
			} else if (up + 1 < emb[l + 1].size()) {
				up++;
			} else {
				break;
			}
		}

		// add the wrap-around edge
		if (radial) {
			if (wrap_back || (!wrap_forw && randomNumber(0, 1) == 1)) {
				if (!wrap_forw) {
					G.newEdge(emb[l].back(), emb[l + 1].front());
				}
			} else {
				OGDF_ASSERT(!wrap_back);
				G.newEdge(emb[l].front(), emb[l + 1].back());
			}
		}
	}
}

void pruneEdges(Graph& G, int max_edges, int min_deg) {
	vector<edge> edges;
	for (edge e : G.edges) {
		edges.push_back(e);
	}
	mt19937 mt(randomSeed());
	shuffle(edges.begin(), edges.end(), mt);
	for (edge e : edges) {
		if (e->source()->degree() > min_deg && e->target()->degree() > min_deg) {
			G.delEdge(e);
		}
		if (G.numberOfEdges() <= max_edges) {
			break;
		}
	}
}

}
