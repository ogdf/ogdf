/** \file
 * \author Simon D. Fink <ogdf@niko.fink.bayern>
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
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/Logger.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/basic.h>
#include <ogdf/cluster/ClusterGraph.h>
#include <ogdf/cluster/sync_plan/utils/Preprocess.h>

#include <ostream>

namespace ogdf::sync_plan {

using namespace preprocess;

bool preprocessClusterGraph(ClusterGraph& C, Graph& G) {
	bool modified = false;
	bool progress = true;
	while (progress) {
		progress = removeSmallNodes(C, G) || unsplitDeg2Nodes(C, G)
				|| disconnectedClustersToStars(C, G) || removeSmallClusters(C, G);
		modified = modified || progress;
	}
	return modified;
}

bool canPreprocessClusterGraph(const ClusterGraph& C, const Graph& G) {
	return !findSmallNodes(C, G).empty() || !findDeg2Nodes(C, G).empty()
			|| !findDisconnectedClusters(C, G).empty() || !findSmallClusters(C, G).empty();
}

namespace preprocess {
ogdf::Logger preprocessLog;

SList<node> findSmallNodes(const ClusterGraph& C, const Graph& G) {
	SList<node> toRemove;
	for (node v : G.nodes) {
		if (v->degree() == 0) {
			toRemove.pushBack(v);
		} else if (v->degree() == 1) {
			node u = v->firstAdj()->twinNode();
			cluster cv = C.clusterOf(v);
			cluster cu = C.clusterOf(u);

			if (cv == cu) {
				toRemove.pushBack(v);
			} else {
				for (adjEntry adj : u->adjEntries) {
					node w = adj->twinNode();
					cluster cw = C.clusterOf(w);
					if (w == v || cw == cu) {
						continue;
					}

					List<cluster> path;
					C.commonClusterPath(u, w, path);
					if (path.search(cv).valid()) {
						toRemove.pushBack(v);
						break;
					}
				}
			}
		}
	}
	return toRemove;
}

bool removeSmallNodes(const ClusterGraph& C, Graph& G) {
	SList<node> toRemove = findSmallNodes(C, G);
	preprocessLog.lout() << "Remove " << toRemove.size() << " deg-0/1 vertices" << std::endl;
	for (node vDel : toRemove) {
		G.delNode(vDel);
	}
	return !toRemove.empty();
}

SList<node> findDeg2Nodes(const ClusterGraph& C, const Graph& G) {
	SList<node> toRemove;
	NodeArray<bool> marked(G, false); // mark nodes that shall be unsplit
	for (node v : G.nodes) {
		if (v->degree() != 2) {
			continue;
		}
		cluster cv = C.clusterOf(v);

		// the two neighbors
		node u = v->firstAdj()->twinNode();
		node w = v->lastAdj()->twinNode();

		if (marked[u] || marked[w]) {
			continue;
		}

		List<cluster> path;
		C.commonClusterPath(u, w, path);
		if (path.search(cv).valid()) {
			marked[v] = true;
			toRemove.pushBack(v);
		}
	}
	return toRemove;
}

bool unsplitDeg2Nodes(const ClusterGraph& C, Graph& G) {
	SList<node> toRemove = findDeg2Nodes(C, G);
	preprocessLog.lout() << "Unsplit " << toRemove.size() << " deg-2 vertices" << std::endl;
	for (node vDel : toRemove) {
		node u = vDel->firstAdj()->twinNode();
		node w = vDel->lastAdj()->twinNode();
		if (u != w && G.searchEdge(u, w) == nullptr) {
			if (vDel->firstAdj()->isSource()) {
				G.reverseEdge(vDel->firstAdj()->theEdge());
			}
			if (!vDel->lastAdj()->isSource()) {
				G.reverseEdge(vDel->lastAdj()->theEdge());
			}
			G.unsplit(vDel);
		} else {
			G.delNode(vDel);
		}
	}
	return !toRemove.empty();
}

SList<cluster> findDisconnectedClusters(const ClusterGraph& C, const Graph& G,
		ClusterArray<node>* centre) {
	SList<cluster> toRemove;
	for (cluster c : C.clusters) {
		if (c->cCount() > 0 || c->nCount() < 3) {
			continue;
		}

		bool replaceByStar = true;
		node w = nullptr; // node with deg > 1
		for (node v : c->nodes) {
			for (adjEntry adj : v->adjEntries) {
				if (C.clusterOf(adj->twinNode()) == c) {
					replaceByStar = false;
					break;
				}
			}
			if (v->degree() > 1) {
				if (w == nullptr) {
					w = v;
				} else {
					replaceByStar = false;
				}
			}
			if (!replaceByStar) {
				break;
			}
		}

		if (replaceByStar) {
			if (w == nullptr) {
				w = *c->nBegin();
			}
			if (centre) {
				(*centre)[c] = w;
			}
			toRemove.pushBack(c);
		}
	}
	return toRemove;
}

bool disconnectedClustersToStars(ClusterGraph& C, Graph& G) {
	ClusterArray<node> centre(C, nullptr);
	SList<cluster> toRemove = findDisconnectedClusters(C, G, &centre);
	preprocessLog.lout() << "Replace " << toRemove.size() << " clusters by stars" << std::endl;
	for (cluster c : toRemove) {
		node w = centre[c];
		for (node v : c->nodes) {
			if (v != w) {
				OGDF_ASSERT(v->degree() == 1);
				G.newEdge(v, w);
			}
		}
		C.delCluster(c);
	}
	return !toRemove.empty();
}

SList<cluster> findSmallClusters(const ClusterGraph& C, const Graph& G) {
	SList<cluster> toRemove;
	cluster r = C.rootCluster();
	for (cluster c : C.clusters) {
		if (c == r) {
			continue;
		}
		if ((c->nCount() == 2 && c->cCount() == 0) || (c->cCount() + c->nCount() <= 1)
				|| (c->parent() == r && r->cCount() + r->nCount() == 1)) {
			toRemove.pushBack(c);
		}
	}
	return toRemove;
}

bool removeSmallClusters(ClusterGraph& C, Graph& G) {
	SList<cluster> toRemove = findSmallClusters(C, G);
	preprocessLog.lout() << "Remove " << toRemove.size() << " small clusters" << std::endl;
	for (cluster c : toRemove) {
		if (c->nCount() == 2 && c->cCount() == 0) {
			node v = *c->nBegin();
			node w = *c->nBegin().succ();
			if (G.searchEdge(v, w) == nullptr) {
				if (C.adjAvailable() && c->adjEntries.size() >= 2) {
					adjEntry prev = c->adjEntries.front();
					adjEntry succ = nullptr;
					for (adjEntry adj : c->adjEntries) {
						if (prev->theNode() != adj->theNode()) {
							succ = adj;
							break;
						}
						prev = adj;
					}
					if (succ) {
						G.newEdge(prev, succ->cyclicPred(), Direction::after);
					} else {
						G.newEdge(v, w);
					}
				} else {
					G.newEdge(v, w);
				}
			}
		}
		C.delCluster(c);
	}
	return !toRemove.empty();
}

}
}
