/** \file
 * \brief Calculate minimum cut value for a given Graph.
 *
 * \author Sascha Alder
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

#include <ogdf/graphalg/MinimumCutNagamochiIbaraki.h>

namespace ogdf {

MinimumCutNagamochiIbaraki::MinimumCutNagamochiIbaraki(bool p, bool preprocessing, Level logLevel)
	: Logger(LogMode::Log, logLevel) {
	pr = p;
	m_preprocess = preprocessing;
#ifdef OGDF_DEBUG
	lout(Level::Minor) << "MinimumCutNagamochiIbaraki init done" << std::endl;
#endif
}

MinimumCutNagamochiIbaraki::~MinimumCutNagamochiIbaraki() {
#ifdef OGDF_DEBUG
	lout(Level::High) << "MinimumCutNagamochiIbaraki done" << std::endl;
	lout(Level::High) << std::endl;
	if (barLambda > 1 && size > 2) {
		lout(Level::High) << "Remaining Nodes: ";
		for (const auto& v : allNodes) {
			lout(Level::High) << v->index() << " ";
		}
	}
#endif
}

void MinimumCutNagamochiIbaraki::init(const Graph& G) {
	m_GC.init(G);
	N = m_GC.numberOfNodes();
	M = m_GC.numberOfEdges();
	hiddenEdges = new Graph::HiddenEdgeSet(m_GC);
	size = m_GC.numberOfNodes();
	edgeCapacity.resize(M);
	degree.resize(N);
	setid.resize(N);
	allNodes.reserve(N);
	for (const auto& v : m_GC.nodes) {
		allNodes.insert(v);
	}
}

void MinimumCutNagamochiIbaraki::contractClusters(const std::vector<clusterstruct>& clusters) {
	//contract each cluster (0...level-1)
	//in each round at least one cluster exists  --> otherwise the algorithm is wrong
	//  level is level of lastcluster +1 --> level-1 is last cluster (most right node in mao)
	for (const auto& cl : clusters) {
#ifdef OGDF_DEBUG
		lout(Level::Minor) << "CLUSTER ";
		for (auto& v : cl.clusterNodes) {
			lout(Level::Minor) << v->index() << " ";
		}
		lout(Level::Minor) << std::endl;
		lout(Level::Minor) << "  S: " << cl.clusterhead->index() << std::endl;
#endif
		//compute neighbourhood
		//faster than always trying to use searchEdge()
		const node& s = cl.clusterhead;
		//contract nodes
		contract(s, cl.clusterNodes, setid[s->index()], clusters);
	}

	return;
}

//safe for Nagamochi-Ibaraki
void MinimumCutNagamochiIbaraki::contract(const node& s, const ListPure<node>& toContract,
		const int& clusterLevel, const std::vector<clusterstruct>& clusters) {
	const auto& s_index = s->index();
	auto& s_degree = degree[s_index];

	std::unordered_map<node, edge> neigh;
	neigh.reserve(size / 4);

	ListPure<edge> edgeList;
	s->adjEdges(edgeList);
	node opposite;
	for (edge e : edgeList) {
		opposite = e->opposite(s);
		if (setid[opposite->index()] == clusterLevel) {
			//nodes are in the same cluster
			s_degree = s_degree - edgeCapacity[e->index()];
			hiddenEdges->hide(e);
#ifdef OGDF_DEBUG
			lout() << " del (" << e->source()->index() << "," << e->target()->index()
				   << ") new degree " << s_degree << std::endl;
			lout() << " node " << opposite->index() << " setid " << setid[opposite->index()]
				   << std::endl;
#endif
		} else {
			//opposite is in some other cluster
			neigh.insert(std::make_pair(opposite, e));
		}
	}
	edge e;
	adjEntry succ, adj;

	std::unordered_map<node, edge>::iterator opposite_edge, pi_edge;
	for (node t : toContract) {
#ifdef OGDF_DEBUG
		lout() << "  contract " << s->index() << " " << t->index() << std::endl;
#endif

		adj = t->firstAdj();
		while (adj != nullptr) {
			succ = adj->succ();
			e = getAdjEdge(adj, t, opposite);

#ifdef OGDF_DEBUG
			lout() << "    looking at edge (" << e->source()->index() << "," << e->target()->index()
				   << ") with cap: " << edgeCapacity[e->index()] << std::endl;
#endif
			const auto& opposite_index = opposite->index();
			const auto& opposite_level = setid[opposite_index];
			if (opposite_level == -1) {
				//opposite doesnt have cluster
				s_degree += edgeCapacity[e->index()];
#ifdef OGDF_DEBUG
				lout() << "    -- OPPOSITE HAS NO CLUSTER --" << std::endl;
				lout() << "      adding edge (" << s->index() << "," << opposite->index()
					   << ") with cap: " << edgeCapacity[e->index()] << std::endl;
				lout() << "      new degree " << s_degree << " of node " << s->index() << std::endl;
#endif
				opposite_edge = neigh.find(opposite);
				if (opposite_edge != neigh.end()) {
					//Edge Exists
					edgeCapacity[(opposite_edge->second)->index()] += edgeCapacity[e->index()];
#ifdef OGDF_DEBUG
					lout(Level::Minor) << "     edge (" << opposite_edge->second->source() << ","
									   << opposite_edge->second->target() << ") cap "
									   << edgeCapacity[opposite_edge->second->index()] << std::endl;
#endif
					//hide (opposite,t)
					hiddenEdges->hide(e);
				} else {
					// no edge (s,opposite) exists  -> just change t to s from (t,opposite)
					if (e->source() == t) {
						m_GC.moveSource(e, s);
					} else {
						m_GC.moveTarget(e, s);
					}
#ifdef OGDF_DEBUG
					lout() << "     create edge (" << e->source()->index() << ","
						   << e->target()->index() << ")" << std::endl;
#endif
					neigh.insert(std::make_pair(opposite, e));
				}
			} else if (opposite_level != clusterLevel) {
				// opposite is some other cluster
				s_degree += edgeCapacity[e->index()];
				const auto& pi = clusters[opposite_level].clusterhead;
#ifdef OGDF_DEBUG
				lout() << "  -- OPPOSITE IN OTHER CLUSTER --" << std::endl;
				lout() << "      adding edge (" << s->index() << "," << pi->index() << ")"
					   << std::endl;
				lout() << "      node " << s->index() << " degree " << s_degree << std::endl;
#endif

				if (pi != opposite) {
					//change degree of cluster head if the opposite is in the cluster
					degree[pi->index()] += edgeCapacity[e->index()];
				}
				pi_edge = neigh.find(pi);
				if (pi_edge != neigh.end()) {
					//Edge (s,pi) exists
					edgeCapacity[(pi_edge->second)->index()] += edgeCapacity[e->index()];
#ifdef OGDF_DEBUG
					lout() << "     edge e (" << e->source()->index() << "," << e->target()->index()
						   << ") with capacity " << edgeCapacity[e->index()] << std::endl;
					lout() << "     edge pi (" << (pi_edge->second)->source()->index() << ","
						   << (pi_edge->second)->target()->index() << ") with capacity "
						   << edgeCapacity[e->index()] << std::endl;
					lout() << "     node " << pi << " degree " << degree[pi->index()] << std::endl;
#endif
					//hide o(t,pi)
					hiddenEdges->hide(e);
				} else {
					//edge doesnt exist -> create new
					m_GC.moveSource(e, s);
					m_GC.moveTarget(e, pi);
#ifdef OGDF_DEBUG
					lout() << "     create  new edge (" << e->source()->index() << ","
						   << e->target()->index() << ") with index " << e->index() << "/"
						   << edgeCapacity.capacity() << std::endl;
					lout() << "       capacity " << edgeCapacity[e->index()] << std::endl;
					lout() << "     node " << pi->index() << " degree " << degree[pi->index()]
						   << std::endl;
#endif
					neigh.insert(std::make_pair(pi, e));
				}
			} else {
				//same cluster ->hide
				hiddenEdges->hide(e);
			}
			adj = succ;
		}
		size--;
		allNodes.erase(t);
	}
	updateLambda(s_degree);

	//reset setid
	setid[s->index()] = -1;
}

void MinimumCutNagamochiIbaraki::updateClusters(const node& head, const node& currentNode,
		std::vector<clusterstruct>& clusters, int& level) {
	//add to cluster
	if (setid[head->index()] == -1) {
		//cluster doesnt exist -> need new one
		setid[head->index()] = level;

		//create new cluster
		clusters.emplace_back(clusterstruct {ListPure<node>(), head});

		//add node to cluster
		(clusters.rbegin())->clusterNodes.emplaceBack(currentNode);
		setid[currentNode->index()] = level;

		//need new level for next possible clusters
		level++;
	} else {
		//cluster already exists
		//setid equals level of clusterhead
		setid[currentNode->index()] = level - 1;
		(clusters.rbegin())->clusterNodes.emplaceBack(currentNode);
	}
#ifdef OGDF_DEBUG
	lout(Level::Minor)
			<< "      After Cluster update pi: " << clusters[level - 1].clusterhead->index()
			<< std::endl;
#endif
}

//LValues[maxAdj] for L
void MinimumCutNagamochiIbaraki::fillL(const int& maxAdj, ListPure<node>& unviewed, BoundedList& L,
		std::vector<adjInfo>& adjToViewed) {
	//repair List if it is empty
	//should have at least one element in list
	auto it = unviewed.begin();

	while (it != unviewed.end()) {
		const auto& adjValue = adjToViewed[(*it)->index()].adj;
		if (adjValue == -1) {
			auto temp = it;
			it++;
			unviewed.del(temp);
			continue;
		}
		if (adjValue == maxAdj) {
			L.nodesInList++;
			adjToViewed[(*it)->index()].placeInL = L.list.emplaceBack(*it);
		}
		it++;
	}
}

void MinimumCutNagamochiIbaraki::deleteFromL(BoundedList& L, ListIterator<node>& placeInL) {
	if (placeInL.operator!=(nullptr)) {
		L.list.del(placeInL);
		L.nodesInList--;
		placeInL = nullptr;
	}
	L.size--;
}

node MinimumCutNagamochiIbaraki::getFirstNode(BoundedList& L) {
	L.size--;
	L.nodesInList--;
	return L.list.popFrontRet();
}

node MinimumCutNagamochiIbaraki::MAOComputation(const node& s) {
	std::vector<adjInfo> adjToViewed;
	adjToViewed.resize(N);

	//lesser is better than greater for key comparer
	std::priority_queue<int> LValues;
	std::unordered_map<int, BoundedList> L;
	L.rehash(size / 2);
	node adjNodeToViewed;

	ListPure<node> unviewed;
	for (auto v : allNodes) {
		unviewed.emplaceBack(v);
	}

	node lastViewed = s;
	auto viewedSize = 1;
	adjToViewed[lastViewed->index()].adj = -1;
	auto maxAdj = 0;

	auto level = 0;

	std::vector<clusterstruct> clusters;
	clusters.reserve(OGDF_MINCUTNI_CLUSTERSIZE);

	LValues.emplace(0);
	L.emplace(std::make_pair(0, BoundedList {ListPure<node>(), 0, size - 1}));

	edge e;
	while (true) {
		//adjust adjacency values -> Iterate over all neighbours
#ifdef OGDF_DEBUG
		lout(Level::Minor) << "lastViewed: " << lastViewed->index() << std::endl;
#endif
		for (const auto& adj : lastViewed->adjEntries) {
			//neighbour node
			e = getAdjEdge(adj, lastViewed, adjNodeToViewed);

			const auto adjNode_index = adjNodeToViewed->index();
			auto& currentAdj = adjToViewed[adjNode_index].adj;
			//search end in remaining nodes
#ifdef OGDF_DEBUG
			lout(Level::Minor) << "  adjNode:" << adjNode_index << std::endl;
			lout(Level::Minor) << "    adjValue: " << currentAdj << std::endl;
#endif
			if (currentAdj != -1) {
				//Delete Node from LValues
				deleteFromL(L[currentAdj], adjToViewed[adjNode_index].placeInL);

				currentAdj += edgeCapacity[e->index()];

#ifdef OGDF_DEBUG
				lout(Level::Minor) << "    updated adjValue: " << currentAdj << std::endl;
#endif
				if (L.find(currentAdj) == L.end()) {
					// L doesnt exist yet
					LValues.emplace(currentAdj);
					L.emplace(std::make_pair(currentAdj, BoundedList {ListPure<node>(), 0, 0}));
					Math::updateMax(maxAdj, currentAdj);
				}

				auto& L_currentAdj = L[currentAdj];
				//add it to new place if  the size of L isnt above maxListSize
				L_currentAdj.size++;
				if (L_currentAdj.nodesInList < OGDF_MINCUTNI_MAXLISTSIZE) {
					adjToViewed[adjNode_index].placeInL =
							L_currentAdj.list.emplaceBack(adjNodeToViewed);

					L_currentAdj.nodesInList++;
				}

				if (setid[adjNode_index] == -1 && currentAdj >= barLambda) {
					//add node to V and edge to E_lambda,o
					updateClusters(lastViewed, adjNodeToViewed, clusters, level);
#ifdef OGDF_DEBUG
					lout(Level::Minor) << "      add node " << adjNode_index << " to V" << std::endl;
					lout(Level::Minor)
							<< "        lastViewed id " << setid[lastViewed->index()] << std::endl;
					lout(Level::Minor) << "        adjNode id " << setid[adjNode_index] << std::endl;
#endif
				}
			}
		}

#ifdef OGDF_DEBUG
		lout(Level::Minor) << "AFTER ADJ UPDATE" << std::endl;
		lout(Level::Minor) << "  currentNode: " << lastViewed->index() << std::endl;
		lout(Level::Minor) << "  maxAdj: " << maxAdj << std::endl;
		lout(Level::Minor) << "  size: " << size << " viewedSize: " << viewedSize << std::endl;
		lout(Level::Minor) << std::endl;
#endif

		if (L[maxAdj].nodesInList == 0) {
#ifdef OGDF_DEBUG
			lout() << "update" << std::endl;
			lout() << "  L " << maxAdj << std::endl;
			lout() << "  viewedSize " << viewedSize << std::endl;
#endif
			fillL(maxAdj, unviewed, L[maxAdj], adjToViewed);
		}

		//add new node to viewed and delete it from L
		lastViewed = getFirstNode(L[maxAdj]);
		viewedSize++;

		//mark lastViewed as viewed
		adjToViewed[lastViewed->index()].adj = -1;

		if (viewedSize == size) {
			//all nodes viewed
			break;
		}

		//check if L[maxAdj] is empty
		while (L[maxAdj].size == 0) {
			//remove maxAdj from L*
			L.erase(maxAdj);
			LValues.pop();
			// maxAdj is end of the List
			maxAdj = LValues.top();
		}
	}

#ifdef OGDF_DEBUG
	bool isEmpty = true;
	for (int i = level; i--;) {
		if (!clusters[i].clusterNodes.empty()) {
			isEmpty = false;
		}
	}
	OGDF_ASSERT(!isEmpty);
#endif
	contractClusters(clusters);

	return clusters.rbegin()->clusterhead;
}

const int& MinimumCutNagamochiIbaraki::minCutWeighted(const Graph& G,
		const std::vector<int>& capacity) {
	init(G);
	for (edge e : G.edges) {
		edgeCapacity[m_GC.copy(e)->index()] = capacity[e->index()];
#ifdef OGDF_DEBUG
		lout(Level::Minor) << "Orig Edge (" << e->source()->index() << "," << e->target()->index()
						   << "): " << capacity[e->index()] << std::endl;
		lout(Level::Minor) << "Edge (" << e->source()->index() << "," << e->target()->index()
						   << "): " << edgeCapacity[e->index()] << std::endl;
#endif
	}
	const int& value {minCutWeighted()};
	delete hiddenEdges;
	return value;
}

const int& MinimumCutNagamochiIbaraki::minCutWeighted() {
	if (N < 2) {
		barLambda = std::numeric_limits<int>::max();
		return barLambda;
	}

	//test if edge capacities are allowed
	List<edge> cap0edges;
	for (const auto& e : m_GC.edges) {
		//wrong capacity if less-equal 0
		if (edgeCapacity[e->index()] == 0) {
			cap0edges.pushBack(e);
		}
		OGDF_ASSERT(edgeCapacity[e->index()] >= 0);
	}
	for (edge e : cap0edges) {
		m_GC.delEdge(e);
	}
	if (!isConnected(m_GC)) {
		barLambda = 0;
		return barLambda;
	}

	//compute degree and barLambda
	fill_n(degree.begin(), N, 0);

	for (const auto& adj : m_GC.firstNode()->adjEntries) {
		barLambda += edgeCapacity[adj->theEdge()->index()];
	}
	for (const auto& v : m_GC.nodes) {
		for (auto adj : v->adjEntries) {
			degree[v->index()] += edgeCapacity[adj->theEdge()->index()];
		}

		if (degree[v->index()] < barLambda) {
			barLambda = degree[v->index()];
		}
	}

	minCut();
	return barLambda;
}

const int& MinimumCutNagamochiIbaraki::minCutUnweighted(const Graph& G) {
	init(G);
	if (N < 2) {
		barLambda = std::numeric_limits<int>::max();
		return barLambda;
	}
	if (!isConnected(m_GC)) {
		barLambda = 0;
		return barLambda;
	}

	//set edge capacities
	fill_n(edgeCapacity.begin(), M, 1);

	//compute degree
	barLambda = m_GC.firstNode()->degree();
	for (const auto& v : m_GC.nodes) {
		const auto& deg = v->degree();
		degree[v->index()] = deg;
		if (deg < barLambda) {
			barLambda = deg;
		}
	}

	minCut();
	return barLambda;
}

//NAGAMOCHI-IBARAKI Section
void MinimumCutNagamochiIbaraki::minCut() {
	fill_n(setid.begin(), N, -1);

#ifdef OGDF_DEBUG
	lout(Level::High) << "lambda " << barLambda << std::endl;
	OGDF_ASSERT(barLambda > 0);
	lout(Level::Minor) << "min cut init done" << std::endl;
	lout(Level::Minor) << "size: " << size << std::endl;
#endif

	//Preprocess
	node lastContracted = m_GC.firstNode();
	switch (barLambda) {
	case 0:
		return;
	case 1:
		return;
	default: {
		//use PR Preprocess if wished
		if (m_preprocess) {
			for (auto e : m_GC.edges) {
				const auto e_index = e->index();
				const auto s = e->source();
				const auto t = e->target();
				if (PRTest1(e_index) || PRTest2(e_index, s->index(), t->index())) {
					PRPass1_2(s);
					lastContracted = s;
					break;
				}
			}
		}
		break;
	}
	}

	// MinimumCutNagamochiIbaraki MAIN ALG
	while (size > 2) {
#ifdef OGDF_DEBUG
		lout(Level::Minor) << "new round" << std::endl;
		lout(Level::Minor) << "-----------------------------------------" << std::endl;
		for (edge e : m_GC.edges) {
			lout() << "  edge e: (" << e->source()->index() << "," << e->target()->index()
				   << ") capacity: " << edgeCapacity[e->index()] << std::endl;
		}
		lout() << "size: " << size << std::endl;
		lout() << "barLambda: " << barLambda << std::endl;
#endif
		NIRounds++;
		lastContracted = MAOComputation(lastContracted);

		//use Padberg-Rinaldi if wished
		if (pr) {
			if (size <= 2) {
				break;
			}
			PRPass1_2(lastContracted);
		}
	}
#ifdef OGDF_DEBUG
	lout(Level::Minor) << "Mincut done" << std::endl;
#endif
	return;
}

//test lastContracted until PR doesnt contract
//contract all adjacent nodes of lastContracted into lastContracted if PRPass succeeds
void MinimumCutNagamochiIbaraki::PRPass1_2(const node& lastContracted) {
	node t, opposite;
	adjEntry adj, succ, t_adj, t_succ;
	const int& lastContracted_index = lastContracted->index();
	auto& parent_degree = degree[lastContracted_index];
	edge e, t_e;

	std::unordered_map<node, edge> neigh;
	neigh.reserve(size / 4);
	for (const auto& s_adj : lastContracted->adjEntries) {
		e = getAdjEdge(s_adj, lastContracted, t);
		neigh.insert(std::make_pair(t, e));
	}
	OGDF_ASSERT(!neigh.empty());

	do {
#ifdef OGDF_DEBUG
		lout(Level::Minor) << "------PR 1 2 Round------" << std::endl;
		lout(Level::Minor) << "last contracted: " << lastContracted_index << std::endl;
		lout(Level::Minor) << "  neigh: ";
		for (auto neighnode : neigh) {
			auto e_t = neighnode.second;
			auto v_t = neighnode.first;
			lout(Level::Minor) << "    ([" << v_t->index() << "], (" << e_t->source()->index()
							   << "," << e_t->target()->index() << ") ) ";
		}
		lout(Level::Minor) << std::endl;
#endif
		auto shrunk = 0;
		auto contracted = 0;
		adj = lastContracted->firstAdj();
		//first iterate over all neighbours  and save them in neigh
		while (adj != nullptr) {
			e = getAdjEdge(adj, lastContracted, t);
			const auto& t_index = t->index();
			const auto& e_index = e->index();
			succ = adj->succ();
#ifdef OGDF_DEBUG
			lout(Level::Minor) << "edge e: (" << e->source()->index() << "," << e->target()->index()
							   << ")" << std::endl;
			if (succ != nullptr) {
				lout(Level::Minor) << "next edge e: (" << succ->theEdge()->source()->index() << ","
								   << succ->theEdge()->target()->index() << ")" << std::endl;
			}
#endif
			//PR Pass 1
			if (PRTest1(e_index)) {
#ifdef OGDF_DEBUG
				lout() << "  PR Pass 1 success" << std::endl;
#endif
				shrunk++;
			} else if (PRTest2(e_index, lastContracted_index, t_index)) {
				//PR Pass 2
#ifdef OGDF_DEBUG
				lout() << "  PR Pass 2 success" << std::endl;
#endif
				shrunk++;
			}

			//PR Test 1 or 2 succeeded
			if (shrunk > contracted) {
				neigh.erase(t);
				parent_degree -= edgeCapacity[e_index];
				contracted = shrunk;
				hiddenEdges->hide(e);

				t_adj = t->firstAdj();
#ifdef OGDF_DEBUG
				lout(Level::Minor) << "    current t: " << t->index() << std::endl;
				lout(Level::Minor) << "    node " << lastContracted_index << " with degree "
								   << parent_degree << std::endl;
#endif
				while (t_adj != nullptr) {
					t_succ = t_adj->succ();
					t_e = getAdjEdge(t_adj, t, opposite);
#ifdef OGDF_DEBUG
					lout(Level::Minor) << "      edge e: (" << t_e->source()->index() << ","
									   << t_e->target()->index() << ")" << std::endl;
					if (t_succ != nullptr) {
						lout(Level::Minor)
								<< "      next edge e: (" << t_succ->theEdge()->source()->index()
								<< "," << t_succ->theEdge()->target()->index() << ")" << std::endl;
					}
#endif

					parent_degree += edgeCapacity[t_e->index()];
					if (neigh.find(opposite) == neigh.end()) {
						//s,opposite doesnt exist
						if (t_e->source() == t) {
							m_GC.moveSource(t_e, lastContracted);
						} else {
							m_GC.moveTarget(t_e, lastContracted);
						}
#ifdef OGDF_DEBUG
						lout(Level::Minor) << "      Adding edge: (" << t_e->source()->index() << ","
										   << t_e->target()->index() << ") to neigh" << std::endl;
#endif

						neigh.insert(std::make_pair(opposite, t_e));
					} else {
						//s,opposite does exist
#ifdef OGDF_DEBUG
						lout(Level::Minor) << "      edge e: (" << t_e->source()->index() << ","
										   << t_e->target()->index() << ")" << std::endl;
#endif
						edgeCapacity[neigh[opposite]->index()] += edgeCapacity[t_e->index()];
						hiddenEdges->hide(t_e);
					}

					t_adj = t_succ;
				}
				size--;
				updateLambda(parent_degree);
#ifdef OGDF_DEBUG
				lout(Level::Minor) << "  After Contract" << std::endl;
				lout(Level::Minor) << "    node " << lastContracted_index << " with degree "
								   << parent_degree << std::endl;
				lout(Level::Minor) << "    size " << size << std::endl;
#endif
				allNodes.erase(t);
			}
			if (size == 2) {
				break;
			}
			adj = succ;
		}
		prRounds++;
		if (shrunk < OGDF_MINCUTNI_PRTHR) {
			break;
		}
	} while (size > 2);
}

}
