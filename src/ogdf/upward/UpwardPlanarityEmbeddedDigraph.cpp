/*
 * $Revision: 3504 $
 *
 * last checkin:
 *   $Author: beyer $
 *   $Date: 2013-05-16 14:49:39 +0200 (Do, 16. Mai 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of the internal class UpwardPlanarityEmbeddedDigraph.
 * 			Implements a algorithm testing upward planarity
 * 			of a given embedded Digraph G.
 *
 * \author Robert Zeranski
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.txt in the root directory of the OGDF installation for details.
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
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/

#include <ogdf/internal/upward/UpwardPlanarityEmbeddedDigraph.h>

namespace ogdf {

UpwardPlanarityEmbeddedDigraph::UpwardPlanarityEmbeddedDigraph(const Graph &H):
	//initialising the variables
	m_G(H),
	m_s(m_G.firstNode()),
	m_t(m_G.firstNode()),
	m_combEmb(m_G),
	m_A(m_combEmb,0),
	m_assignedSourcesAndSinks(m_combEmb),
	m_correspondingSourceOrSink(m_B,NULL),
	m_correspondingFace(m_B,NULL),
	m_correspondingFaceNode(m_combEmb),
	m_correspondingEdge(m_B,NULL)
{
}

//DFS for calculating a feasible augmentation path
void UpwardPlanarityEmbeddedDigraph::getPath(Stack<node> &st, EdgeArray<int> &capacity, EdgeArray<int> &flow) {
	node u = m_s;
	NodeArray<bool> dfi(m_B,false);
	dfi[u] = true;
	st.push(u);
	while (!st.empty() && u!=m_t) {
		u=st.pop();
		adjEntry adj;
		forall_adj(adj,u) {
			if (adj->theEdge()->target() == u) continue;
			node x = adj->twinNode();
			edge e = adj->theEdge();
			if (dfi[x] || flow[e] >= capacity[e]) continue;
			dfi[x] = true;
			st.push(u);
			st.push(x);
			u = x;
			break;
		}
	}
}

//finds the minimum possible increasing flow of the augmentation path
int UpwardPlanarityEmbeddedDigraph::getMin(Stack<node> stack, EdgeArray<int> &capacity, EdgeArray<int> &flow) {
	int min = -1;
	while (!stack.empty()) {
		node u = stack.pop();
		if (stack.empty()) break;
		node v = stack.top();
		adjEntry adj_u;
		forall_adj(adj_u, v) {
			if (adj_u->theEdge()->target() == u) break;
		}
		edge k = adj_u->theEdge();
		if (capacity[k]- flow[k] < min || min == -1) min = capacity[k]- flow[k];
	}
	return min;
}

//tests whether a flow of power r is feasible
bool UpwardPlanarityEmbeddedDigraph::isFlow(EdgeArray<int> &capacity, EdgeArray<int> &flow, const int r) {
	if (r == 0) return true;
	int currentFlow = 0;
	bool rFlow = false;
	edge e;
	EdgeArray<edge> R(m_B);
	adjEntry adj;
	forall_edges(e,m_B) {
		bool check = true;
		node u = e->source();
		node v = e->target();
		forall_adj(adj,v) {
			if (adj->theEdge()->target() == u) {
				check = false;
				edge k = adj->theEdge();
				R[e] = k;
				R[k] = e;
			}
		}
		if (check) {
			//insert backedges
			edge f = m_B.newEdge(v,u);
			capacity[f] = 0;
			R[e] = f;
			R[f] = e;
			flow[f] = 0;
		}
	}
	Stack<node> stack;
	while (!rFlow) {
		//find the augmentation path
		getPath(stack,capacity,flow);
		//calculate the value for one augmentation step
		int min = getMin(stack,capacity,flow);
		if (stack.empty()) break;
		while (!stack.empty()) {
				//increase the flow for the augmentation path
				node u = stack.pop();
				if (stack.empty()) break;
				node v = stack.top();
				adjEntry adj_u;
				forall_adj(adj_u, v) {
					if (adj_u->theEdge()->target() == u) break;
				}
				edge k = adj_u->theEdge();
				flow[k] = flow[k] + min;
				flow[R[k]] = -flow[k];
		}
		currentFlow = 0;
		forall_adj(adj,m_s) {
			//calculate the current flow in B
			if (adj->theEdge()->target() == m_s) continue;
			edge k = adj->theEdge();
			currentFlow = currentFlow + flow[k];
		}
		if (currentFlow >= r) rFlow = true;
	}
	return rFlow;
}

//constructs a flow-network corresponding to graph m_G
void UpwardPlanarityEmbeddedDigraph::constructNetwork(EdgeArray<int> &capacity, EdgeArray<int> &flow) {
	node v,w;
	edge e;
	face f;
	//super-source
	node s = m_B.newNode();
	//super-sink
	node t = m_B.newNode();
	forall_nodes(v,m_G) {
		//assign for all sources and sinks of m_G a corresponding node in the flow-network m_B
		if ((v->indeg() == 0) || (v->outdeg() == 0)) {
			w = m_B.newNode();
			m_correspondingSourceOrSink[w] = v;
			e = m_B.newEdge(s,w);
			capacity[e] = 1;
		}
	}
	forall_faces(f,m_combEmb) {
		//assign for all faces of m_G a corresponding node in the flow-network m_B
		w = m_B.newNode();
		m_correspondingFace[w] = f;
		m_correspondingFaceNode[f] = w;
		e = m_B.newEdge(w,t);
		m_correspondingEdge[w] = e;
		capacity[e] = m_A[f]-1;
	}
	forall_nodes(v,m_B) {
		//insert edges between nodes in m_B if the corresponding nodes in m_G are a source or sink of the corresponding face in m_G
		if (m_correspondingSourceOrSink[v] != NULL) {
			forall_nodes(w,m_B) {
				if (m_correspondingFace[w] != NULL) {
					f = m_correspondingFace[w];
					ListIterator<node> it = m_assignedSourcesAndSinks[f].begin();
					while (it.valid()) {
						if (*it == m_correspondingSourceOrSink[v]) {
							e = m_B.newEdge(v,w);
							capacity[e] = 1;
						}
						it++;
					}
				}
			}
		}
	}
}

//tests whether G is upward-planar (fixed embedding)
void UpwardPlanarityEmbeddedDigraph::isUpwardPlanarEmbedded(bool val, List<adjEntry> &possibleExternalFaces) {
		node v,w;
		face f;
		EdgeArray<int> capacity(m_B,0);
		EdgeArray<int> flow(m_B,0);
		//stack of feasible external faces in an upward-planar drawing
		adjEntry adj;
		List<edge> list;
		edge e;
		forall_faces(f,m_combEmb) {
			//compute the number m_A of angles in face f and the corresponding sources and sinks of f
			adj = f->firstAdj();
			forall_face_adj(adj,f) {
				e = adj->theEdge();
				if (!list.empty()) {
					if (e->target() == list.back()->target()) {
						m_A[f]++;
						w = e->target();
						m_assignedSourcesAndSinks[f].pushBack(w);
					}
					else if (e->source() == list.back()->source()) {
						m_A[f]++;
						w = e->source();
						m_assignedSourcesAndSinks[f].pushBack(w);
					}
				}
				list.pushBack(e);
			}
			if ((list.front()->target() == list.back()->target()) )  {
				m_A[f]++;
				w = list.front()->target();
				m_assignedSourcesAndSinks[f].pushBack(w);
			}
			else if (list.front()->source() == list.back()->source()) {
				m_A[f]++;
				w = list.front()->source();
				m_assignedSourcesAndSinks[f].pushBack(w);
			}
			m_A[f] = m_A[f]/2;
			list.clear();
		}
		//construct flow-network m_B corresponding to m_G
		constructNetwork(capacity,flow);
		forall_nodes(v,m_B) {
			if (v->index() == 0) m_s = v;
			if (v->index() == 1) m_t = v;
			if (v->index() > 1) break;
		}
		int r = 0;
		forall_nodes(v,m_G) {
			//r = number of sources and sinks in G
			if (v->indeg() == 0 || v->outdeg() == 0) r++;
		}
		//tests whether the network permits a flow of power r-2 without determining the external face
		if (isFlow(capacity,flow,r-2)) {
			forall_faces(f,m_combEmb) {
				//tests whether the network permits a flow of power r for the choice of face f as external face
				v = m_correspondingFaceNode[f];
				e = m_correspondingEdge[v];
				capacity[e] = m_A[f] + 1;
				EdgeArray<int> capacityCopy(capacity);
				EdgeArray<int> flowCopy(flow);
				if (isFlow(capacityCopy,flowCopy,r)) {
					possibleExternalFaces.pushBack(f->firstAdj());
					if (val) break;
				}
				capacity[e] = m_A[f]-1;
			}
		}
}

//tests whether G is upward-planar (fixed embedding)
//returns the set of feasible external faces (represented by the first AdjEntry)
bool UpwardPlanarityEmbeddedDigraph::isUpwardPlanarEmbedded() {
	List<adjEntry> possibleExternalFaces;
	isUpwardPlanarEmbedded(true, possibleExternalFaces);
	return !possibleExternalFaces.empty();
}

//tests whether G is upward-planar (fixed embedding)
bool UpwardPlanarityEmbeddedDigraph::isUpwardPlanarEmbedded(List<adjEntry> &possibleExternalFaces) {
	isUpwardPlanarEmbedded(false, possibleExternalFaces);
	return !possibleExternalFaces.empty();
}


}
