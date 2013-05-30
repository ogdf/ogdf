/*
 * $Revision: 3504 $
 *
 * last checkin:
 *   $Author: beyer $
 *   $Date: 2013-05-16 14:49:39 +0200 (Do, 16. Mai 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of class PlanarSPQRTree
 *
 * \author Carsten Gutwenger
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


#include <ogdf/decomposition/PlanarSPQRTree.h>
#include <ogdf/basic/extended_graph_alg.h>


namespace ogdf {

//-------------------------------------------------------------------
//                           PlanarSPQRTree
//-------------------------------------------------------------------

//
// initialization: additionally embeds skeleton graphs or adpots embedding
// given by original graph
void PlanarSPQRTree::init(bool isEmbedded)
{
	m_finished = true;
	if (isEmbedded) {
		adoptEmbedding();

	} else {

		node v;
		forall_nodes(v,tree())
			planarEmbed(skeleton(v).getGraph());
	}
}


void PlanarSPQRTree::adoptEmbedding()
{
	OGDF_ASSERT_IF(dlExtendedChecking, originalGraph().representsCombEmbedding());

	// ordered list of adjacency entries (for one original node) in all
	// skeletons (where this node occurs)
	NodeArray<SListPure<adjEntry> > adjEdges(tree());
	// copy in skeleton of current original node
	NodeArray<node> currentCopy(tree(),0);
	NodeArray<adjEntry> lastAdj(tree(),0);
	SListPure<node> current; // currently processed nodes

	node vOrig;
	forall_nodes(vOrig,originalGraph())
	{
		adjEntry adjOrig;
		forall_adj(adjOrig,vOrig)
		{
			edge            eOrig = adjOrig->theEdge();
			const Skeleton &S     = skeletonOfReal(eOrig);
			edge            eCopy = copyOfReal(eOrig);

			adjEntry adjCopy = (S.original(eCopy->source()) == vOrig) ?
				eCopy->adjSource() : eCopy->adjTarget();

			setPosInEmbedding(adjEdges,currentCopy,lastAdj,current,S,adjCopy);
		}

		SListConstIterator<node> it;
		for(it = current.begin(); it.valid(); ++it) {
			node vT = *it;

			skeleton(vT).getGraph().sort(currentCopy[vT],adjEdges[vT]);

			adjEdges[vT].clear();
			currentCopy[vT] = 0;
		}

		current.clear();
	}
}


void PlanarSPQRTree::setPosInEmbedding(
	NodeArray<SListPure<adjEntry> > &adjEdges,
	NodeArray<node> &currentCopy,
	NodeArray<adjEntry> &lastAdj,
	SListPure<node> &current,
	const Skeleton &S,
	adjEntry adj)
{
	node vT = S.treeNode();

	adjEdges[vT].pushBack(adj);

	node vCopy = adj->theNode();
	node vOrig = S.original(vCopy);

	if(currentCopy[vT] == 0) {
		currentCopy[vT] = vCopy;
		current.pushBack(vT);

		adjEntry adjVirt;
		forall_adj(adjVirt,vCopy) {
			edge eCopy = S.twinEdge(adjVirt->theEdge());
			if (eCopy == 0) continue;
			if (adjVirt == adj) {
				lastAdj[vT] = adj;
				continue;
			}

			const Skeleton &STwin = skeleton(S.twinTreeNode(adjVirt->theEdge()));

			adjEntry adjCopy = (STwin.original(eCopy->source()) == vOrig) ?
				eCopy->adjSource() : eCopy->adjTarget();

			setPosInEmbedding(adjEdges,currentCopy,lastAdj,current,
				STwin, adjCopy);
		}

	} else if (lastAdj[vT] != 0 && lastAdj[vT] != adj) {
		adjEntry adjVirt = lastAdj[vT];
		edge eCopy = S.twinEdge(adjVirt->theEdge());

		const Skeleton &STwin = skeleton(S.twinTreeNode(adjVirt->theEdge()));

		adjEntry adjCopy = (STwin.original(eCopy->source()) == vOrig) ?
			eCopy->adjSource() : eCopy->adjTarget();

		setPosInEmbedding(adjEdges,currentCopy,lastAdj,current,
			STwin, adjCopy);

		lastAdj[vT] = 0;
	}

}


//
// embed original graph according to embedding of skeletons
//
// The procedure also handles the case when some (real or virtual)
// edges are reversed (used in upward-planarity algorithms)
void PlanarSPQRTree::embed(Graph &G)
{
	OGDF_ASSERT(&G == &originalGraph());

	const Skeleton &S = skeleton(rootNode());
	const Graph &M = S.getGraph();

	node v;
	forall_nodes(v,M)
	{
		node vOrig = S.original(v);
		SListPure<adjEntry> adjEdges;

		adjEntry adj;
		forall_adj(adj,v) {
			edge e = adj->theEdge();
			edge eOrig = S.realEdge(e);

			if (eOrig != 0) {
				adjEntry adjOrig = (vOrig == eOrig->source()) ?
					eOrig->adjSource() : eOrig->adjTarget();
				OGDF_ASSERT(adjOrig->theNode() == S.original(v));
				adjEdges.pushBack(adjOrig);

			} else {
				node wT    = S.twinTreeNode(e);
				edge eTwin = S.twinEdge(e);
				expandVirtualEmbed(wT,
					(vOrig == skeleton(wT).original(eTwin->source())) ?
					eTwin->adjSource() : eTwin->adjTarget(),
					adjEdges);
			}
		}

		G.sort(vOrig,adjEdges);
	}

	edge e;
	forall_adj_edges(e,rootNode()) {
		node wT = e->target();
		if (wT != rootNode())
			createInnerVerticesEmbed(G, wT);
	}
}


void PlanarSPQRTree::expandVirtualEmbed(node vT,
	adjEntry adjVirt,
	SListPure<adjEntry> &adjEdges)
{
	const Skeleton &S = skeleton(vT);

	node v = adjVirt->theNode();
	node vOrig = S.original(v);

	adjEntry adj;
	for (adj = adjVirt->cyclicSucc(); adj != adjVirt; adj = adj->cyclicSucc())
	{
		edge e = adj->theEdge();
		edge eOrig = S.realEdge(e);

		if (eOrig != 0) {
			adjEntry adjOrig = (vOrig == eOrig->source()) ?
				eOrig->adjSource() : eOrig->adjTarget();
			OGDF_ASSERT(adjOrig->theNode() == S.original(v));
			adjEdges.pushBack(adjOrig);

		} else {
			node wT    = S.twinTreeNode(e);
			edge eTwin = S.twinEdge(e);
			expandVirtualEmbed(wT,
				(vOrig == skeleton(wT).original(eTwin->source())) ?
				eTwin->adjSource() : eTwin->adjTarget(),
				adjEdges);
		}
	}
}


void PlanarSPQRTree::createInnerVerticesEmbed(Graph &G, node vT)
{
	const Skeleton &S = skeleton(vT);
	const Graph& M = S.getGraph();

	node src = S.referenceEdge()->source();
	node tgt = S.referenceEdge()->target();

	node v;
	forall_nodes(v,M)
	{
		if (v == src || v == tgt) continue;

		node vOrig = S.original(v);
		SListPure<adjEntry> adjEdges;

		adjEntry adj;
		forall_adj(adj,v) {
			edge e = adj->theEdge();
			edge eOrig = S.realEdge(e);

			if (eOrig != 0) {
				adjEntry adjOrig = (vOrig == eOrig->source()) ?
					eOrig->adjSource() : eOrig->adjTarget();
				OGDF_ASSERT(adjOrig->theNode() == S.original(v));
				adjEdges.pushBack(adjOrig);
			} else {
				node wT    = S.twinTreeNode(e);
				edge eTwin = S.twinEdge(e);
				expandVirtualEmbed(wT,
					(vOrig == skeleton(wT).original(eTwin->source())) ?
					eTwin->adjSource() : eTwin->adjTarget(),
					adjEdges);
			}
		}

		G.sort(vOrig,adjEdges);
	}

	edge e;
	forall_adj_edges(e,vT) {
		node wT = e->target();
		if (wT != vT)
			createInnerVerticesEmbed(G, wT);
	}
}



//
// basic update functions for manipulating embeddings

//   reversing the skeleton of an R- or P-node
void PlanarSPQRTree::reverse(node vT)
{
	skeleton(vT).getGraph().reverseAdjEdges();
}


//   swapping two adjacency entries in the skeleton of a P-node
void PlanarSPQRTree::swap(node vT, adjEntry adj1, adjEntry adj2)
{
	OGDF_ASSERT(typeOf(vT) == PNode);

	Graph &M = skeleton(vT).getGraph();

	M.swapAdjEdges(adj1,adj2);
	M.swapAdjEdges(adj1->twin(),adj2->twin());
}


//   swapping two edges in the skeleton of a P-node
void PlanarSPQRTree::swap(node vT, edge e1, edge e2)
{
	OGDF_ASSERT(typeOf(vT) == PNode);

	if (e1->source() == e2->source())
		swap(vT,e1->adjSource(),e2->adjSource());
	else
		swap(vT,e1->adjSource(),e2->adjTarget());
}


//
// number of possible embeddings of original graph
//
double PlanarSPQRTree::numberOfEmbeddings(node vT) const
{
	double num = 1.0;

	switch(typeOf(vT)) {
	case RNode:
		num = 2; break;
	case PNode:
		//node vFirst = skeleton(vT).getGraph().firstNode();
		for (int i = skeleton(vT).getGraph().firstNode()->degree()-1; i >= 2; --i)
			num *= i;
		break;
	case SNode:
		break;
	}

	edge e;
	forall_adj_edges(e,vT) {
		node wT = e->target();
		if(wT != vT)
			num *= numberOfEmbeddings(wT);
	}

	return num;
}



//
// randomly embed skeleton graphs
//
void PlanarSPQRTree::randomEmbed()
{
	node vT;
	forall_nodes(vT,tree()) {
		if (typeOf(vT) == RNode) {
			int doReverse = randomNumber(0,1);

			if (doReverse == 1)
				reverse(vT);

		} else if (typeOf(vT) == PNode) {
			const Skeleton &S = skeleton(vT);
			adjEntry adjRef = S.referenceEdge()->adjSource();

			SList<adjEntry> adjEdges;
			adjEntry adj;
			for (adj = adjRef->cyclicSucc(); adj != adjRef; adj = adj->cyclicSucc())
				adjEdges.pushBack(adj);

			adjEdges.permute();

			adj = adjRef->cyclicSucc();
			SListConstIterator<adjEntry> it;
			for (it = adjEdges.begin(); it.valid(); ++it)
			{
				adjEntry adjNext = *it;
				if (adjNext != adj) {
					swap(vT,adj,adjNext);
					adj = adjNext;
				}
				adj = adj->cyclicSucc();
			}
		}
	}
}

//***************************************************************
//***************************************************************
//**														   **
//** Methods to enumerate all embeddings of the original graph **
//**														   **
//***************************************************************
//***************************************************************

double PlanarSPQRTree::numberOfNodeEmbeddings(node vT) {

	double num = 1.0;

	switch(typeOf(vT)) {
		case RNode:
			num = 2;
			break;
		case PNode:
			for (int i = skeleton(vT).getGraph().firstNode()->degree()-1; i >= 2; --i)
				num *= i;
			break;
		case SNode:
			break;
	}

	return num;
}

void PlanarSPQRTree::embed(node &vT, int x) {

	OGDF_ASSERT(x >= 0);
	OGDF_ASSERT(x < numberOfNodeEmbeddings(vT));

	//if it is a P-node
	if (typeOf(vT) == PNode) {
		//encode the id of the permutation
		int id = x;
		//number of elements of the permutation
		const int p = skeleton(vT).getGraph().numberOfEdges() - 1;
		//sequence for the edges
		Array<int> seq(p);
		//base
		int b = 2;
		//compute the correct sequence of the edges
		for (int i = 0; i < p - 1; i++) {
			seq[i] = id % b;
			id = id / b;
			b++;
		}
		//swap the sequence
		int low  = 0;
		int high = p-2;
		while (low < high) {
			int t = seq[low];
			seq[low] = seq[high]; seq[high] = t;
			low++; high--;
		}
		seq[p-1] = 0;
		//encode the sequence to the permutation
		Array<int> list(p);
		Array<int> permutation(p);
		Array<bool> set(0,p-1,false);
		for (int i = 0; i < p; i++) list[i] = i;
		for (int i = 0; i < p; i++) {
			int e = seq[i];
			int cnt = 0;
			int ind;
			for (ind = 0; ind < p; ind++) {
					if (!set[ind]) {
							if (cnt == e) break;
							cnt++;
					}
			}
			permutation[ind] = list[i];
			set[ind] = true;
		}
		//adopt the permutation to the skeleton
		//first vertex
		List<adjEntry> order;
		node nP = skeleton(vT).getGraph().firstNode();
		skeleton(vT).getGraph().adjEntries(nP,order);
		TargetComparer<AdjElement,AdjElement> comp;
		order.quicksort(comp);
		Array<adjEntry> normalized(p);
		List<adjEntry> newOrder;
		newOrder.pushBack(order.popFrontRet());
		for (int i = 0; i < p; i++) normalized[i] = order.popFrontRet();
		for (int i = 0; i < p; i++) newOrder.pushBack(normalized[permutation[i]]);
		skeleton(vT).getGraph().sort(nP,newOrder);
		//second vertex
		List<adjEntry> newOrderLast;
		ListIterator<adjEntry> it = newOrder.begin();
		while (it.valid()) {
			newOrderLast.pushFront((*it)->twin());
			it++;
		}
		skeleton(vT).getGraph().sort(skeleton(vT).getGraph().lastNode(),newOrderLast);
		//P-node ends here
	}
	//if it is a R-node
	if (typeOf(vT) == RNode) {
		node nP = skeleton(vT).getGraph().firstNode();
		if (x == 0 && nP->firstAdj()->index() > nP->lastAdj()->index()) reverse(vT);
		if (x == 1 && nP->firstAdj()->index() < nP->lastAdj()->index()) reverse(vT);
		//R-node ends here
	}

}

//
// compute the first embedding of the original graph
//
void PlanarSPQRTree::firstEmbedding(Graph &G)
{

	OGDF_ASSERT(&G == &originalGraph());

	m_finished = false;
	node vT;
	forall_nodes(vT,tree()) {
		firstEmbedding(vT);
	}
	embed(G);
}


//
// compute the next embedding of the original graph
//
bool PlanarSPQRTree::nextEmbedding(Graph &G)
{

	OGDF_ASSERT(&G == &originalGraph());

	//if there is at least one new embedding: compute it using
	//nextEmbedding(n) to the first node of the SPQR-tree (represented by a list-iterator)
	List<node> nodes;
	tree().allNodes(nodes);
	if(!m_finished && nextEmbedding(nodes.begin())) {
		//compute the actual embedding of the original graph
		embed(G);
		return true;
	}
	m_finished = true;
	return false;
}

//
// computes the first embedding of the skeleton of node vT
//
void PlanarSPQRTree::firstEmbedding(node &vT)
{
	//if vT is a R-node
	if (typeOf(vT) == RNode) {
		//if the R-node were reversed in former steps
		//then reverse it to its original embedding
		node nP = skeleton(vT).getGraph().firstNode();
		if (nP->firstAdj()->index() > nP->lastAdj()->index()) reverse(vT);
	}
	//if vT is a P-node
	if (typeOf(vT) == PNode) {
		//then sort the adjEntries by their indices
		//first vertex
		List<adjEntry> order;
		node nP = skeleton(vT).getGraph().firstNode();
		skeleton(vT).getGraph().adjEntries(nP,order);
		TargetComparer<AdjElement,AdjElement> comp;
		order.quicksort(comp);
		skeleton(vT).getGraph().sort(nP,order);
		//second vertex
		List<adjEntry> newOrderLast;
		ListIterator<adjEntry> it = order.begin();
		while (it.valid()) {
			newOrderLast.pushFront((*it)->twin());
			it++;
		}
		skeleton(vT).getGraph().sort(skeleton(vT).getGraph().lastNode(),newOrderLast);
	}

}

void PlanarSPQRTree::reverse(node &vP, adjEntry &first, adjEntry &last)
{
	//swap the first and the last adjEntry
	adjEntry it_f = first;
	adjEntry it_l = last;
	swap(vP,it_f,it_l);
	adjEntry temp = it_f;
	it_f = it_l->succ();
	it_l = temp->pred();
	//while there are swapable adjEntries: swap it, i.e.
	//if left == right or left->pred() == right->succ() then stop
	while (it_f != it_l && it_l->succ() != it_f) {
			swap(vP,it_f,it_l);
			temp = it_f;
			it_f = it_l->succ();
			it_l = temp->pred();
	}
}


//
// computes the next embedding of the skeleton of node vT
//
bool PlanarSPQRTree::nextEmbedding(node &vT)
{
	//if vT is a R-node
	if (typeOf(vT) == RNode) {
		node nP = skeleton(vT).getGraph().firstNode();
		//compute the next embedding (might be the first embedding)
		//of the sekeleton of vT
		reverse(vT);
		//if next embedding = first embedding then return false
		//otherwise return true
		return nP->firstAdj()->index() > nP->lastAdj()->index();
	}
	//if vT is a P-node
	if (typeOf(vT) == PNode) {
		//take on of the two nodes of the skeleton of vT
		node nP = skeleton(vT).getGraph().firstNode();
		//if its degree is two then there is only one embedding
		if (nP->degree() < 3) return false;
		//otherwise compute the next embedding by computing the next
		//permutation of the adjEntries of nP -- excluding the first adjEntry
		adjEntry last;
		adjEntry it = nP->lastAdj();
		//compute the largest adjEntry s.t. its index is smaller than its successor
		//largest means that the adjEntry is near the last adjEntry
		while (it->index() < it->pred()->index()) it = it->pred();
		//if such adjEntry was not found then all permutations were computed -- thus
		//compute the first embedding by reversing the order of the adjEntries
		//beginning at the second adjEntry
		if (it == nP->firstAdj()->succ()) {
				last = nP->lastAdj();
				reverse(vT,it,last);
				return false;
		}
		//otherwise compute the next permutation (embedding)
		it = it->pred();
		adjEntry it_max = nP->lastAdj();
		//compute the largest adjEntry s.t. its index is smaller than the index of
		//the adjEntry computed above -- the existence is guaranteed
		while (it->index() > it_max->index()) it_max = it_max->pred();
		//swap both adjEntries found above
		swap(vT,it,it_max);
		//reverse the order of the AdjEntries between the found adjEnrty
		//and the last adjEntry of nP
		last = nP->lastAdj();
		it_max = it_max->succ();
		if (it_max != NULL && it_max != last) reverse(vT,it_max,last);
		//the next permutation (embedding) was computed
		return true;
	}
	//if vT is a S-node there is only one embedding
	return false;
}

//
// computes the next embedding of the skeleton of node *it
//
bool PlanarSPQRTree::nextEmbedding(ListIterator<node> it)
{
	//if the last embedding of the actual SPQR-node *it was computed: compute
	//the first embedding of *it and the next embedding of *it++
	//otherwise: compute the next embedding of *it
	if (!nextEmbedding(*it)) {
		it++;
		//if there is a valid successor of *it in the list of P- and R-nodes: compute
		//its next embedding
		//otherwise: all embeddings are computed
		if (it.valid()) return nextEmbedding(it);
		else return false;
	}
	return true;
}


} // end namespace ogdf
