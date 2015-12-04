/** \file
 * \brief Implementation of GraphCopySimple and GraphCopy classes
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


#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/FaceSet.h>
#include <ogdf/basic/extended_graph_alg.h>


namespace ogdf {


//---------------------------------------------------------
// GraphCopySimple
// simple graph copies (no support for edge splitting)
//---------------------------------------------------------

GraphCopySimple::GraphCopySimple(const Graph &G)
{
	init(G);
}

GraphCopySimple::GraphCopySimple(const GraphCopySimple &GC) : Graph()
{
	*this = GC;
}

void GraphCopySimple::init(const Graph &G)
{
	m_pGraph = &G;

	Graph::construct(G,m_vCopy,m_eCopy);

	m_vOrig.init(*this,nullptr);
	m_eOrig.init(*this,nullptr);

	for(node v : G.nodes)
		m_vOrig[m_vCopy[v]] = v;

	for(edge e : G.edges)
		m_eOrig[m_eCopy[e]] = e;
}


GraphCopySimple &GraphCopySimple::operator=(const GraphCopySimple &GC)
{
	NodeArray<node> vCopy;
	EdgeArray<edge> eCopy;

	Graph::assign(GC,vCopy,eCopy);
	initGC(GC,vCopy,eCopy);

	return *this;
}


void GraphCopySimple::initGC(const GraphCopySimple &GC,
	NodeArray<node> &vCopy,
	EdgeArray<edge> &eCopy)
{
	m_pGraph = GC.m_pGraph;

	m_vOrig.init(*this,nullptr); m_eOrig.init(*this,nullptr);
	m_vCopy.init(*m_pGraph,nullptr); m_eCopy.init(*m_pGraph,nullptr);

	for(node v : GC.nodes) {
		node w = GC.m_vOrig[v];
		m_vOrig[vCopy[v]] = w;
		if(w != nullptr) {
			m_vCopy[w] = vCopy[v];
		}
	}

	for(edge e : GC.edges) {
		edge eOrig = GC.m_eOrig[e];
		m_eOrig[eCopy[e]] = eOrig;
		if (eOrig)
			m_eCopy[eOrig] = eCopy[e];
	}
}


//---------------------------------------------------------
// GraphCopy
// graph copies (support for edge splitting)
//---------------------------------------------------------

GraphCopy::GraphCopy(const Graph &G)
{
	init(G);
}


GraphCopy::GraphCopy(const GraphCopy &GC) : Graph()
{
	*this = GC;
}


void GraphCopy::initGC(const GraphCopy &GC,
	NodeArray<node> &vCopy,
	EdgeArray<edge> &eCopy)
{
	createEmpty(*GC.m_pGraph);

	for(node v : GC.nodes)
		m_vOrig[vCopy[v]] = GC.original(v);

	for(edge e : GC.edges)
		m_eOrig[eCopy[e]] = GC.original(e);

	for (node v : nodes) {
		node w = m_vOrig[v];
		if (w != nullptr)
			m_vCopy[w] = v;
	}

	for(edge e : m_pGraph->edges) {
		ListConstIterator<edge> it;
		for (edge ei : GC.m_eCopy[e])
			m_eIterator[eCopy[ei]] = m_eCopy[e].pushBack(eCopy[ei]);
	}
}


void GraphCopy::init(const Graph &G)
{
	m_pGraph = &G;

	EdgeArray<edge> eCopy;
	Graph::construct(G,m_vCopy,eCopy);

	m_vOrig.init(*this,nullptr);
	m_eOrig.init(*this,nullptr);
	m_eCopy.init(G);
	m_eIterator.init(*this,nullptr);

	for(node v : G.nodes)
		m_vOrig[m_vCopy[v]] = v;

	for(edge e : G.edges) {
		m_eIterator[eCopy[e]] = m_eCopy[e].pushBack(eCopy[e]);
		m_eOrig[eCopy[e]] = e;
	}
}


void GraphCopy::createEmpty(const Graph &G)
{
	m_pGraph = &G;

	m_vCopy.init(G,nullptr);
	m_eCopy.init(G);
	m_vOrig.init(*this,nullptr);
	m_eOrig.init(*this,nullptr);
	m_eIterator.init(*this,nullptr);
}


void GraphCopy::initByCC(const CCsInfo &info, int cc, EdgeArray<edge> &eCopy)
{
	Graph::constructInitByCC(info, cc, m_vCopy, eCopy);

	for(int i = info.startNode(cc); i < info.stopNode(cc); ++i) {
		node v = info.v(i);
		m_vOrig[m_vCopy[v]] = v;
	}

	for(int i = info.startEdge(cc); i < info.stopEdge(cc); ++i) {
		edge e = info.e(i);
		m_eIterator[eCopy[e]] = m_eCopy[e].pushBack(eCopy[e]);
		m_eOrig[eCopy[e]] = e;
	}
}


void GraphCopy::initByNodes(const List<node> &origNodes, EdgeArray<edge> &eCopy)
{
#ifdef OGDF_DEBUG
	Stack<node> *stack = new Stack<node>();
	List<node> *copyList = new List<node>(origNodes);
	while(!copyList->empty()){
		stack->push(copyList->popFrontRet());
		node v;
		while(!stack->empty()){
			v = stack->pop();
			for(adjEntry adj = v->firstAdj(); adj != nullptr; adj = adj->succ()){
				node w = adj->twinNode();
				if(copyList->search(w).valid()){
					stack->push(w);
					copyList->del(copyList->search(w));
				}
				if(!origNodes.search(w).valid()){
					delete stack;
					delete copyList;
					throw PreconditionViolatedException();
				}
			}
		}
	}
	delete stack;
	delete copyList;
#endif

	Graph::constructInitByNodes(*m_pGraph,origNodes,m_vCopy,eCopy);

	for(node v : origNodes)
	{
		m_vOrig[m_vCopy[v]] = v;

		for(adjEntry adj : v->adjEntries) {
			if ((adj->index() & 1) == 0) {
				edge e = adj->theEdge();
				//
				// edge ec = eCopy[e];
				//
				m_eIterator[eCopy[e]] = m_eCopy[e].pushBack(eCopy[e]);
				m_eOrig[eCopy[e]] = e;
			}
		}
	}
}


void GraphCopy::initByActiveNodes(
	const List<node> &nodeList,
	const NodeArray<bool> &activeNodes,
	EdgeArray<edge> &eCopy)
{
	Graph::constructInitByActiveNodes(nodeList, activeNodes, m_vCopy, eCopy);

	for(node v : nodeList)
	{
		m_vOrig[m_vCopy[v]] = v;

		for(adjEntry adj : v->adjEntries) {
			if ((adj->index() & 1) == 0) {
				edge e = adj->theEdge();
				//
				// edge ec = eCopy[e];
				//
				OGDF_ASSERT(m_eCopy[e].size() == 0)
				if (activeNodes[e->opposite(v)])
				{
					m_eIterator[eCopy[e]] = m_eCopy[e].pushBack(eCopy[e]);
					m_eOrig[eCopy[e]] = e;
				}
			}
		}
	}

}


GraphCopy &GraphCopy::operator=(const GraphCopy &GC)
{
	NodeArray<node> vCopy;
	EdgeArray<edge> eCopy;

	Graph::assign(GC,vCopy,eCopy);
	initGC(GC,vCopy,eCopy);

	return *this;
}


void GraphCopy::setOriginalEmbedding()
{
	if(m_pGraph->numberOfNodes() != numberOfNodes() || m_pGraph->numberOfEdges() != numberOfEdges()){
		throw PreconditionViolatedException();
	}
	for(node v : m_pGraph->nodes)
	{
		if (m_vCopy[v] != nullptr)
		{
			List<adjEntry> newAdjOrder;
			newAdjOrder.clear();

			for(adjEntry adjOr : v->adjEntries)
			{
				if (m_eCopy[adjOr->theEdge()].size() > 0){
					//we have outgoing adjEntries for all
					//incoming and outgoing edges, check the direction
					//to find the correct copy adjEntry
					bool outEdge = (adjOr == (adjOr->theEdge()->adjSource()));

					if(chain(adjOr->theEdge()).size() != 1){
						throw PreconditionViolatedException();
					}
					edge cEdge = chain(adjOr->theEdge()).front();
					adjEntry cAdj = (outEdge ? cEdge->adjSource() : cEdge->adjTarget());
					newAdjOrder.pushBack(cAdj);
				}
				else
				{
					throw PreconditionViolatedException();
				}
			}
			sort(copy(v), newAdjOrder);
		}
		else
		{
			throw PreconditionViolatedException();
		}
	}
}


edge GraphCopy::split(edge e)
{
	edge eNew  = Graph::split(e);
	edge eOrig = m_eOrig[e];

	if ((m_eOrig[eNew] = eOrig) != nullptr) {
		m_eIterator[eNew] = m_eCopy[eOrig].insert(eNew,m_eIterator[e],after);
	}

	return eNew;
}


void GraphCopy::unsplit(edge eIn, edge eOut)
{
	edge eOrig = m_eOrig[eOut];

	// update chain of eOrig if eOrig exists
	if (eOrig != nullptr) {
		m_eCopy[eOrig].del(m_eIterator[eOut]);
	}

	Graph::unsplit(eIn,eOut);
}


edge GraphCopy::newEdge(edge eOrig)
{
#ifdef OGDF_DEBUG
	if(eOrig == nullptr){
		throw PreconditionViolatedException();
	}
	if(eOrig->graphOf() != m_pGraph){
		throw PreconditionViolatedException();
	}
	if(!m_eCopy[eOrig].empty()){
		throw PreconditionViolatedException();
	} // no support for edge splitting!
#endif

	edge e = Graph::newEdge(m_vCopy[eOrig->source()], m_vCopy[eOrig->target()]);
	m_eCopy[m_eOrig[e] = eOrig].pushBack(e);

	return e;
}


//inserts edge preserving the embedding
//todo: rename adjEnd to show the symmetric character
edge GraphCopy::newEdge(node v, adjEntry adjEnd, edge eOrig, CombinatorialEmbedding &E)
{
#ifdef OGDF_DEBUG
	if(v == nullptr){
		throw PreconditionViolatedException();
	}
	if(adjEnd == nullptr){
		throw PreconditionViolatedException();
	}
	if(v->graphOf() != this){
		throw PreconditionViolatedException();
	}
	if(adjEnd->graphOf() != this){
		throw PreconditionViolatedException();
	}
	if(&E.getGraph() != this){
		throw PreconditionViolatedException();
	}
	if(!m_eCopy[eOrig].empty()){
		throw PreconditionViolatedException();
	}
#endif

	//check which direction is correct
	edge e;
	if (original(v) == eOrig->source())
		e = E.splitFace(v, adjEnd);
	else
		e = E.splitFace(adjEnd, v);
	m_eIterator[e] = m_eCopy[eOrig].pushBack(e);
	m_eOrig[e] = eOrig;

	return e;
}//newedge


void GraphCopy::setEdge(edge eOrig, edge eCopy){
#ifdef OGDF_DEBUG
	if(eOrig == nullptr){
		throw PreconditionViolatedException();
	}
	if(eOrig->graphOf() != m_pGraph){
		throw PreconditionViolatedException();
	}
	if(eCopy == nullptr){
		throw PreconditionViolatedException();
	}
	if(eCopy->graphOf() != this){
		throw PreconditionViolatedException();
	}
	if(eCopy->target() != m_vCopy[eOrig->target()]){
		throw PreconditionViolatedException();
	}
	if(eCopy->source() != m_vCopy[eOrig->source()]){
		throw PreconditionViolatedException();
	}
	if(!m_eCopy[eOrig].empty()){
		throw PreconditionViolatedException();
	}
#endif

	m_eCopy[m_eOrig[eCopy] = eOrig].pushBack(eCopy);
}


void GraphCopy::insertEdgePathEmbedded(
	edge eOrig,
	CombinatorialEmbedding &E,
	const SList<adjEntry> &crossedEdges)
{
	if(m_eCopy[eOrig].size() != 0){
		FaceSetPure fsp(E);
		removeEdgePathEmbedded(E, eOrig, fsp);
	}
	m_eCopy[eOrig].clear();

	adjEntry adjSrc, adjTgt;
	SListConstIterator<adjEntry> it = crossedEdges.begin();
	SListConstIterator<adjEntry> itLast = crossedEdges.rbegin();

	// iterate over all adjacency entries in crossedEdges except for first
	// and last
	adjSrc = *it;
	for(++it; it != itLast; ++it)
	{
		adjEntry adj = *it;
		// split edge
		node u = E.split(adj->theEdge())->source();

		// determine target adjacency entry and source adjacency entry
		// in the next iteration step
		adjTgt = u->firstAdj();
		adjEntry adjSrcNext = adjTgt->succ();

		if (adjTgt != adj->twin())
			swap(adjTgt,adjSrcNext);

		// insert a new edge into the face
		edge eNew = E.splitFace(adjSrc,adjTgt);
		m_eIterator[eNew] = m_eCopy[eOrig].pushBack(eNew);
		m_eOrig[eNew] = eOrig;

		adjSrc = adjSrcNext;
	}

	// insert last edge
	edge eNew = E.splitFace(adjSrc,*it);
	m_eIterator[eNew] = m_eCopy[eOrig].pushBack(eNew);
	m_eOrig[eNew] = eOrig;
}


void GraphCopy::insertEdgePath(edge eOrig, const SList<adjEntry> &crossedEdges)
{
	if(m_eCopy[eOrig].size() != 0){
		removeEdgePath(eOrig);
	}
	node v = copy(eOrig->source());

	SListConstIterator<adjEntry> it;
	for(adjEntry adj : crossedEdges)
	{
		node u = split(adj->theEdge())->source();

		edge eNew = newEdge(v,u);
		m_eIterator[eNew] = m_eCopy[eOrig].pushBack(eNew);
		m_eOrig[eNew] = eOrig;

		v = u;
	}

	edge eNew = newEdge(v,copy(eOrig->target()));
	m_eIterator[eNew] = m_eCopy[eOrig].pushBack(eNew);
	m_eOrig[eNew] = eOrig;
}

void GraphCopy::insertEdgePath(node srcOrig, node tgtOrig, const SList<adjEntry> &crossedEdges)
{
	node v = copy(srcOrig);

	for(adjEntry adj : crossedEdges)
	{
		node u = split(adj->theEdge())->source();

		edge eNew = newEdge(v,u);
	//	m_eIterator[eNew] = m_eCopy[eOrig].pushBack(eNew);
		m_eOrig[eNew] = nullptr;

		v = u;
	}

	edge eNew = newEdge(v,copy(tgtOrig));
	//m_eIterator[eNew] = m_eCopy[eOrig].pushBack(eNew);
	m_eOrig[eNew] = nullptr;
}

//inserts crossing between two copy edges already in PlanRep
//returns newly introduced copy edge of crossed edge
//the crossing edge parameter is changed to allow iteration
//over an edge's crossings in the edge direction
//the parameter topDown describes he following:
// if the crossingEdge is running horizontally from left to right,
// is the crossedEdge direction top->down?
edge GraphCopy::insertCrossing(
	edge& crossingEdge,
	edge crossedEdge,
	bool topDown)
	//const SList<edge> &crossedCopies)
{
	//we first split the crossed edge
	edge e = split(crossedEdge);
	edge eOrig = original(crossingEdge);
	adjEntry adSource = crossingEdge->adjSource();

	//now we delete the crossing copy edge and replace it
	//by two edges adjacent to the crossing vertex
	//we have to consider the copy ordering of the
	//original edge
	//we have to keep the correct order of the adjacency entries
	//because even without a combinatorial embedding, the
	//ordering of the edges may already be fixed
	//Problem: wie erkennt man die Reihenfolge am split?
	//Man muss die Richtung der gekreuzten Kante kennen
	//=>Parameter, und hier adjSource und adjTarget trennen
	edge eNew1, eNew2;
	if (topDown)
	{
		//case 1: crossingEdge runs top-down
		eNew1 = newEdge(adSource, e->adjSource());
		eNew2 = newEdge(e->adjSource()->cyclicPred(),
			crossingEdge->adjTarget()->cyclicPred());
	}
	else
	{
		//case 2: crossingEdge runs bottom-up
		eNew1 = newEdge(adSource, e->adjSource()->cyclicPred());
		eNew2 = newEdge(e->adjSource(), crossingEdge->adjTarget()->cyclicPred());
	}//else bottom up
	//insert new edge after old entry
	m_eIterator[eNew1] = m_eCopy[eOrig].insert(eNew1, m_eIterator[crossingEdge]);
	m_eOrig[eNew1] = eOrig;
	m_eIterator[eNew2] = m_eCopy[eOrig].insert(eNew2, m_eIterator[eNew1]);
	m_eOrig[eNew2] = eOrig;
	//now we delete the input copy edge
	m_eCopy[eOrig].del(m_eIterator[crossingEdge]);
	Graph::delEdge(crossingEdge);
	crossingEdge = eNew2;

	return e;//eNew2;
}

void GraphCopy::delEdge(edge e)
{
	edge eOrig = m_eOrig[e];

	Graph::delEdge(e);
	if (eOrig == nullptr)	return;

#ifdef OGDF_DEBUG
	if(m_eCopy[eOrig].size() != 1){
		throw PreconditionViolatedException();
	}
#endif
	m_eCopy[eOrig].clear();
}


void GraphCopy::delNode(node v)
{
	node w = m_vOrig[v];
	if (w != nullptr) m_vCopy[w] = nullptr;

	Graph::delNode(v);
}


void GraphCopy::removeEdgePathEmbedded(
	CombinatorialEmbedding &E,
	edge eOrig,
	FaceSetPure &newFaces)
{
	const List<edge> &path = m_eCopy[eOrig];
#ifdef OGDF_DEBUG
	ListConstIterator<edge> testIt = path.begin();
	for(++testIt; testIt.valid(); ++testIt){
		node v = (*testIt)->source();
		if(v->degree() != 4){
			throw PreconditionViolatedException();
		}
		if(original(v->firstAdj()->theEdge()) != original(v->lastAdj()->pred()->theEdge())){
			throw PreconditionViolatedException();
		}
		if(original(v->lastAdj()->theEdge()) != original(v->firstAdj()->succ()->theEdge())){
			throw PreconditionViolatedException();
		}
	}
#endif

	ListConstIterator<edge> it = path.begin();

	newFaces.insert(E.joinFacesPure(*it));
	Graph::delEdge(*it);

	for(++it; it.valid(); ++it)
	{
		edge e = *it;
		node u = e->source();

		newFaces.remove(E.rightFace(e->adjSource()));
		newFaces.remove(E.rightFace(e->adjTarget()));

		newFaces.insert(E.joinFacesPure(e));
		Graph::delEdge(e);

		edge eIn = u->firstAdj()->theEdge();
		edge eOut = u->lastAdj()->theEdge();
		if (eIn->target() != u)
			swap(eIn,eOut);

		E.unsplit(eIn,eOut);
	}

	m_eCopy[eOrig].clear();
}


void GraphCopy::removeEdgePath(edge eOrig)
{
	const List<edge> &path = m_eCopy[eOrig];
#ifdef OGDF_DEBUG
	ListConstIterator<edge> testIt = path.begin();
	for(++testIt; testIt.valid(); ++testIt){
		node v = (*testIt)->source();
		if(v->degree() != 4){
			throw PreconditionViolatedException();
		}
	}
#endif
	ListConstIterator<edge> it = path.begin();

	Graph::delEdge(*it);

	for(++it; it.valid(); ++it)
	{
		edge e = *it;
		node u = e->source();

		Graph::delEdge(e);

		edge eIn = u->firstAdj()->theEdge();
		edge eOut = u->lastAdj()->theEdge();
		if (eIn->target() != u)
			swap(eIn,eOut);

		unsplit(eIn,eOut);
	}

	m_eCopy[eOrig].clear();
}


void GraphCopy::removeUnnecessaryCrossing(
	adjEntry adjA1,
	adjEntry adjA2,
	adjEntry adjB1,
	adjEntry adjB2)
{
	node v = adjA1->theNode();

	if(adjA1->theEdge()->source() == v)
		moveSource(adjA1->theEdge(), adjA2->twin(), before);
	else
		moveTarget(adjA1->theEdge(), adjA2->twin(), before);

	if(adjB1->theEdge()->source() == v)
		moveSource(adjB1->theEdge(), adjB2->twin(), before);
	else
		moveTarget(adjB1->theEdge(), adjB2->twin(), before);

	edge eOrigA = original(adjA1->theEdge());
	edge eOrigB = original(adjB1->theEdge());

	if (eOrigA != nullptr)
		m_eCopy[eOrigA].del(m_eIterator[adjA2->theEdge()]);
	if (eOrigB != nullptr)
		m_eCopy[eOrigB].del(m_eIterator[adjB2->theEdge()]);

	Graph::delEdge(adjB2->theEdge());
	Graph::delEdge(adjA2->theEdge());

	delNode(v);
}


bool GraphCopy::embed()
{
	return planarEmbed(*this);
}


void GraphCopy::removePseudoCrossings()
{
	node v, vSucc;
	for(v = firstNode(); v != nullptr; v = vSucc)
	{
		vSucc = v->succ();

		if (original(v) != nullptr || v->degree() != 4)
			continue;

		adjEntry adj1 = v->firstAdj();
		adjEntry adj2 = adj1->succ();
		adjEntry adj3 = adj2->succ();
		adjEntry adj4 = adj3->succ();

		if(original(adj1->theEdge()) == original(adj2->theEdge()))
			removeUnnecessaryCrossing(adj1,adj2,adj3,adj4);
		else if (original(adj2->theEdge()) == original(adj3->theEdge()))
			removeUnnecessaryCrossing(adj2,adj3,adj4,adj1);
	}
}


bool GraphCopy::consistencyCheck() const
{
	if (Graph::consistencyCheck() == false) {
		return false;
	}

	const Graph &G = *m_pGraph;

	for (node vG : G.nodes) {
		node v = m_vCopy[vG];
#ifdef OGDF_DEBUG
		if (v && v->graphOf() != this)
			return false;
#endif
		if (v && m_vOrig[v] != vG)
			return false;
	}

	for (node v : nodes) {
		node vG = m_vOrig[v];
#ifdef OGDF_DEBUG
		if (vG && vG->graphOf() != &G)
			return false;
#endif
		if (vG && m_vCopy[vG] != v)
			return false;
	}

	for (edge eG : G.edges) {
		const List<edge> &path = m_eCopy[eG];
		ListConstIterator<edge> it;
		for (edge e : path) {
#ifdef OGDF_DEBUG
			if (e->graphOf() != this)
				return false;
#endif
			if (m_eOrig[e] != eG)
				return false;
		}
	}

#ifdef OGDF_DEBUG
	for (edge e : edges) {
		edge eG = m_eOrig[e];
		if (eG
		 && eG->graphOf() != &G)
			return false;
	}
#endif

	return true;
}



} // end namespace ogdf
