
#include <ogdf/internal/planarity/FixEdgeInserterCore.h>
#include <ogdf/basic/FaceSet.h>
#include <ogdf/basic/Queue.h>


namespace ogdf {

	//---------------------------------------------------------
	// FEICrossingsBucket
	// bucket function for sorting edges by decreasing number
	// of crossings
	class FEICrossingsBucket : public BucketFunc<edge>
	{
		const PlanRepLight *m_pPG;

	public:
		FEICrossingsBucket(const PlanRepLight *pPG) :
			m_pPG(pPG) { }

		int getBucket(const edge &e) {
			return -m_pPG->chain(e).size();
		}
	};


	void FixEdgeInserterCore::init(CombinatorialEmbedding &E)
	{
		m_dual.clear();
		m_primalAdj.init(m_dual);
		m_nodeOf.init(E);
	}

	void FixEdgeInserterCore::cleanup()
	{
		delete m_newFaces;
		delete m_delFaces;

		m_nodeOf.init();
		m_primalAdj.init();
		m_dual.clear();
	}

	void FixEdgeInserterUMLCore::init(CombinatorialEmbedding &E)
	{
		FixEdgeInserterCore::init(E);
		m_primalIsGen.init(m_dual,false);
	}

	void FixEdgeInserterUMLCore::cleanup()
	{
		m_primalIsGen.init();
		FixEdgeInserterCore::cleanup();
	}


	//--------------------------------------------------------------------
	// actual algorithm call
	//--------------------------------------------------------------------
	Module::ReturnType FixEdgeInserterCore::call(
		const Array<edge> &origEdges,
		bool keepEmbedding,
		RemoveReinsertType rrPost,
		double percentMostCrossed)
	{
		double T;
		usedTime(T);

		Module::ReturnType retValue = Module::retFeasible;
		m_runsPostprocessing = 0;

		if(!keepEmbedding) m_pr.embed();
		OGDF_ASSERT(m_pr.representsCombEmbedding() == true);

		if (origEdges.size() == 0)
			return Module::retOptimal;  // nothing to do

		// initialization
		CombinatorialEmbedding E(m_pr);  // embedding of PG

		init(E);
		constructDual(E);

		// m_delFaces and m_newFaces are used by removeEdge()
		// if we can't allocate memory for them, we throw an exception
		if (rrPost != rrNone) {
			m_delFaces = new FaceSetSimple(E);
			if (m_delFaces == 0)
				OGDF_THROW(InsufficientMemoryException);

			m_newFaces = new FaceSetPure(E);
			if (m_newFaces == 0) {
				delete m_delFaces;
				OGDF_THROW(InsufficientMemoryException);
			}

		// no postprocessing -> no removeEdge()
		} else {
			m_delFaces = 0;
			m_newFaces = 0;
		}

		SListPure<edge> currentOrigEdges;
		if(rrPost == rrIncremental) {
			edge e;
			forall_edges(e,m_pr)
				currentOrigEdges.pushBack(m_pr.original(e));
		}

		// insertion of edges
		bool doIncrementalPostprocessing =
			( rrPost == rrIncremental || rrPost == rrIncInserted );
		for(int i = origEdges.low(); i <= origEdges.high(); ++i)
		{
			edge eOrig = origEdges[i];
			storeTypeOfCurrentEdge(eOrig);
			//int eSubGraph = 0;  // edgeSubGraphs-data of eOrig
			//if(edgeSubGraphs!=0) eSubGraph = (*edgeSubGraphs)[eOrig];

			SList<adjEntry> crossed;
			if(m_pCost != 0) {
				findWeightedShortestPath(E, eOrig, crossed);
			} else {
				findShortestPath(E, eOrig, crossed);
			}

			insertEdge(E, eOrig, crossed);

			if(doIncrementalPostprocessing) {
				currentOrigEdges.pushBack(eOrig);

				bool improved;
				do {
					++m_runsPostprocessing;
					improved = false;

					SListConstIterator<edge> itRR;
					for(itRR = currentOrigEdges.begin(); itRR.valid(); ++itRR)
					{
						edge eOrigRR = *itRR;

						int pathLength = (m_pCost != 0) ? costCrossed(eOrigRR) : (m_pr.chain(eOrigRR).size() - 1);
						if (pathLength == 0) continue; // cannot improve

						removeEdge(E, eOrigRR);

						storeTypeOfCurrentEdge(eOrigRR);

						// try to find a better insertion path
						SList<adjEntry> crossed;
						if(m_pCost != 0) {
							findWeightedShortestPath(E, eOrigRR, crossed);
						} else {
							findShortestPath(E, eOrigRR, crossed);
						}

						// re-insert edge (insertion path cannot be longer)
						insertEdge(E, eOrigRR, crossed);

						int newPathLength = (m_pCost != 0) ? costCrossed(eOrigRR) : (m_pr.chain(eOrigRR).size() - 1);
						OGDF_ASSERT(newPathLength <= pathLength);

						if(newPathLength < pathLength)
							improved = true;
					}
				} while (improved);
			}
		}

		if(!doIncrementalPostprocessing) {
			// postprocessing (remove-reinsert heuristc)
			const int m = m_pr.original().numberOfEdges();
			SListPure<edge> rrEdges;

			switch(rrPost)
			{
			case rrAll:
			case rrMostCrossed:
				for(int i = m_pr.startEdge(); i < m_pr.stopEdge(); ++i)
					rrEdges.pushBack(m_pr.e(i));
				break;

			case rrInserted:
				for(int i = origEdges.low(); i <= origEdges.high(); ++i)
					rrEdges.pushBack(origEdges[i]);
				break;

			case rrNone:
			case rrIncremental:
			case rrIncInserted:
				break;
			}

			// marks the end of the interval of rrEdges over which we iterate
			// initially set to invalid iterator which means all edges
			SListConstIterator<edge> itStop;

			bool improved;
			do {
				// abort postprocessing if time limit reached
				if (m_timeLimit >= 0 && m_timeLimit <= usedTime(T)) {
					retValue = Module::retTimeoutFeasible;
					break;
				}

				++m_runsPostprocessing;
				improved = false;

				if(rrPost == rrMostCrossed)
				{
					FEICrossingsBucket bucket(&m_pr);
					rrEdges.bucketSort(bucket);

					const int num = int(0.01 * percentMostCrossed * m);
					itStop = rrEdges.get(num);
				}

				SListConstIterator<edge> it;
				for(it = rrEdges.begin(); it != itStop; ++it)
				{
					edge eOrig = *it;

					int pathLength = (m_pCost != 0) ? costCrossed(eOrig) : (m_pr.chain(eOrig).size() - 1);
					if (pathLength == 0) continue; // cannot improve

					removeEdge(E, eOrig);

					storeTypeOfCurrentEdge(eOrig);

					// try to find a better insertion path
					SList<adjEntry> crossed;
					if(m_pCost != 0) {
						findWeightedShortestPath(E, eOrig, crossed);
					} else {
						findShortestPath(E, eOrig, crossed);
					}

					// re-insert edge (insertion path cannot be longer)
					insertEdge(E, eOrig, crossed);

					// we cannot find a shortest path that is longer than before!
					int newPathLength = (m_pCost != 0) ? costCrossed(eOrig) : (m_pr.chain(eOrig).size() - 1);
					OGDF_ASSERT(newPathLength <= pathLength);

					if(newPathLength < pathLength)
						improved = true;
				}
			} while (improved);
		}

		// verify computed planarization
		OGDF_ASSERT(m_pr.representsCombEmbedding());

		// free resources
		cleanup();

		return retValue;
	}


	static edge crossedEdge(adjEntry adj)
	{
		edge e = adj->theEdge();

		adj = adj->cyclicSucc();
		while(adj->theEdge() == e)
			adj = adj->cyclicSucc();

		return adj->theEdge();
	}


	int FixEdgeInserterCore::costCrossed(edge eOrig) const
	{
		int c = 0;

		const List<edge> &L = m_pr.chain(eOrig);

		ListConstIterator<edge> it = L.begin();
		if(m_pSubgraph != 0) {
			for(++it; it.valid(); ++it) {
				int counter = 0;
				edge e = m_pr.original(crossedEdge((*it)->adjSource()));
				for(int i = 0; i < 32; i++)
					if((*m_pSubgraph)[eOrig] & (*m_pSubgraph)[e] & (1<<i))
						counter++;
				c += counter * (*m_pCost)[e];
			}

		} else {
			for(++it; it.valid(); ++it) {
				c += (*m_pCost)[m_pr.original(crossedEdge((*it)->adjSource()))];
			}
		}

		return c;
	}


	//---------------------------------------------------------
	// construct dual graph
	// assumes that m_pDual, m_primalAdj and m_nodeOf are already constructed
	//
	void FixEdgeInserterCore::constructDual(const CombinatorialEmbedding &E)
	{
		// insert a node in the dual graph for each face in E
		face f;
		forall_faces(f,E)
			m_nodeOf[f] = m_dual.newNode();


		// Insert an edge into the dual graph for each adjacency entry in E.
		// The edges are directed from the left face to the right face.
		node v;
		forall_nodes(v,m_pr)
		{
			adjEntry adj;
			forall_adj(adj,v)
			{
				// Do not insert edges into dual if crossing the original edge
				// is forbidden
				if(m_pForbidden && (*m_pForbidden)[m_pr.original(adj->theEdge())] == true)
					continue;

				node vLeft  = m_nodeOf[E.leftFace (adj)];
				node vRight = m_nodeOf[E.rightFace(adj)];

				m_primalAdj[m_dual.newEdge(vLeft,vRight)] = adj;
			}
		}

		// Augment the dual graph by two new vertices. These are used temporarily
		// when searching for a shortest path in the dual graph.
		m_vS = m_dual.newNode();
		m_vT = m_dual.newNode();
	}


	//---------------------------------------------------------
	// construct dual graph, marks dual edges corresponding to generalization
	// in m_primalIsGen
	// assumes that m_pDual, m_primalAdj and m_nodeOf are already constructed
	//
	void FixEdgeInserterUMLCore::constructDual(const CombinatorialEmbedding &E)
	{
		// insert a node in the dual graph for each face in E
		face f;
		forall_faces(f,E)
			m_nodeOf[f] = m_dual.newNode();


		// Insert an edge into the dual graph for each adjacency entry in E.
		// The edges are directed from the left face to the right face.
		node v;
		forall_nodes(v,m_pr)
		{
			adjEntry adj;
			forall_adj(adj,v)
			{
				node vLeft  = m_nodeOf[E.leftFace (adj)];
				node vRight = m_nodeOf[E.rightFace(adj)];

				edge e = m_dual.newEdge(vLeft,vRight);
				m_primalAdj[e] = adj;

				// mark dual edges corresponding to generalizations
				if (m_pr.typeOf(adj->theEdge()) == Graph::generalization)
					m_primalIsGen[e] = true;
			}
		}

		// Augment the dual graph by two new vertices. These are used temporarily
		// when searching for a shortest path in the dual graph.
		m_vS = m_dual.newNode();
		m_vT = m_dual.newNode();
	}


	//---------------------------------------------------------
	// finds a shortest path in the dual graph augmented by s and t (represented
	// by m_vS and m_vT); returns list of crossed adjacency entries (corresponding
	// to used edges in the dual) in crossed.
	//

	void FixEdgeInserterCore::appendCandidates(QueuePure<edge> &queue, node v)
	{
		edge e;
		forall_adj_edges(e,v) {
			if(v == e->source())
				queue.append(e);
		}
	}

	void FixEdgeInserterUMLCore::appendCandidates(QueuePure<edge> &queue, node v)
	{
		edge e;
		forall_adj_edges(e,v) {
			if(v == e->source() &&
				(m_typeOfCurrentEdge != Graph::generalization || m_primalIsGen[e] == false))
			{
				queue.append(e);
			}
		}
	}

	void FixEdgeInserterCore::findShortestPath(const CombinatorialEmbedding &E, edge eOrig, SList<adjEntry> &crossed)
	{
		node s = m_pr.copy(eOrig->source());
		node t = m_pr.copy(eOrig->target());
		OGDF_ASSERT(s != t);

		NodeArray<edge> spPred(m_dual,0);
		QueuePure<edge> queue;
		int oldIdCount = m_dual.maxEdgeIndex();

		// augment dual by edges from s to all adjacent faces of s ...
		adjEntry adj;
		forall_adj(adj,s) {
			// starting edges of bfs-search are all edges leaving s
			edge eDual = m_dual.newEdge(m_vS, m_nodeOf[E.rightFace(adj)]);
			m_primalAdj[eDual] = adj;
			queue.append(eDual);
		}

		// ... and from all adjacent faces of t to t
		forall_adj(adj,t) {
			edge eDual = m_dual.newEdge(m_nodeOf[E.rightFace(adj)], m_vT);
			m_primalAdj[eDual] = adj;
		}

		// actual search (using bfs on directed dual)
		for( ; ;)
		{
			// next candidate edge
			edge eCand = queue.pop();
			node v = eCand->target();

			// leads to an unvisited node?
			if (spPred[v] == 0)
			{
				// yes, then we set v's predecessor in search tree
				spPred[v] = eCand;

				// have we reached t ...
				if (v == m_vT)
				{
					// ... then search is done.
					// constructed list of used edges (translated to crossed
					// adjacency entries in PG) from t back to s (including first
					// and last!)

					do {
						edge eDual = spPred[v];
						crossed.pushFront(m_primalAdj[eDual]);
						v = eDual->source();
					} while(v != m_vS);

					break;
				}

				// append next candidate edges to queue (all edges leaving v)
				appendCandidates(queue, v);
			}
		}


		// remove augmented edges again
		while ((adj = m_vS->firstAdj()) != 0)
			m_dual.delEdge(adj->theEdge());

		while ((adj = m_vT->firstAdj()) != 0)
			m_dual.delEdge(adj->theEdge());

		m_dual.resetEdgeIdCount(oldIdCount);
	}


	int FixEdgeInserterCore::getCost(edge e, int stSubgraph) const
	{
		edge eOrig = m_pr.original(e);
		if(m_pSubgraph == 0)
			return (eOrig == 0) ? 0 : (*m_pCost)[eOrig];

		int edgeCost = 0;
		if(eOrig != 0) {
			for(int i = 0; i < 32; i++) {
				if((((*m_pSubgraph)[eOrig] & (1 << i)) != 0) && ((stSubgraph & (1 << i)) != 0))
					edgeCost++;
			}
			edgeCost *= (*m_pCost)[eOrig];
			edgeCost *= 10000;
			if(edgeCost == 0)
				edgeCost = 1;
		}
		return edgeCost;
	}


	//---------------------------------------------------------
	// finds a weighted shortest path in the dual graph augmented by s and t
	// (represented by m_vS and m_vT) using edges weights given by costOrig;
	// returns list of crossed adjacency entries (corresponding
	// to used edges in the dual) in crossed.
	//
	// running time: O(|dual| + L + C),
	//   where L ist the weighted length of the insertion path and C the
	//   maximum cost of an edge
	//

	void FixEdgeInserterCore::appendCandidates(
		Array<SListPure<edge> > &nodesAtDist, EdgeArray<int> &costDual, int maxCost, node v, int currentDist)
	{
		edge e;
		forall_adj_edges(e,v) {
			if(v == e->source()) {
				int listPos = (currentDist + costDual[e]) % maxCost;
				nodesAtDist[listPos].pushBack(e);
			}
		}
	}

	void FixEdgeInserterUMLCore::appendCandidates(
		Array<SListPure<edge> > &nodesAtDist, EdgeArray<int> &costDual, int maxCost, node v, int currentDist)
	{
		edge e;
		forall_adj_edges(e,v) {
			if(v == e->source() &&
				(m_typeOfCurrentEdge != Graph::generalization || m_primalIsGen[e] == false))
			{
				int listPos = (currentDist + costDual[e]) % maxCost;
				nodesAtDist[listPos].pushBack(e);
			}
		}
	}

	void FixEdgeInserterCore::findWeightedShortestPath(const CombinatorialEmbedding &E, edge eOrig, SList<adjEntry> &crossed)
	{
		node s = m_pr.copy(eOrig->source());
		node t = m_pr.copy(eOrig->target());
		OGDF_ASSERT(s != t);

		int eSubgraph = (m_pSubgraph != 0) ? (*m_pSubgraph)[eOrig] : 0;

		EdgeArray<int> costDual(m_dual, 0);
		int maxCost = 0;
		edge eDual;
		forall_edges(eDual, m_dual) {
			int c = getCost(m_primalAdj[eDual]->theEdge(), eSubgraph);
			costDual[eDual] = c;
			if (c > maxCost)
				maxCost = c;
		}

		++maxCost;
		Array<SListPure<edge> > nodesAtDist(maxCost);

		NodeArray<edge> spPred(m_dual,0);

		int oldIdCount = m_dual.maxEdgeIndex();

		// augment dual by edges from s to all adjacent faces of s ...
		adjEntry adj;
		forall_adj(adj,s) {
			// starting edges of bfs-search are all edges leaving s
			edge eDual = m_dual.newEdge(m_vS, m_nodeOf[E.rightFace(adj)]);
			m_primalAdj[eDual] = adj;
			nodesAtDist[0].pushBack(eDual);
		}

		// ... and from all adjacent faces of t to t
		forall_adj(adj,t) {
			edge eDual = m_dual.newEdge(m_nodeOf[E.rightFace(adj)], m_vT);
			m_primalAdj[eDual] = adj;
		}

		// actual search (using extended bfs on directed dual)
		int currentDist = 0;

		for( ; ; )
		{
			// next candidate edge
			while(nodesAtDist[currentDist % maxCost].empty())
				++currentDist;

			edge eCand = nodesAtDist[currentDist % maxCost].popFrontRet();
			node v = eCand->target();

			// leads to an unvisited node?
			if (spPred[v] == 0)
			{
				// yes, then we set v's predecessor in search tree
				spPred[v] = eCand;

				// have we reached t ...
				if (v == m_vT)
				{
					// ... then search is done.
					// constructed list of used edges (translated to crossed
					// adjacency entries in PG) from t back to s (including first
					// and last!)

					do {
						edge eDual = spPred[v];
						crossed.pushFront(m_primalAdj[eDual]);
						v = eDual->source();
					} while(v != m_vS);

					break;
				}

				// append next candidate edges to queue (all edges leaving v)
				appendCandidates(nodesAtDist, costDual, maxCost, v, currentDist);
			}
		}

		// remove augmented edges again
		while ((adj = m_vS->firstAdj()) != 0)
			m_dual.delEdge(adj->theEdge());

		while ((adj = m_vT->firstAdj()) != 0)
			m_dual.delEdge(adj->theEdge());

		m_dual.resetEdgeIdCount(oldIdCount);
	}


	//---------------------------------------------------------
	// inserts edge e according to insertion path crossed.
	// updates embeding and dual graph
	//

	void FixEdgeInserterCore::insertEdgesIntoDual(const CombinatorialEmbedding &E, adjEntry adjSrc)
	{
		face f = E.rightFace(adjSrc);
		node vRight = m_nodeOf[f];

		adjEntry adj1 = f->firstAdj(), adj = adj1;
		do {
			if(m_pForbidden && (*m_pForbidden)[m_pr.original(adj->theEdge())] == true)
				continue;

			node vLeft = m_nodeOf[E.leftFace(adj)];

			edge eLR = m_dual.newEdge(vLeft,vRight);
			m_primalAdj[eLR] = adj;

			edge eRL = m_dual.newEdge(vRight,vLeft);
			m_primalAdj[eRL] = adj->twin();
		}
		while((adj = adj->faceCycleSucc()) != adj1);

		// the other face adjacent to *itEdge ...
		f = E.rightFace(adjSrc->twin());
		vRight = m_nodeOf[f];

		adj1 = f->firstAdj();
		adj = adj1;
		do {
			if(m_pForbidden && (*m_pForbidden)[m_pr.original(adj->theEdge())] == true)
				continue;

			node vLeft = m_nodeOf[E.leftFace(adj)];

			edge eLR = m_dual.newEdge(vLeft,vRight);
			m_primalAdj[eLR] = adj;

			edge eRL = m_dual.newEdge(vRight,vLeft);
			m_primalAdj[eRL] = adj->twin();
		}
		while((adj = adj->faceCycleSucc()) != adj1);
	}

	void FixEdgeInserterUMLCore::insertEdgesIntoDual(const CombinatorialEmbedding &E, adjEntry adjSrc)
	{
		face f = E.rightFace(adjSrc);
		node vRight = m_nodeOf[f];

		adjEntry adj1 = f->firstAdj(), adj = adj1;
		do {
			node vLeft = m_nodeOf[E.leftFace(adj)];

			edge eLR = m_dual.newEdge(vLeft,vRight);
			m_primalAdj[eLR] = adj;

			edge eRL = m_dual.newEdge(vRight,vLeft);
			m_primalAdj[eRL] = adj->twin();

			if(m_pr.typeOf(adj->theEdge()) == Graph::generalization)
				m_primalIsGen[eLR] = m_primalIsGen[eRL] = true;
		}
		while((adj = adj->faceCycleSucc()) != adj1);

		// the other face adjacent to *itEdge ...
		f = E.rightFace(adjSrc->twin());
		vRight = m_nodeOf[f];

		adj1 = f->firstAdj();
		adj = adj1;
		do {
			node vLeft = m_nodeOf[E.leftFace(adj)];

			edge eLR = m_dual.newEdge(vLeft,vRight);
			m_primalAdj[eLR] = adj;

			edge eRL = m_dual.newEdge(vRight,vLeft);
			m_primalAdj[eRL] = adj->twin();

			if(m_pr.typeOf(adj->theEdge()) == Graph::generalization)
				m_primalIsGen[eLR] = m_primalIsGen[eRL] = true;
		}
		while((adj = adj->faceCycleSucc()) != adj1);
	}

	void FixEdgeInserterCore::insertEdge(CombinatorialEmbedding &E, edge eOrig, const SList<adjEntry> &crossed)
	{
		// remove dual nodes on insertion path
		SListConstIterator<adjEntry> it;
		for(it = crossed.begin(); it != crossed.rbegin(); ++it) {
			m_dual.delNode(m_nodeOf[E.rightFace(*it)]);
		}

		// update primal
		m_pr.insertEdgePathEmbedded(eOrig,E,crossed);

		// insert new face nodes into dual
		const List<edge> &path = m_pr.chain(eOrig);
		ListConstIterator<edge> itEdge;
		for(itEdge = path.begin(); itEdge.valid(); ++itEdge)
		{
			adjEntry adj = (*itEdge)->adjSource();
			m_nodeOf[E.leftFace (adj)] = m_dual.newNode();
			m_nodeOf[E.rightFace(adj)] = m_dual.newNode();
		}

		// insert new edges into dual
		for(itEdge = path.begin(); itEdge.valid(); ++itEdge)
			insertEdgesIntoDual(E, (*itEdge)->adjSource());
	}

	//---------------------------------------------------------
	// removes edge eOrig; updates embedding and dual graph
	//

	void FixEdgeInserterCore::insertEdgesIntoDualAfterRemove(const CombinatorialEmbedding &E, face f)
	{
		node vRight = m_nodeOf[f];

		adjEntry adj1 = f->firstAdj(), adj = adj1;
		do {
			if(m_pForbidden && (*m_pForbidden)[m_pr.original(adj->theEdge())] == true)
				continue;

			node vLeft = m_nodeOf[E.leftFace(adj)];

			edge eLR = m_dual.newEdge(vLeft,vRight);
			m_primalAdj[eLR] = adj;

			edge eRL = m_dual.newEdge(vRight,vLeft);
			m_primalAdj[eRL] = adj->twin();
		}
		while((adj = adj->faceCycleSucc()) != adj1);
	}

	void FixEdgeInserterUMLCore::insertEdgesIntoDualAfterRemove(const CombinatorialEmbedding &E, face f)
	{
		node vRight = m_nodeOf[f];

		adjEntry adj1 = f->firstAdj(), adj = adj1;
		do {
			node vLeft = m_nodeOf[E.leftFace(adj)];

			edge eLR = m_dual.newEdge(vLeft,vRight);
			m_primalAdj[eLR] = adj;

			edge eRL = m_dual.newEdge(vRight,vLeft);
			m_primalAdj[eRL] = adj->twin();

			if(m_pr.typeOf(adj->theEdge()) == Graph::generalization)
			{
				m_primalIsGen[eLR] = m_primalIsGen[eRL] = true;
			}
		}
		while((adj = adj->faceCycleSucc()) != adj1);
	}

	void FixEdgeInserterCore::removeEdge(CombinatorialEmbedding &E, edge eOrig)
	{
		const List<edge> &path = m_pr.chain(eOrig);
		ListConstIterator<edge> itEdge;
		for(itEdge = path.begin(); itEdge.valid(); ++itEdge)
		{
			adjEntry adj = (*itEdge)->adjSource();
			m_delFaces->insert(E.leftFace  (adj));
			m_delFaces->insert(E.rightFace (adj));
		}

		// delete all corresponding nodes in dual
		SListConstIterator<face> itsF;
		for(itsF = m_delFaces->faces().begin(); itsF.valid(); ++itsF)
			m_dual.delNode(m_nodeOf[*itsF]);

		m_delFaces->clear();

		// remove edge path from PG
		m_pr.removeEdgePathEmbedded(E,eOrig,*m_newFaces);

		// update dual
		// insert new nodes
		ListConstIterator<face> itF;
		for(itF = m_newFaces->faces().begin(); itF.valid(); ++itF) {
			m_nodeOf[*itF] = m_dual.newNode();
		}

		// insert new edges into dual
		for(itF = m_newFaces->faces().begin(); itF.valid(); ++itF)
			insertEdgesIntoDualAfterRemove(E, *itF);

		m_newFaces->clear();
	}

}
