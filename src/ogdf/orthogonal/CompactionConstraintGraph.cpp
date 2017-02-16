/** \file
 * \brief Implementation of class CompactionConstraintGraphBase.
 *
 * Represents base class for CompactionConstraintGraph<ATYPE>
 *
 * \author Carsten Gutwenger
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


#include <ogdf/orthogonal/CompactionConstraintGraph.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/basic/extended_graph_alg.h>


namespace ogdf {


// constructor for orthogonal representation
// builds constraint graph with basic arcs
CompactionConstraintGraphBase::CompactionConstraintGraphBase(
	const OrthoRep &OR,
	const PlanRep &PG,
	OrthoDir arcDir,
	int costGen,
	int costAssoc,
	bool align
) :
	m_pathNode(OR),
	m_edgeToBasicArc(OR,nullptr)
{
	OGDF_ASSERT(&PG == &(const Graph &)OR);

	m_path        .init(*this);
	m_cost        .init(*this,costAssoc);
	m_type        .init(*this,ConstraintEdgeType::BasicArc);
	m_verticalGen .init(PG, false);
	m_verticalArc .init(*this, false);
	m_border      .init(*this, false);
	m_alignmentArc.init(*this, false);
	m_pathToEdge  .init(*this, nullptr);
	m_originalEdge.init(*this, nullptr);

	m_pPR       = &PG;//only used in detecting cage visibility arcs
	m_pOR       = &OR;
	m_align     = align;
	m_arcDir    = arcDir;
	m_oppArcDir = OR.oppDir(arcDir);
	m_edgeCost[static_cast<int>(Graph::EdgeType::generalization)] = costGen;
	m_edgeCost[static_cast<int>(Graph::EdgeType::association)   ] = costAssoc;

	for(edge e : PG.edges)
	{
		if ((PG.typeOf(e) == Graph::EdgeType::generalization) && (!PG.isExpansionEdge(e)))
			m_verticalGen[e] = true;
	}//for

	insertPathVertices(PG);
	insertBasicArcs(PG);
}


// insert vertex for each segment
void CompactionConstraintGraphBase::insertPathVertices(const PlanRep &PG)
{
	NodeArray<node> genOpposite(PG,nullptr);

	for(node v : PG.nodes)
	{
		const OrthoRep::VertexInfoUML *vi = m_pOR->cageInfo(v);
		if (vi == nullptr || PG.typeOf(v) == Graph::NodeType::generalizationMerger) continue;

		adjEntry adjGen = vi->m_side[static_cast<int>(m_arcDir)   ].m_adjGen;
		adjEntry adjOpp = vi->m_side[static_cast<int>(m_oppArcDir)].m_adjGen;
		if (adjGen != nullptr && adjOpp != nullptr)
		{
			node v1 = adjGen->theNode();
			node v2 = adjOpp->theNode();
			genOpposite[genOpposite[v1] = v2] = v1;
		}
	}

	//TODO: hier abfragen, ob Kantensegment einer Originalkante
	//und diese hat Multikantenstatus => man muss sich die Abstaende am Rand merken
	//und auf alle Segmente anwenden

	NodeArray<bool> visited(PG,false);

	for(node v : PG.nodes)
	{
		if (!visited[v]) {
			node pathVertex = newNode();

			dfsInsertPathVertex(v, pathVertex, visited, genOpposite);

			//test for multi sep
			//if exact two nodes in path, edge between is original edge segment,
			//save this to recall the minsep for all segments later over original
			//muss man hier originaledge als rueckgabe reintun???
			if ((m_path[pathVertex].size() == 2) && m_pathToEdge[pathVertex])
			{

			}//if original segment
			else m_pathToEdge[pathVertex] = nullptr;
		}
	}
}


// insert all graph vertices into segment pathVertex
void CompactionConstraintGraphBase::dfsInsertPathVertex(
	node v,
	node pathVertex,
	NodeArray<bool> &visited,
	const NodeArray<node> &genOpposite)
{
	visited[v] = true;
	m_path[pathVertex].pushFront(v);
	m_pathNode[v] = pathVertex;

	for(adjEntry adj : v->adjEntries)
	{
		OrthoDir dirAdj = m_pOR->direction(adj);
		OGDF_ASSERT(dirAdj != OrthoDir::Undefined);

		if (dirAdj != m_arcDir && dirAdj != m_oppArcDir) {
			//for multiedges, only useful if only one edge considered on path
			//maybe zero if no original edge exists
			if (!m_pathToEdge[pathVertex])
			{
				//only reset later for multi edge segments
				m_pathToEdge[pathVertex] = m_pPR->original(adj->theEdge());
				//used for all vertices

			}

			node w = adj->theEdge()->opposite(v);
			if (!visited[w])
				dfsInsertPathVertex(w, pathVertex, visited, genOpposite);
		}
	}

	node w = genOpposite[v];
	if (w != nullptr && !visited[w])
		dfsInsertPathVertex(w, pathVertex, visited, genOpposite);
}



//
// insert an arc for each edge with direction m_arcDir
void CompactionConstraintGraphBase::insertBasicArcs(const PlanRep &PG)
{
	const Graph &G = *m_pOR;

	for(node v : G.nodes)
	{
		node start = m_pathNode[v];

		for(adjEntry adj : v->adjEntries) {
			if (m_pOR->direction(adj) == m_arcDir) {
				edge e = newEdge(start, m_pathNode[adj->theEdge()->opposite(v)]);
				m_edgeToBasicArc[adj] = e;

				m_cost[e] = m_edgeCost[static_cast<int>(PG.typeOf(adj->theEdge()))];

				//try to pull nodes up in hierarchies
				if ( (PG.typeOf(adj->theEdge()) == Graph::EdgeType::generalization) &&
					(PG.typeOf(adj->theEdge()->target()) == Graph::NodeType::generalizationExpander) &&
					!(PG.isExpansionEdge(adj->theEdge()))
					)
				{
					if (m_align)
					{
						//got to be higher than vertexarccost*doublebendfactor
						m_cost[e] = 4000*m_cost[e]; //use parameter later corresponding
						m_alignmentArc[e] = true;
					}//if align
					//to compconsgraph::doublebendfactor
					else m_cost[e] = 2*m_cost[e];
				}

				//set generalization type
				if (verticalGen(adj->theEdge())) m_verticalArc[e] = true;
				//set onborder
				if (PG.isDegreeExpansionEdge(adj->theEdge()))
				{
					edge borderE = adj->theEdge();
					node v1 = borderE->source();
					node v2 = borderE->target();
					m_border[e] = ((v1->degree()>2) && (v2->degree()>2) ? 2 : 1);
				}

			}
		}
	}
}


// embeds constraint graph such that all sources and sinks lie in a common
// face
void CompactionConstraintGraphBase::embed()
{
	NodeArray<bool> onExternal(*this,false);
	const CombinatorialEmbedding &E = *m_pOR;
	face fExternal = E.externalFace();

	for(adjEntry adj : fExternal->entries)
		onExternal[m_pathNode[adj->theNode()]] = true;

	// compute lists of sources and sinks
	SList<node> sources, sinks;

	for(node v : nodes) {
		if (onExternal[v]) {
			if (v->indeg() == 0)
				sources.pushBack(v);
			if (v->outdeg() == 0)
				sinks.pushBack(v);
		}
	}

	// determine super source and super sink
	node s,t;
	if (sources.size() > 1)
	{
		s = newNode();
		for (node v : sources)
			newEdge(s,v);
	}
	else
		s = sources.front();

	if (sinks.size() > 1)
	{
		t = newNode();
		for (node v : sinks)
			newEdge(v,t);
	}
	else
		t = sinks.front();

	edge st = newEdge(s,t);

	bool isPlanar = planarEmbed(*this);
	if (!isPlanar) OGDF_THROW(AlgorithmFailureException);


	delEdge(st);
	if (sources.size() > 1)
		delNode(s);
	if (sinks.size() > 1)
		delNode(t);
}


// computes topological numbering on the segments of the constraint graph.
// Usage: If used on the basic (and vertex size) arcs, the numbering can be
//   used in order to serve as sorting criteria for respecting the given
//   embedding, e.g., when computing visibility arcs and allowing edges
//   with length 0.
void CompactionConstraintGraphBase::computeTopologicalSegmentNum(
	NodeArray<int> &topNum)
{
	NodeArray<int> indeg(*this);
	StackPure<node> sources;

	for(node v : nodes) {
		topNum[v] = 0;
		indeg[v] = v->indeg();
		if(indeg[v] == 0)
			sources.push(v);
	}

	while(!sources.empty())
	{
		node v = sources.pop();

		for(adjEntry adj : v->adjEntries) {
			edge e = adj->theEdge();
			if(e->source() != v) continue;

			node w = e->target();

			if (topNum[w] < topNum[v] + 1)
				topNum[w] = topNum[v] + 1;

			if (--indeg[w] == 0)
				sources.push(w);
		}
	}
}

// remove "arcs" from visibArcs which we already have in the constraint graph
// (as basic arcs)
void CompactionConstraintGraphBase::removeRedundantVisibArcs(
	SListPure<Tuple2<node,node> > &visibArcs)
{
	// bucket sort list of all edges
	SListPure<edge> all;
	allEdges(all);
	parallelFreeSort(*this,all);

	// bucket sort visibArcs
	struct : public BucketFunc<Tuple2<node,node>> {
		int getBucket(const Tuple2<node,node> &t) override {
			return t.x1()->index();
		}
	} bucketSrc;
	visibArcs.bucketSort(0,maxNodeIndex(),bucketSrc);

	struct : public BucketFunc<Tuple2<node,node>> {
		int getBucket(const Tuple2<node,node> &t) override {
			return t.x2()->index();
		}
	} bucketTgt;
	visibArcs.bucketSort(0,maxNodeIndex(),bucketTgt);

	// now, in both lists, arcs are sorted by increasing target index,
	// and arcs with the same target index by increasing source index.
	SListConstIterator<edge> itAll = all.begin();
	SListIterator<Tuple2<node,node> > it, itNext, itPrev;

	// for each arc in visibArcs, we check if it is also contained in list all
	for(it = visibArcs.begin(); it.valid(); it = itNext)
	{
		// required since we delete from the list we traverse
		itNext = it.succ();
		int i = (*it).x1()->index();
		int j = (*it).x2()->index();

		// skip all arcs with smaller target index
		while(itAll.valid() && (*itAll)->target()->index() < j)
			++itAll;

		// no more arcs => no more duplicates, so return
		if (!itAll.valid()) break;

		// if target index is j, we also skip all arcs with target index i
		// and source index smaller than i
		while(itAll.valid() && (*itAll)->target()->index() == j && (*itAll)->source()->index() < i)
			++itAll;

		// no more arcs => no more duplicates, so return
		if (!itAll.valid()) break;

		// if (i,j) is already present, we delete it from visibArcs
		if ((*itAll)->source()->index() == i &&
			(*itAll)->target()->index() == j)
		{
			//visibArcs.del(it);
			if (itPrev.valid())
				visibArcs.delSucc(itPrev);
			else
				visibArcs.popFront();
		} else
			itPrev = it;
	}//for visibArcs

	// CHECK for
	//special treatment for cage visibility
	//two cases: input node cage: just compare arbitrary node
	//           merger cage: check first if there are mergers
	itPrev = nullptr;
	for(it = visibArcs.begin(); it.valid(); it = itNext)
	{

		itNext = it.succ();

		OGDF_ASSERT(!m_path[(*it).x1()].empty());
		OGDF_ASSERT(!m_path[(*it).x1()].empty());

		node boundRepresentant1 = m_path[(*it).x1()].front();
		node boundRepresentant2 = m_path[(*it).x2()].front();
		node en1 = m_pPR->expandedNode(boundRepresentant1);
		node en2 = m_pPR->expandedNode(boundRepresentant2);
		//do not allow visibility constraints in fixed cages
		//due to non-planarity with middle position constraints

		if ( ( en1 && en2 ) && ( en1 == en2) )
		{
			if (itPrev.valid()) visibArcs.delSucc(itPrev);
			else visibArcs.popFront();
		}
		else
		{
			//check if its a genmergerspanning vis arc, merge cases later
			node firstn = nullptr, secondn = nullptr;
			for (node n : m_path[(*it).x1()])
			{
				node en = m_pPR->expandedNode(n);
				if (!en) continue;
				if (!(m_pPR->typeOf(n) == Graph::NodeType::generalizationExpander)) continue;
				else { firstn = en; break; }
			}//for
			for (node n : m_path[(*it).x2()])
			{
				node en = m_pPR->expandedNode(n);
				if (!en) continue;
				if (!(m_pPR->typeOf(n) == Graph::NodeType::generalizationExpander)) continue;
				else { secondn = en; break; }
			}//for
			if ((firstn && secondn) && (firstn == secondn))
			{
				if (itPrev.valid()) visibArcs.delSucc(itPrev);
				else visibArcs.popFront();
			}
			else itPrev = it;
		}
	}//for visibArcs

}



// output in gml-format with special edge colouring
// arcs with cost 0 are green, other arcs red
void CompactionConstraintGraphBase::writeGML(const char *filename) const
{
	ofstream os(filename);
	writeGML(os);
}
void CompactionConstraintGraphBase::writeGML(const char *filename, NodeArray<bool> one) const
{
	ofstream os(filename);
	writeGML(os, one);
}

void CompactionConstraintGraphBase::writeGML(ostream &os) const
{
	const Graph &G = *this;

	NodeArray<int> id(*this);
	int nextId = 0;

	os.setf(ios::showpoint);
	os.precision(10);

	os << "Creator \"ogdf::CompactionConstraintGraphBase::writeGML\"\n";
	os << "graph [\n";
	os << "  directed 1\n";

	for(node v : G.nodes) {
		os << "  node [\n";

		os << "    id " << (id[v] = nextId++) << "\n";

		os << "    graphics [\n";
		os << "      x 0.0\n";
		os << "      y 0.0\n";
		os << "      w 30.0\n";
		os << "      h 30.0\n";

		os << "      fill \"#FFFF00\"\n";
		os << "    ]\n"; // graphics

		os << "  ]\n"; // node
	}


	for(edge e : G.edges) {
		os << "  edge [\n";

		os << "    source " << id[e->source()] << "\n";
		os << "    target " << id[e->target()] << "\n";

#if 0
		// show edge lengths as edge lables (not yet supported)
		os << "    label \"";
		writeLength(os,e);
		os << "\"\n";
#endif

		os << "    graphics [\n";

		os << "      type \"line\"\n";
		os << "      arrow \"last\"\n";
		switch(m_type[e])
		{
		case ConstraintEdgeType::BasicArc: // red
			os << "      fill \"#FF0000\"\n";
			break;
		case ConstraintEdgeType::VertexSizeArc: // blue
			os << "      fill \"#0000FF\"\n";
			break;
		case ConstraintEdgeType::VisibilityArc: // green
			os << "      fill \"#00FF00\"\n";
			break;
		case ConstraintEdgeType::ReducibleArc: // rose
			os << "      fill \"#FF00FF\"\n";
			break;
		case ConstraintEdgeType::FixToZeroArc: //violett
			os << "      fill \"#3F00FF\"\n";
			break;
		default:
			OGDF_ASSERT(false);
		}

		os << "    ]\n"; // graphics

#if 0
		os << "    LabelGraphics [\n";
		os << "      type \"text\"\n";
		os << "      fill \"#000000\"\n";
		os << "      anchor \"w\"\n";
		os << "    ]\n";
#endif

		os << "  ]\n"; // edge
	}

	os << "]\n"; // graph
}//writegml

void CompactionConstraintGraphBase::writeGML(ostream &os, NodeArray<bool> one) const
{
	const Graph &G = *this;

	NodeArray<int> id(*this);
	int nextId = 0;

	os.setf(ios::showpoint);
	os.precision(10);

	os << "Creator \"ogdf::CompactionConstraintGraphBase::writeGML\"\n";
	os << "graph [\n";
	os << "  directed 1\n";

	for(node v : G.nodes) {
		os << "  node [\n";

		os << "    id " << (id[v] = nextId++) << "\n";

		os << "    graphics [\n";
		os << "      x 0.0\n";
		os << "      y 0.0\n";
		os << "      w 30.0\n";
		os << "      h 30.0\n";
		if (one[v]) {
			os << "      fill \"#FF0F0F\"\n";
		} else {
			os << "      fill \"#FFFF00\"\n";
		}
		os << "    ]\n"; // graphics

		os << "  ]\n"; // node
	}


	for(edge e : G.edges) {
		os << "  edge [\n";

		os << "    source " << id[e->source()] << "\n";
		os << "    target " << id[e->target()] << "\n";

#if 0
		// show edge lengths as edge lables (not yet supported)
		os << "    label \"";
		writeLength(os,e);
		os << "\"\n";
#endif

		os << "    graphics [\n";

		os << "      type \"line\"\n";
		os << "      arrow \"last\"\n";
		switch(m_type[e])
		{
		case ConstraintEdgeType::BasicArc: // red
			os << "      fill \"#FF0000\"\n";
			break;
		case ConstraintEdgeType::VertexSizeArc: // blue
			os << "      fill \"#0000FF\"\n";
			break;
		case ConstraintEdgeType::VisibilityArc: // green
			os <<       "fill \"#00FF00\"\n";
			break;
		case ConstraintEdgeType::ReducibleArc: // rose
			os << "      fill \"#FF00FF\"\n";
			break;
		case ConstraintEdgeType::FixToZeroArc: //violett
			os << "      fill \"#3F00FF\"\n";
			break;
		default:
			OGDF_ASSERT(false);
		}

		os << "    ]\n"; // graphics

#if 0
		os << "    LabelGraphics [\n";
		os << "      type \"text\"\n";
		os << "      fill \"#000000\"\n";
		os << "      anchor \"w\"\n";
		os << "    ]\n";
#endif

		os << "  ]\n"; // edge
	}

	os << "]\n"; // graph
}


} // end namespace ogdf
