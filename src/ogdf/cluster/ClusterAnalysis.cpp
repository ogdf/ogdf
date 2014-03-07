#include <ogdf/cluster/ClusterAnalysis.h>
#include <ogdf/cluster/ClusterArray.h>
#include <ogdf/basic/Queue.h>
#include <ogdf/basic/DisjointSets.h>
#include <ogdf/basic/exceptions.h>

#include <ogdf/basic/Logger.h>

// Comment on use of ClusterArrays:
// We would like to save some space by only reserving one slot
// per existing cluster instead of maxClusterIndex() slots,
// which might be larger. However, we then would
// need to store an index with each cluster in a struct here,
// minimizing the effect again. Note ClusterArrays are static
// here, a change in the ClusterGraph won't be detected
// by ClusterAnalysis after initialization.
// Where multiple indexing is done, the static version
// of ClusterArrays is used with successive cluster index
// numbers computed locally, to save memory. When only a single
// indexed structure is needed this is clearly an overhead.
// Therefore the ClusterArrays are then created with size maxClusterIndex().


namespace ogdf {

//Needs to be the largest int allowed, as it is used as default,
//and an update is done for smaller values
const int ClusterAnalysis::IsNotActiveBound = numeric_limits<int>::max();
const int ClusterAnalysis::DefaultIndex = -1;

//Constructor
ClusterAnalysis::ClusterAnalysis(const ClusterGraph& C, bool oalists, bool indyBags) : m_C(&C), m_oanum(0), m_ianum(0),
		m_bags(0), m_storeoalists(oalists), m_lcaEdges(0), m_indyBags(indyBags), m_numIndyBags(-1), m_indyBagRoots(0)
{
	init();
	computeBags();
	if (m_indyBags)
	{
		computeIndyBags();
	}
}
ClusterAnalysis::ClusterAnalysis(const ClusterGraph& C, bool indyBags): m_C(&C), m_oanum(0), m_ianum(0),
			m_bags(0), m_storeoalists(true), m_lcaEdges(0), m_indyBags(indyBags),m_numIndyBags(-1), m_indyBagRoots(0)
{
	init();
	computeBags();
	// Even though it looks like we could compute bags and
	// indyBags in one pass, we do not look at each vertex
	// in each cluster during bag computation, as this would
	// be inefficient (use Union-Find instead). However, for
	// independent bags, we need to identify bags wo outeractive
	// vertices per cluster (TODO there may be shortcuts however, e.g.
	// if number of outeractive vertices is zero for a cluster, then
	// every bag is independent).
	if (m_indyBags)
	{
		computeIndyBags();
	}
}
ClusterAnalysis::~ClusterAnalysis()
{
	cleanUp();
}

void ClusterAnalysis::cleanUp()
{
	if (m_oanum != 0) delete m_oanum;
	if (m_ianum != 0) delete m_ianum;
	if (m_bags != 0) delete m_bags;
	if (m_lcaEdges != 0) delete m_lcaEdges;
	if (m_storeoalists) delete m_oalists;
	node v;
	forall_nodes(v,m_C->constGraph())
	{
		delete m_bagindex[v];
	}
	if (m_indyBags) delete m_indyBagRoots;
}

int ClusterAnalysis::outerActive(cluster c) {return (*m_oanum)[c];}
int ClusterAnalysis::innerActive(cluster c) {return (*m_ianum)[c];}
int ClusterAnalysis::numberOfBags(cluster c) {return (*m_bags)[c];}
List<node>& ClusterAnalysis::oaNodes(cluster c) {
	OGDF_ASSERT(m_storeoalists);
	return (*m_oalists)[c];
}

bool ClusterAnalysis::isOuterActive(node v, cluster c)
{
	return (*m_oactive[v])[c] > 0;
}
bool ClusterAnalysis::isInnerActive(node v, cluster c)
{
	return (*m_iactive[v])[c] > 0;
}
List<edge>& ClusterAnalysis::lcaEdges(cluster c)
{
	return (*m_lcaEdges)[c];
}
int ClusterAnalysis::bagIndex(node v, cluster c) {return ((*m_bagindex[v])[c]);}

int ClusterAnalysis::indyBagIndex(node v) {
	if (!m_indyBags)
	{
		OGDF_THROW_PARAM(AlgorithmFailureException,afcIllegalParameter);
	}
	return m_indyBagNumber[v];

}

cluster ClusterAnalysis::indyBagRoot(int i) {
	if (!m_indyBags)
	{
		OGDF_THROW_PARAM(AlgorithmFailureException,afcIllegalParameter);
	}
	return m_indyBagRoots[i];
}
//we fill all arrays that store the inner/outer activity status
void ClusterAnalysis::init() {
//does it make sense to split into detecting inner / outer activity?
//is it faster to run once over the edges, with bidirectional detection
//and a check if the edge is processed already, or is it at least almost
//always as fast when testing from both directions?

	//first outactive vertices
	const Graph &G = m_C->constGraph();
	m_iactive.init(G);
	m_oactive.init(G);
	m_ialevel.init(G, IsNotActiveBound);
	m_oalevel.init(G, IsNotActiveBound);

	if (m_oanum != 0) delete m_oanum;
	if (m_ianum != 0) delete m_ianum;
	if (m_bags != 0) delete m_bags;
	m_oanum = new ClusterArray<int>(*m_C, 0);
	m_ianum = new ClusterArray<int>(*m_C, 0);
	m_bags = new ClusterArray<int>(*m_C, 0);
	m_lcaEdges = new ClusterArray<List<edge> >(*m_C);
	if (m_storeoalists)
		m_oalists = new ClusterArray<List<node> >(*m_C);

	//We don't want to set dynamic depths update for clusters in m_C,
	//therefore we just compute the values here
	//top-down run through the cluster tree, depth 0 for the root
	ClusterArray<int> cdepth(*m_C);
	cluster c = m_C->rootCluster();
	cdepth[c] = 0;
	Queue<cluster> cq;
	ListConstIterator<cluster> it = c->cBegin();
	while (it.valid())
	{
		cq.append(*it);
		it++;
	}

	while (!cq.empty())
	{
		cluster cc = cq.pop();
		cdepth[cc] = cdepth[cc->parent()]+1;
		it = cc->cBegin();
		while (it.valid())
		{
			cq.append(*it);
			it++;
		}
	}

	node v;
	edge e;

	//store that we already visited e, as we don't have a static lookup
	//for the paths, running the search from both directions is slower.
	EdgeArray<bool> visited(G, false);
	forall_nodes(v, G)
	{
		// See comment on use of ClusterArrays above
		m_iactive[v] = new ClusterArray<int>(*m_C,0,m_C->maxClusterIndex()+1);
		m_oactive[v] = new ClusterArray<int>(*m_C,0,m_C->maxClusterIndex()+1);
	}
	forall_nodes(v, G)
	{
		forall_adj_edges(e,v)
		{
			if (!visited[e])
			{
				node w = e->opposite(v);
				List<cluster> el; // result cluster list of path in cluster tree T between v,w
				cluster c1,c2;   // ancestors of lca(v,w) on path, we don't really need them here

				cluster lca = m_C->commonClusterAncestorsPath(v, w, c1, c2, el);
				OGDF_ASSERT(el.size() > 0)
				ListIterator<cluster> ctit =  el.begin();//m_C->clusterOf(v);

				//run over the path, set activity status (vertices are
				//active for a cluster if adjacent edge crosses the border

				//clusters before lca are left, i.e. v is outer active
				//clusters behind lca are entered, i.e. v is inner active
				while (ctit.valid() && (*ctit) != lca)
				{
					(*m_oactive[v])[(*ctit)]++;
					(*m_iactive[w])[(*ctit)]++;

					//only count vertices a single time
					if ((*m_oactive[v])[(*ctit)] == 1) (*m_oanum)[(*ctit)]++;
					if ((*m_iactive[w])[(*ctit)] == 1) (*m_ianum)[(*ctit)]++;

					//update the activity levels
					//could do this just for the last in the line...
					int clevel = cdepth[(*ctit)];
					if (m_ialevel[w] > clevel)
						m_ialevel[w] = clevel;
					if (m_oalevel[v] > clevel)
						m_oalevel[v] = clevel;

					ctit++;
				}

				OGDF_ASSERT((*ctit) == lca);
				//vertices are never active wrt lca
				//we store however the corresponding edges
				//for later use in bag detection
				(*m_lcaEdges)[lca].pushBack(e);
				ctit++;

				while (ctit.valid())
				{
					(*m_iactive[v])[(*ctit)]++;
					(*m_oactive[w])[(*ctit)]++;

					if ((*m_iactive[v])[(*ctit)] == 1) (*m_ianum)[(*ctit)]++;
					if ((*m_oactive[w])[(*ctit)] == 1) (*m_oanum)[(*ctit)]++;

					//update the activity levels
					int clevel = cdepth[(*ctit)];
					//could do this just for the last in the line...
					if (m_ialevel[v] > clevel)
						m_ialevel[v] = clevel;
					if (m_oalevel[w] > clevel)
						m_oalevel[w] = clevel;

					ctit++;
				}

				visited[e] = true;

#ifdef OGDF_DEBUG
				cout << "Edge "<< v << " " << w <<"\n";
#endif
			}
		}

	}
#ifdef OGDF_DEBUG
	forall_nodes(v, G)
	{
		cout << "Knoten "<<v<<" ist";
		List<cluster> ol;
		List<cluster> il;
		cluster c;
		forall_clusters(c, *m_C)
		{
			if ((*m_iactive[v])[c]>0) il.pushBack(c);
			if ((*m_oactive[v])[c]>0) ol.pushBack(c);
		}
		cout << " inneractive for ";
		ListIterator<cluster> it = il.begin();
		while (it.valid())
		{
			cout << (*it) << ", ";
			it++;
		}
		cout << "\n";
		cout << " outeractive for ";
		it = ol.begin();
		while (it.valid())
		{
			cout << (*it) << ", ";
			it++;
		}
		cout << "\n";
	}

#endif
}

// Runs through a list of vertices (starting with \a the one nodeIT points to)
// which is expected to be a full list of cluster vertices in \a c. Depending on
// outer activity and bag index number of the vertices, independent bags
// are detected and a corresponding index is assigned accordingly for each vertex.
void ClusterAnalysis::partitionCluster(ListConstIterator<node> & nodeIt, cluster c,
		HashArray<int, List<node> > & bagNodes, HashArray<int, bool> & indyBag,
		Skiplist<int*> & indexNumbers, Array<cluster> & bagRoots) {
	// Run through all vertices in c
	while (nodeIt.valid())
	{
		node v = *nodeIt;
		//if vertex is outeractive, the containing bag loses its status
		// Nodes that are already processed can never be outeractive
		// as we traverse bottom up. They are skipped as they are part
		// of an independent bag in an already processed cluster.
		if (m_indyBagNumber[v] == DefaultIndex)
		{
			int ind = bagIndex(v,c);
			if (isOuterActive(v,c))
			{
				indyBag[ind]= false;
			} else //don't need to add the index if vertex is oactive
			{
				if (!indexNumbers.isElement(&ind))
				{
					indexNumbers.add(new int(ind));
				}
				bagNodes[ind].pushBack(v); //store vertex in index' list
			}
		}
#ifdef OGDF_DEBUG
		if (m_indyBagNumber[v] != DefaultIndex)
			OGDF_ASSERT(!isOuterActive(v,c));
#endif

		nodeIt++;
	}
	// Now we have all indexes of bags that don't solely contain oactive vertices.
	// For each index we check if the bag still has independency status,
	// in this case we have found an independent bag and can remove all its
	// vertices (mark them).

	//cout << "Checking skiplist\n";
	SkiplistIterator<int*> its = indexNumbers.begin();
	while (its.valid())
	{
		int bind = *(*its);
		//cout << "Found index "<< bind <<"\n";
		if (indyBag[bind])
		{
#ifdef OGDF_DEBUG
			Logger::slout()  << "Found independent bag with "<< bagNodes[bind].size() << "vertices\n";
#endif

			ListConstIterator<node> itn = bagNodes[bind].begin();
			while (itn.valid())
			{
				// Assign the final index number
				m_indyBagNumber[*itn] = m_numIndyBags;
				itn++;
			}
			bagRoots[m_numIndyBags] = c;
			m_numIndyBags++;
		}

		delete (*its);
		its++;
	}
}//partitionCluster
// For each cluster we check if we can identify an independent
// bag, which might be useful for clustered planarity testing.
// compute independent bag affiliation for all vertices,
// store result in m_indyBagNumber, and set m_numIndyBags accordingly.
// We could interweave this with computeBags to be more efficient but
// this would it make the code more complex and in case we don't need
// the IndyBags we would do additional work in vain.
void ClusterAnalysis::computeIndyBags() {
	m_numIndyBags = 0; //used both to count the bags and to store current
					   //indyBag index number for vertex assignment (i.e., starts with 1)
	const Graph &G = m_C->constGraph();

	// Store the root cluster of each indyBag
	if (m_indyBagRoots) delete m_indyBagRoots;
	// Intermediate storage during computation, maximum of #vertices possible
#ifdef OGDF_DEBUG
	Array<cluster> bagRoots(0,G.numberOfNodes(),0);
#else
	Array<cluster> bagRoots(G.numberOfNodes());
#endif

	// Store indyBag affiliation. Every vertex will get a number != -1 (DefaultIndex,
	// as in the worst case the whole graph is an indyBag (in root cluster).
	// Once assigned, the number won't change during the processing.
	m_indyBagNumber.init(G, DefaultIndex);

	// We run bottom up over all clusters (to find the minimum inclusion).
	// For each vertex, we use the outer activity and bag index information.
	// In case we find a bag without outeractive vertices it is a IndyBag.
	// Already processed vertices are simply marked by an indyBag index entry different to -1.
	//NodeArray<bool> processed(G, false); //mark vertices when IndyBag detected

	// We do not have the sets of vertices for all bags, as only the bag index
	// has been stored for a cluster.
	// Detect the current leaf clusters for bottom up traversal.
	List<cluster> ccleafs;
	ClusterArray<int> unprocessedChildren(*m_C); //processing below: compute bags
	cluster c;
	forall_clusters(c, *m_C)
	{
		if (c->cCount() == 0) ccleafs.pushBack(c);
		unprocessedChildren[c] = c->cCount();
	}
	OGDF_ASSERT(!ccleafs.empty());

	// Run through all clusters, leaves first.
	while (!ccleafs.empty()){
		// We cannot store the following information over the whole
		// graph even though the bag index (which is one out of the
		// set of vertex set ids in union find) would allow this.
		// The bag index is defined per cluster.
		// However, when moving up in the cluster tree indy information
		// is not monotone, we therefore would need to update it accordingly.
		// (Bag which is not outeractive will not become outeractive, but
		// may get a part of an outeractive bag with the same id, and an
		// outeractive bag might become enclosed).
		HashArray<int, bool> indyBag(true); //true if bag with index i does
											//not have outeractive vertices
		// We want to store all vertices for each index that may be a
		// potential indyBag index. We could add these in the entry stored
		// in our index Skiplist, but then we need a comparison of the
		// respecting class objects, slowing down the processing.
		HashArray<int, List<node> > bagNodes; //holds vertices for index i

		// We need a data structure that holds all indexes used, allows
		// to search for them and to iterate through them. Could use
		// insertion sorted dynamically allocated array here, but as
		// SkipLists are implemented in OGDF, which provide comparable
		// efficiency, we use them.
		Skiplist<int*> indexNumbers; //holds all index numbers in c

		cluster c = ccleafs.popFrontRet();

		List<node> nodes;
		ListConstIterator<node> it;
		// Process leaves separately. As long as we don't do something
		// special for non-leaves, we could just use getClusterNodes instead.
		if (c->cCount() == 0) {
			it = c->nBegin();
		}
		else {
			// At this point all child clusters of c have been processed.
			// This is somewhat slow, as we repeatedly collect the vertices
			// Todo: Instead, we could clone the clusternodes code
			// here and process them as we see them.

			c->getClusterNodes(nodes);

			//cout <<"Processing cluster with "<<nodes.size()<<"\n";

			it = nodes.begin();
		}
		// Run through all vertices in c
		partitionCluster(it, c, bagNodes, indyBag, indexNumbers, bagRoots);
		// Now we update the status of the parent cluster and,
		// in case all its children are processed, add it to
		// the process queue.
		if (c != m_C->rootCluster())
		{
			OGDF_ASSERT(unprocessedChildren[c->parent()] > 0)
			unprocessedChildren[c->parent()]--;
			if (unprocessedChildren[c->parent()] == 0) ccleafs.pushBack(c->parent());
		}
	}//while clusters
	// Copying into smaller array is a bit slower than just reserving the whole #v array,
	// but we do it anyway.
	m_indyBagRoots = new cluster[m_numIndyBags];
	for (int k = 0; k < m_numIndyBags; k++)
	{
		OGDF_ASSERT(bagRoots[k] != 0);
		m_indyBagRoots[k] = bagRoots[k];
	}
#ifdef OGDF_DEBUG
/*
	List<node> rnodes;
	m_C->rootCluster()->getClusterNodes(rnodes);

	ListConstIterator<node> rit = rnodes.begin();

	while (rit.valid())
	{
		cout << "Root bag index: "<<bagIndex(*rit, m_C->rootCluster())<<"\n";
		cout << "Indy bag index: "<<m_indyBagNumber[*rit]<<"\n";
		rit++;
	}
*/

	Skiplist<int*> ibind;
	node v;
	forall_nodes(v, G)
	{
		int i = m_indyBagNumber[v];
		if (!ibind.isElement(&i))
		{
			ibind.add(new int(i));
		}
		cout << "numIndyBags: "<<m_numIndyBags<<" i: "<<i<<"\n";
		OGDF_ASSERT(i != DefaultIndex);
		OGDF_ASSERT((i >= 0) && (i < m_numIndyBags));

	}
	int storedBags = 0;
	SkiplistIterator<int*> its = ibind.begin();
	while (its.valid())
	{
		storedBags++;
		delete (*its);
		its++;
	}
	Logger::slout() << m_numIndyBags<< " independent bags detected, "<<storedBags<<" stored\n";
	OGDF_ASSERT(m_numIndyBags==storedBags);
#endif
}

//compute bag affiliation for all vertices
//store result in m_bagindex
void ClusterAnalysis::computeBags() {
	const Graph &G = m_C->constGraph();

	// Storage structure for results
	m_bagindex.init(G);
	// We use Union-Find for chunks and bags
	DisjointSets<> uf;
	NodeArray<int> setid(G); // Index mapping for union-find
	//node* nn = new node[G.numberOfNodes()]; // dito

	node v;

	// Every cluster gets its index
	ClusterArray<int> cind(*m_C);
	// We store the lists of cluster vertices
	List<node>* clists = new List<node>[m_C->numberOfClusters()];
	cluster c;
	int i = 0;

	// Store index and detect the current leaf clusters
	List<cluster> ccleafs;
	ClusterArray<int> unprocessedChildren(*m_C); //processing below: compute bags
	forall_clusters(c, *m_C)
	{
		cind[c] = i++;
		if (c->cCount() == 0) ccleafs.pushBack(c);
		unprocessedChildren[c] = c->cCount();
	}


	// Now we run through all vertices, storing them in the parent lists,
	// at the same time, we initialize m_bagindex
	forall_nodes(v, G)
	{
		// setid is constant in the following
		setid[v] = uf.makeSet();
		// Each vertex v gets its own ClusterArray that stores v's bag index per cluster.
		// See comment on use of ClusterArrays above
		m_bagindex[v] = new ClusterArray<int>(*m_C,DefaultIndex, m_C->maxClusterIndex()+1);//m_C->numberOfClusters());
		cluster c = m_C->clusterOf(v);
		// Push vertices in parent list
		clists[cind[c]].pushBack(v);
	}

	// Now each clist contains the direct vertex descendants
	// We process the clusters bottom-up, compute the chunks
	// of the leafs first. At each level, for a cluster the
	// vertex lists of all children are concatenated
	// (could improve this by having an array of size(#leafs)
	// and concatenating only at child1), then the bags are
	// updated as follows: chunks may be linked by exactly
	// the edges with lca(c) ie the ones in m_lcaEdges[c],
	// and bags may be built by direct child clusters that join chunks.
	// While concatenating the vertex lists, we can check
	// for the vertices in each child if the uf number is the same
	// as the one of a first initial vertex, otherwise we join.

	// First, lowest level clusters are processed: All chunks are bags


	OGDF_ASSERT(!ccleafs.empty());

	while (!ccleafs.empty()){
		cluster c = ccleafs.popFrontRet();
		Skiplist<int* > cbags; //Stores bag indexes ocurring in c
		if (m_storeoalists){
			//no outeractive vertices detected so far
			(*m_oalists)[c].clear();
		}

		//process leafs separately
		if (c->cCount() == 0) {


			//Todo could use lcaEdges list here too, see below
			ListConstIterator<node> it = c->nBegin();
			while (it.valid())
			{
				edge e;
				forall_adj_edges(e,*it)
				{
					node w = e->opposite(*it);
					if (m_C->clusterOf(w) == c)
					{
						uf.link(uf.find(setid[*it]),uf.find(setid[w]));
					}
				}

				it++;
			}
			// Now all chunks in the leaf cluster are computed
			// update for parent is done in the else case
			// We store the result:
			ListConstIterator<node> itv = clists[cind[c]].begin();
			while (itv.valid())
			{
				int theid = uf.find(setid[*itv]);
				(*m_bagindex[*itv])[c] = theid;
				if (!cbags.isElement(&theid))
				{
					cbags.add(new int(theid));
				}
				//push into list of outer active vertices
				if (m_storeoalists){
					if (isOuterActive(*itv, c)) ((*m_oalists)[c]).pushBack(*itv);
				}
				itv++;
			}
			//getNumberOfSets would be all sets, we only need the ones in the cluster
			(*m_bags)[c] = cbags.size();
		}
		else {
			// ?We construct the vertex list by concatenating
			// ?the lists of the children to the current list.
			// We need the lists for storing the results efficiently.
			// (Should be slightly faster than to call clusterNodes each time)
			// Bags are either links of chunks by edges with lca==c
			// or links of chunk by child clusters.
			// Edge links
			ListConstIterator<edge> ite = (*m_lcaEdges)[c].begin();
			while (ite.valid())
			{
				uf.link(uf.find(setid[(*ite)->source()]),uf.find(setid[(*ite)->target()]));
				ite++;
			}
			// Cluster links
			ListConstIterator<cluster> itcc = c->cBegin();

			while (itcc.valid())
			{
				//Initial id per child cluster cc: Use value of first
				//vertex, each time we encounter a different value in cc,
				//we link the chunks

				//add (*itcc)'s vertices to c's list
				ListConstIterator<node> itvc = clists[cind[(*itcc)]].begin();
				int inid;
				if (itvc.valid()) inid = uf.find(setid[(*itvc)]);
				while (itvc.valid())
				{
					int theid = uf.find(setid[(*itvc)]);

					if (theid != inid)
						uf.link(inid,theid);
					clists[cind[c]].pushBack(*itvc);
					itvc++;
				}
				itcc++;

			}
			// We store the result:
			ListConstIterator<node> itv = clists[cind[c]].begin();
			while (itv.valid())
			{
				int theid = uf.find(setid[*itv]);
				(*m_bagindex[*itv])[c] = theid;
				if (!cbags.isElement(&theid))
				{
					cbags.add(new int(theid));
				}
				//push into list of outer active vertices
				if (m_storeoalists){
					if (isOuterActive(*itv, c)) ((*m_oalists)[c]).pushBack(*itv);
				}
				itv++;
			}
			(*m_bags)[c] = cbags.size(); //store number of bags of c
		}
		// Now we update the status of the parent cluster and,
		// in case all its children are processed, add it to
		// the process queue.
		if (c != m_C->rootCluster())
		{
			OGDF_ASSERT(unprocessedChildren[c->parent()] > 0)
			unprocessedChildren[c->parent()]--;
			if (unprocessedChildren[c->parent()] == 0) ccleafs.pushBack(c->parent());
		}
	}//while cluster

	// clean up
	delete [] clists;
}


}
