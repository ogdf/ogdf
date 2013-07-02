//*********************************************************
//  Regression test for planarity tests and embeddings
//
//  Tested classes:
//    - BoothLueker
//    - BoyerMyrvold
//
//  Author: Carsten Gutwenger
//*********************************************************

#include <ogdf/basic/graph_generators.h>
#include <ogdf/planarity/BoothLueker.h>
#include <ogdf/planarity/BoyerMyrvold.h>
#include <ogdf/planarity/FastPlanarSubgraph.h>
#include <ogdf/planarity/MaximalPlanarSubgraphSimple.h>

using namespace ogdf;


void removeEdges(Graph &G, const List<edge> &delEdges)
{
	ListConstIterator<edge> it;
	for(it = delEdges.begin(); it.valid(); ++it)
		G.delEdge(*it);
}

void randomizeAdjLists(Graph G)
{
	node v;
	forall_nodes(v,G) {
		List<adjEntry> L;
		G.adjEntries(v,L);
		L.permute();
		G.sort(v,L);
	}
}

void addRandomEdges(Graph &G, int m)
{
	const int n = G.numberOfNodes();

	Array<node> nodes(n);
	node v;
	int i = 0;
	forall_nodes(v,G)
		nodes[i++] = v;

	for(i = 0; i < m; ++i) {
		int i = randomNumber(0,n-1);
		int j = randomNumber(1,n-1);
		G.newEdge(nodes[i],nodes[(i+j)%n]);
	}
}

void addRandomLoops(Graph &G, int m)
{
	const int n = G.numberOfNodes();

	Array<node> nodes(n);
	node v;
	int i = 0;
	forall_nodes(v,G)
		nodes[i++] = v;

	for(i = 0; i < m; ++i) {
		node v = nodes[randomNumber(0,n-1)];
		G.newEdge(v,v);
	}
}

void addRandomMultiEdges(Graph &G, int add)
{
	const int m = G.numberOfEdges();

	Array<edge> edges(m);
	edge e;
	int i = 0;
	forall_edges(e,G)
		edges[i++] = e;

	for(i = 0; i < add; ++i) {
		e = edges[randomNumber(0,m-1)];
		G.newEdge(e->source(), e->target());
	}
}

bool regPlanarityTest()
{
	const int numGraphs = 10;
	const int min_n = 500;
	const int max_n = 2500;
	const int step_n = 50;

	bool result = true;

	Graph G1, G2;
	BoothLueker bl;
	BoyerMyrvold bm;
	MaximalPlanarSubgraphSimple mps;

	cout << "-> Planar biconnected graphs, planarity test... " << endl;
	srand(4711);

	__int64 T;
	__int64 msecBL = 0, msecBM = 0;
	int failsBL = 0, failsBM = 0;

	bool firstBL = true;
	for(int n = min_n; n <= max_n; n += step_n)
	{
		cout << "\r" << n << flush;
		for(int i = 0; i < numGraphs; ++i)
		{
			bool resPM, resBM;
			__int64 tBL, tBM;

			int m = 3*n/2;
			planarBiconnectedGraph(G1, n, m);
			addRandomLoops(G1,10);
			randomizeAdjLists(G1);
			G2 = G1;

			System::usedRealTime(T);
			if(firstBL) {
				resPM = bl.isPlanar(G1);
				tBL = System::usedRealTime(T);

				resBM = bm.isPlanar(G2);
				tBM = System::usedRealTime(T);
			} else {
				resBM = bm.isPlanar(G1);
				tBM = System::usedRealTime(T);

				resPM = bl.isPlanar(G2);
				tBL = System::usedRealTime(T);
			}
			firstBL = !firstBL;

			msecBL += tBL;
			msecBM += tBM;
			if(resPM == false) failsBL++;
			if(resBM == false) failsBM++;

			m = 2*n;
			planarBiconnectedGraph(G1, n, m);
			addRandomLoops(G1,10);
			randomizeAdjLists(G1);
			G2 = G1;

			System::usedRealTime(T);
			if(firstBL) {
				resPM = bl.isPlanar(G1);
				tBL = System::usedRealTime(T);

				resBM = bm.isPlanar(G2);
				tBM = System::usedRealTime(T);
			} else {
				resBM = bm.isPlanar(G1);
				tBM = System::usedRealTime(T);

				resPM = bl.isPlanar(G2);
				tBL = System::usedRealTime(T);
			}
			firstBL = !firstBL;

			msecBL += tBL;
			msecBM += tBM;
			if(resPM == false) failsBL++;
			if(resBM == false) failsBM++;

			m = 5*n/2;
			planarBiconnectedGraph(G1, n, m);
			addRandomLoops(G1,10);
			randomizeAdjLists(G1);
			G2 = G1;

			System::usedRealTime(T);
			if(firstBL) {
				resPM = bl.isPlanar(G1);
				tBL = System::usedRealTime(T);

				resBM = bm.isPlanar(G2);
				tBM = System::usedRealTime(T);
			} else {
				resBM = bm.isPlanar(G1);
				tBM = System::usedRealTime(T);

				resPM = bl.isPlanar(G2);
				tBL = System::usedRealTime(T);
			}
			firstBL = !firstBL;

			msecBL += tBL;
			msecBM += tBM;
			if(resPM == false) failsBL++;
			if(resBM == false) failsBM++;
		}
	}

	if(failsBL + failsBM > 0) result = false;

	cout.precision(3);
	cout << "\r" << std::fixed;
	cout << "  BoothLueker:  " << 0.001*double(msecBL) << " seconds (" << failsBL << " errors)" << endl;
	cout << "  BoyerMyrvold: " << 0.001*double(msecBM) << " seconds (" << failsBM << " errors)" << endl;

	cout << "-> Connected graphs, planarity test... " << endl;
	srand(4711);

	msecBL  = msecBM  = 0;
	failsBL = failsBM = 0;

	firstBL = true;
	for(int n = min_n; n <= max_n; n += step_n)
	{
		cout << "\r" << n << flush;
		for(int i = 0; i < numGraphs; ++i)
		{
			bool resPM, resBM;
			__int64 tBL, tBM;
			int nb = n/10;

			planarCNBGraph(G1, nb, 3*nb/2, 10);
			addRandomLoops(G1,10);
			randomizeAdjLists(G1);
			G2 = G1;

			System::usedRealTime(T);
			if(firstBL) {
				resPM = bl.isPlanar(G1);
				tBL = System::usedRealTime(T);

				resBM = bm.isPlanar(G2);
				tBM = System::usedRealTime(T);
			} else {
				resBM = bm.isPlanar(G1);
				tBM = System::usedRealTime(T);

				resPM = bl.isPlanar(G2);
				tBL = System::usedRealTime(T);
			}
			firstBL = !firstBL;

			msecBL += tBL;
			msecBM += tBM;
			if(resPM == false) failsBL++;
			if(resBM == false) failsBM++;

			planarCNBGraph(G1, nb, 2*nb, 10);
			addRandomLoops(G1,10);
			randomizeAdjLists(G1);
			G2 = G1;

			System::usedRealTime(T);
			if(firstBL) {
				resPM = bl.isPlanar(G1);
				tBL = System::usedRealTime(T);

				resBM = bm.isPlanar(G2);
				tBM = System::usedRealTime(T);
			} else {
				resBM = bm.isPlanar(G1);
				tBM = System::usedRealTime(T);

				resPM = bl.isPlanar(G2);
				tBL = System::usedRealTime(T);
			}
			firstBL = !firstBL;

			msecBL += tBL;
			msecBM += tBM;
			if(resPM == false) failsBL++;
			if(resBM == false) failsBM++;

			planarCNBGraph(G1, nb, 5*nb/2, 10);
			addRandomLoops(G1,10);
			randomizeAdjLists(G1);
			G2 = G1;

			System::usedRealTime(T);
			if(firstBL) {
				resPM = bl.isPlanar(G1);
				tBL = System::usedRealTime(T);

				resBM = bm.isPlanar(G2);
				tBM = System::usedRealTime(T);
			} else {
				resBM = bm.isPlanar(G1);
				tBM = System::usedRealTime(T);

				resPM = bl.isPlanar(G2);
				tBL = System::usedRealTime(T);
			}
			firstBL = !firstBL;

			msecBL += tBL;
			msecBM += tBM;
			if(resPM == false) failsBL++;
			if(resBM == false) failsBM++;
		}
	}

	if(failsBL + failsBM > 0) result = false;

	cout.precision(3);
	cout << "\r";
	cout << "  BoothLueker:  " << 0.001*double(msecBL) << " seconds (" << failsBL << " errors)" << endl;
	cout << "  BoyerMyrvold: " << 0.001*double(msecBM) << " seconds (" << failsBM << " errors)" << endl;

	cout << "-> Planar biconnected graphs, planar embedding... " << endl;
	srand(4711);

	msecBL  = msecBM  = 0;
	failsBL = failsBM = 0;

	firstBL = true;
	for(int n = min_n; n <= max_n; n += step_n)
	{
		SList< KuratowskiWrapper > dummyList;
		cout << "\r" << n << flush;
		for(int i = 0; i < numGraphs; ++i)
		{
			//bool resPM, resBM;
			__int64 tBL, tBM;

			int m = 3*n/2+1;
			planarBiconnectedGraph(G1, n, m);
			addRandomLoops(G1,10);
			addRandomMultiEdges(G1, 50);
			randomizeAdjLists(G1);
			G2 = G1;

			System::usedRealTime(T);
			if(firstBL) {
				bl.planarEmbed(G1);
				tBL = System::usedRealTime(T);

				bm.planarEmbed(G2,dummyList);
				tBM = System::usedRealTime(T);
			} else {
				bm.planarEmbed(G1,dummyList);
				tBM = System::usedRealTime(T);

				bl.planarEmbed(G2);
				tBL = System::usedRealTime(T);
			}

			if(!G1.representsCombEmbedding()) {
				if(firstBL) failsBL++; else failsBM++;
			}
			if(!G2.representsCombEmbedding()) {
				if(!firstBL) failsBL++; else failsBM++;
			}

			msecBL += tBL;
			msecBM += tBM;
			firstBL = !firstBL;

			m = 2*n+1;
			planarBiconnectedGraph(G1, n, m);
			addRandomLoops(G1,10);
			addRandomMultiEdges(G1, 50);
			randomizeAdjLists(G1);
			G2 = G1;

			System::usedRealTime(T);
			if(firstBL) {
				bl.planarEmbed(G1);
				tBL = System::usedRealTime(T);

				bm.planarEmbed(G2,dummyList);
				tBM = System::usedRealTime(T);
			} else {
				bm.planarEmbed(G1,dummyList);
				tBM = System::usedRealTime(T);

				bl.planarEmbed(G2);
				tBL = System::usedRealTime(T);
			}

			if(!G1.representsCombEmbedding()) {
				if(firstBL) failsBL++; else failsBM++;
			}
			if(!G2.representsCombEmbedding()) {
				if(!firstBL) failsBL++; else failsBM++;
			}

			msecBL += tBL;
			msecBM += tBM;
			firstBL = !firstBL;

			m = 5*n/2+1;
			planarBiconnectedGraph(G1, n, m);
			addRandomLoops(G1,10);
			addRandomMultiEdges(G1, 50);
			randomizeAdjLists(G1);
			G2 = G1;

			System::usedRealTime(T);
			if(firstBL) {
				bl.planarEmbed(G1);
				tBL = System::usedRealTime(T);

				bm.planarEmbed(G2,dummyList);
				tBM = System::usedRealTime(T);
			} else {
				bm.planarEmbed(G1,dummyList);
				tBM = System::usedRealTime(T);

				bl.planarEmbed(G2);
				tBL = System::usedRealTime(T);
			}

			if(!G1.representsCombEmbedding()) {
				if(firstBL) failsBL++; else failsBM++;
			}
			if(!G2.representsCombEmbedding()) {
				if(!firstBL) failsBL++; else failsBM++;
			}

			msecBL += tBL;
			msecBM += tBM;
			firstBL = !firstBL;
		}
	}

	if(failsBL + failsBM > 0) result = false;

	cout.precision(3);
	cout << "\r";
	cout << "  BoothLueker:  " << 0.001*double(msecBL) << " seconds (" << failsBL << " errors)" << endl;
	cout << "  BoyerMyrvold: " << 0.001*double(msecBM) << " seconds (" << failsBM << " errors)" << endl;

	cout << "-> Planar connected graphs + edges, planarity test... " << endl;
	srand(4711);

	msecBL = msecBM = 0;
	int fails = 0;
	int nPlanar = 0, nNonPlanar = 0;

	firstBL = true;
	for(int n = min_n; n <= max_n; n += step_n)
	{
		cout << "\r" << n << flush;
		for(int i = 0; i < numGraphs; ++i)
		{
			bool resPM, resBM;
			__int64 tBL, tBM;

			int nb = n/10;
			int addEdges = randomNumber(0,4);
			planarCNBGraph(G1, nb, 3*nb/2, 10);
			addRandomEdges(G1,addEdges);
			addRandomMultiEdges(G1, 50);
			randomizeAdjLists(G1);
			G2 = G1;

			System::usedRealTime(T);
			if(firstBL) {
				resPM = bl.isPlanar(G1);
				tBL = System::usedRealTime(T);

				resBM = bm.isPlanar(G2);
				tBM = System::usedRealTime(T);
			} else {
				resBM = bm.isPlanar(G1);
				tBM = System::usedRealTime(T);

				resPM = bl.isPlanar(G2);
				tBL = System::usedRealTime(T);
			}
			firstBL = !firstBL;

			msecBL += tBL;
			msecBM += tBM;
			if(resPM != resBM)
				fails++;
			else if(resPM)
				nPlanar++;
			else
				nNonPlanar++;

			nb = n/10;
			addEdges = randomNumber(0,3);
			planarCNBGraph(G1, nb, 3*nb/2, 10);
			addRandomEdges(G1,addEdges);
			addRandomMultiEdges(G1, 60);
			randomizeAdjLists(G1);
			G2 = G1;

			System::usedRealTime(T);
			if(firstBL) {
				resPM = bl.isPlanar(G1);
				tBL = System::usedRealTime(T);

				resBM = bm.isPlanar(G2);
				tBM = System::usedRealTime(T);
			} else {
				resBM = bm.isPlanar(G1);
				tBM = System::usedRealTime(T);

				resPM = bl.isPlanar(G2);
				tBL = System::usedRealTime(T);
			}
			firstBL = !firstBL;

			msecBL += tBL;
			msecBM += tBM;
			if(resPM != resBM)
				fails++;
			else if(resPM)
				nPlanar++;
			else
				nNonPlanar++;
		}
	}

	if(fails > 0) result = false;

	cout.precision(3);
	cout << "\r";
	cout << "  planar:       " << nPlanar << endl;
	cout << "  non-planar:   " << nNonPlanar << endl;
	cout << "  fails:        " << fails << endl;
	cout << "  BoothLueker:  " << 0.001*double(msecBL) << " seconds" << endl;
	cout << "  BoyerMyrvold: " << 0.001*double(msecBM) << " seconds" << endl;

	cout << "-> Planar biconnected graphs + 5 edges, planar subgraph... " << endl;
	srand(4711);

	FastPlanarSubgraph fps;
	fps.runs(16);

	msecBL = 0;
	fails = 0;
	int sumDelEdges = 0, ng = 0;

	for(int n = min_n; n <= max_n; n += step_n)
	{
		cout << "\r" << n << flush;
		for(int i = 0; i < numGraphs; ++i)
		{
			List<edge> delEdges;

			int m = 3*n/2;
			planarBiconnectedGraph(G1, n , m); ++ng;
			addRandomMultiEdges(G1,5);
			addRandomEdges(G1,5);
			randomizeAdjLists(G1);

			System::usedRealTime(T);
			fps.call(G1,delEdges);
			msecBL += System::usedRealTime(T);

			sumDelEdges += delEdges.size();
			removeEdges(G1,delEdges);
			if(bm.isPlanar(G1) == false) fails++;

			m = 2*n;
			planarBiconnectedGraph(G1, n , m); ++ng;
			addRandomMultiEdges(G1,5);
			addRandomEdges(G1,5);
			randomizeAdjLists(G1);

			System::usedRealTime(T);
			fps.call(G1,delEdges);
			msecBL += System::usedRealTime(T);

			sumDelEdges += delEdges.size();
			removeEdges(G1,delEdges);
			if(bm.isPlanar(G1) == false) fails++;

			m = 5*n/2;
			planarBiconnectedGraph(G1, n , m); ++ng;
			addRandomMultiEdges(G1,5);
			addRandomEdges(G1,5);
			randomizeAdjLists(G1);

			System::usedRealTime(T);
			fps.call(G1,delEdges);
			msecBL += System::usedRealTime(T);

			sumDelEdges += delEdges.size();
			removeEdges(G1,delEdges);
			if(bm.isPlanar(G1) == false) fails++;
		}
	}

	if(fails > 0) result = false;

	cout.precision(3);
	cout << "\r";
	cout << "  time:    " << 0.001*double(msecBL) << " seconds (avg. " << double(msecBL)/ng << " msec)" << endl;
	cout << "  removed: " << sumDelEdges << " (avg. " << double(sumDelEdges)/ng << ")" << endl;
	cout << "  fails:   " << fails << endl;

	cout << "-> Runtime test: Planar biconnected graphs, planar embedding... " << endl;
	cout << "   nodes  BL time [errors]  BM time [errors]" << endl;
	srand(4711);

	firstBL = true;

	const int nr_min = 1<<8;
	const int nr_max = 1<<14;

	for(int n = nr_min; n <= nr_max; n *= 2)
	{
		cout << n << ":\t" << flush;
		SList< KuratowskiWrapper > dummyList;
		msecBL  = msecBM  = 0;
		failsBL = failsBM = 0;
		for(int i = 0; i < numGraphs; ++i)
		{
			bool resPM, resBM;
			__int64 tBL, tBM;

			int m = 2*n;
			planarBiconnectedGraph(G1, n, m);
			randomizeAdjLists(G1);
			G2 = G1;

			System::usedRealTime(T);
			if(firstBL) {
				resPM = bl.planarEmbed(G1);
				tBL = System::usedRealTime(T);

				resBM = bm.planarEmbed(G2,dummyList);
				tBM = System::usedRealTime(T);
			} else {
				resBM = bm.planarEmbed(G1,dummyList);
				tBM = System::usedRealTime(T);

				resPM = bl.planarEmbed(G2);
				tBL = System::usedRealTime(T);
			}
			if(!G1.representsCombEmbedding()) {
				if(firstBL) failsBL++; else failsBM++;
			}
			if(!G2.representsCombEmbedding()) {
				if(!firstBL) failsBL++; else failsBM++;
			}

			firstBL = !firstBL;

			msecBL += tBL;
			msecBM += tBM;
			if(resPM == false) failsBL++;
			if(resBM == false) failsBM++;
		}
		if(failsBL + failsBM > 0) result = false;
		cout << 0.001*double(msecBL) << " sec [" << failsBL << "] \t";
		cout << 0.001*double(msecBM) << " sec [" << failsBM << "]" << endl;
	}

	return result;
}
