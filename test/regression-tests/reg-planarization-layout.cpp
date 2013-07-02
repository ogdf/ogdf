//*********************************************************
//  Regression test for planarization layout algorithms
//
//  Tested classes:
//    - PlanarizationLayout
//    - PlanarizationGridLayout
//
//  Author: Carsten Gutwenger
//*********************************************************

#include <ogdf/basic/graph_generators.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/planarity/PlanarizationLayout.h>
#include <ogdf/planarity/PlanarizationGridLayout.h>
#include <ogdf/orthogonal/OrthoLayout.h>
#include <ogdf/planarity/SubgraphPlanarizer.h>
#include <ogdf/planarity/FixedEmbeddingInserter.h>
#include <ogdf/planarity/FastPlanarSubgraph.h>
#include <ogdf/planarity/FixedEmbeddingInserter.h>
#include <ogdf/planarity/VariableEmbeddingInserter.h>
#include <ogdf/planarlayout/MixedModelLayout.h>
#include <ogdf/planarlayout/MMCBLocalStretch.h>


using namespace ogdf;


void createAlmostPlanarGraph(Graph &G, int n, int m, int add_m)
{
	planarBiconnectedGraph(G,n,m);

	Array<node> table(n);

	int i = 0;
	node v;
	forall_nodes(v,G)
		table[i++] = v;

	for(i = 0; i < add_m; ++i)
		G.newEdge(table[randomNumber(0,n-1)], table[randomNumber(0,n-1)]);

	makeSimpleUndirected(G);
}


bool regPlanarizationLayout()
{
	const int numGraphs = 5;
	const int min_n = 25;
	const int max_n = 100;
	const int step_n = 25;

	int totalGraphs = numGraphs * (1+(max_n-min_n)/step_n) * 3 * 2;

	PlanarizationLayout pl;
	PlanarizationGridLayout pgl;

	OrthoLayout *ortho = new OrthoLayout;
	pl.setPlanarLayouter(ortho);

	Graph G;
	GridLayout gl(G);
	GraphAttributes GA(G);

	cout << "-> Var inserter, biconnected graphs... " << endl;
	srand(4711);

	VariableEmbeddingInserter *pVarInserter = new VariableEmbeddingInserter;
	pVarInserter->removeReinsert(rrNone);

	SubgraphPlanarizer *pCrossMin = new SubgraphPlanarizer;
	pCrossMin->setSubgraph(new FastPlanarSubgraph);
	pCrossMin->setInserter(pVarInserter);
	pCrossMin->permutations(4);

	pl .setCrossMin(pCrossMin->clone());
	pgl.setCrossMin(pCrossMin);

	__int64 T;
	__int64 msec = 0;
	System::usedRealTime(T);

	for(int n = min_n; n <= max_n; n += step_n)
	{
		cout << "\r" << n << flush;
		for(int i = 0; i < numGraphs; ++i)
		{
			int m = 3*n/2;
			randomBiconnectedGraph(G, n, m);

			System::usedRealTime(T);
			pl .call(GA);
			msec += System::usedRealTime(T);

			makeSimpleUndirected(G);
			System::usedRealTime(T);
			pgl.call(GA);
			msec += System::usedRealTime(T);

			m = 2*n;
			randomBiconnectedGraph(G, n, m);

			System::usedRealTime(T);
			pl .call(GA);
			msec += System::usedRealTime(T);

			makeSimpleUndirected(G);
			System::usedRealTime(T);
			pgl.call(GA);
			msec += System::usedRealTime(T);

			m = 5*n/2;
			randomBiconnectedGraph(G, n, m);

			System::usedRealTime(T);
			pl .call(GA);
			msec += System::usedRealTime(T);

			makeSimpleUndirected(G);
			System::usedRealTime(T);
			pgl.call(GA);
			msec += System::usedRealTime(T);
		}
	}

	cout.precision(3);
	cout << "\r    " << std::fixed << 0.001*double(msec) << " seconds (avg. " <<
		double(msec) / totalGraphs <<
		" ms per graph)" << endl;


	cout << "-> Fix inserter, beautifier (grid), simple graphs... " << endl;
	srand(4711);

	FixedEmbeddingInserter *pFixInserter = new FixedEmbeddingInserter;
	pFixInserter->removeReinsert(rrAll);

	pCrossMin = new SubgraphPlanarizer;
	pCrossMin->setSubgraph(new FastPlanarSubgraph);
	pCrossMin->setInserter(pFixInserter);
	pCrossMin->permutations(1);

	pl .setCrossMin(pCrossMin->clone());
	pgl.setCrossMin(pCrossMin);

	MixedModelLayout *mml = new MixedModelLayout;
	mml->setCrossingsBeautifier(new MMCBLocalStretch);
	pgl.setPlanarLayouter(mml);

	msec = 0;
	for(int n = min_n; n <= max_n; n += step_n)
	{
		cout << "\r" << n << flush;
		for(int i = 0; i < numGraphs; ++i)
		{
			int m = 3*n/2;
			randomSimpleGraph(G, n, m);

			System::usedRealTime(T);
			pl .call(GA);
			pgl.callGrid(G,gl);
			msec += System::usedRealTime(T);

			m = 2*n;
			randomSimpleGraph(G, n, m);

			System::usedRealTime(T);
			pl .call(GA);
			pgl.callGrid(G,gl);
			msec += System::usedRealTime(T);

			m = 5*n/2;
			randomSimpleGraph(G, n, m);

			System::usedRealTime(T);
			pl .call(GA);
			pgl.callGrid(G,gl);
			msec += System::usedRealTime(T);
		}
	}

	cout.precision(3);
	cout << "\r    " << std::fixed << 0.001*double(msec) << " seconds (avg. " <<
		double(msec) / totalGraphs <<
		" ms per graph)" << endl;


	cout << "-> Runtime test, almost planar graphs... " << endl;
	srand(4711);

	const int nr_min = 32;
	const int nr_max = 1024;

	for(int n = nr_min; n <= nr_max; n *= 2)
	{
		msec = 0;
		for(int i = 0; i < numGraphs; ++i)
		{
			createAlmostPlanarGraph(G, n, 3*n/2, 10);

			System::usedRealTime(T);
			pl .call(GA);
			pgl.callGrid(G, gl);

			msec += System::usedRealTime(T);
		}
		cout.precision(3);
		cout << n << ": \t" << std::fixed << 0.001*double(msec) << " seconds (avg. " <<
			double(msec) / (2*numGraphs) << " ms per graph)" << endl;
	}

	return true;
}
