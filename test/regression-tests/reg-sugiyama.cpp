//*********************************************************
//  Regression test for Sugiyama layout
//
//  Tested classes:
//    - SugiyamaLayout
//
//  Author: Carsten Gutwenger
//*********************************************************

#include <ogdf/layered/SugiyamaLayout.h>
#include <ogdf/layered/FastHierarchyLayout.h>
#include <ogdf/layered/MedianHeuristic.h>
#include <ogdf/layered/BarycenterHeuristic.h>
#include <ogdf/layered/OptimalHierarchyLayout.h>
#include <ogdf/layered/OptimalRanking.h>
#include <ogdf/layered/LongestPathRanking.h>
#include <ogdf/layered/GreedyCycleRemoval.h>
#include <ogdf/layered/DfsAcyclicSubgraph.h>

#include <ogdf/basic/graph_generators.h>
#include <ogdf/basic/simple_graph_alg.h>


using namespace ogdf;

void insertGraph(Graph &g, const Graph &g2)
{
	NodeArray<node> map(g2);

	node v;
	forall_nodes(v,g2)
		map[v] = g.newNode();

	edge e;
	forall_edges(e,g2)
		g.newEdge(map[e->source()], map[e->target()]);
}

void createDisconnectedGraph(
	Graph &G,
	int n_max,
	double density_min,
	double density_max,
	int cc)
{
	int n = randomNumber(3,n_max);
	int m = randomNumber((int)(density_min*n),(int)(density_max*n));
	int b = n/25+1;
	planarCNBGraph(G, n/b+1, m/b+1, b);

	for(int i = cc-1; i > 0; --i) {
		Graph g2;
		n = randomNumber(3,n);
		m = randomNumber((int)(density_min*n),(int)(density_max*n));
		b = n/25+1;
		planarCNBGraph(g2, n/b+1, m/b+1, b);
		insertGraph(G,g2);
	}
}


bool regSugiyama()
{
	const int numGraphs = 10;
	const int min_n = 50;
	const int max_n = 250;
	const int step_n = 50;

	int totalGraphs = numGraphs * (1+(max_n-min_n)/step_n) * 3;

	Graph G;
	GraphAttributes GA(G);

	SugiyamaLayout sugi;
	sugi.setLayout(new FastHierarchyLayout);

	cout << "-> Defaults, sT-graphs... " << endl;
	srand(4711);

	__int64 T;
	__int64 msec = 0;
	int cr = 0;

	for(int n = min_n; n <= max_n; n += step_n)
	{
		cout << "\r" << n << flush;
		for(int i = 0; i < numGraphs; ++i)
		{
			int m = 3*n/2;
			randomHierarchy(G, n, m, false, true, false);

			System::usedRealTime(T);
			sugi.call(GA);
			msec += System::usedRealTime(T);
			cr += sugi.numberOfCrossings();

			m = 2*n;
			randomHierarchy(G, n, m, false, true, false);

			System::usedRealTime(T);
			sugi.call(GA);
			msec += System::usedRealTime(T);
			cr += sugi.numberOfCrossings();

			m = 5*n/2;
			randomHierarchy(G, n, m, false, true, false);

			System::usedRealTime(T);
			sugi.call(GA);
			msec += System::usedRealTime(T);
			cr += sugi.numberOfCrossings();
		}
	}

	cout.precision(3);
	cout << "\r    " << std::fixed << 0.001*double(msec) << " s, " <<
		cr << " crossings (avg. " <<
		double(msec) / totalGraphs <<	" ms, " <<
		double(cr)/totalGraphs << " crossings per graph)" << endl;


#ifdef OGDF_LP_SOLVER
	sugi.setLayout(new OptimalHierarchyLayout);
#endif

	OptimalRanking *optr = new OptimalRanking;
	optr->setSubgraph(new GreedyCycleRemoval);

	sugi.setRanking(optr);
	sugi.setCrossMin(new MedianHeuristic);
	sugi.transpose(false);

#ifdef OGDF_LP_SOLVER
	cout << "-> OptimalRanking, Median, OptimalHierarchy, planarCNB-graphs... " << endl;
#else
	cout << "-> OptimalRanking, Median, planarCNB-graphs... " << endl;
#endif
	srand(4711);

	msec = 0; cr = 0;
	for(int n = min_n; n <= max_n; n += step_n)
	{
		int b = n/25;

		cout << "\r" << n << flush;
		for(int i = 0; i < numGraphs; ++i)
		{
			planarCNBGraph(G, 25, 38, b);

			System::usedRealTime(T);
			sugi.call(GA);
			msec += System::usedRealTime(T);
			cr += sugi.numberOfCrossings();

			planarCNBGraph(G, 25, 50, b);

			System::usedRealTime(T);
			sugi.call(GA);
			msec += System::usedRealTime(T);
			cr += sugi.numberOfCrossings();

			planarCNBGraph(G, 25, 63, b);

			System::usedRealTime(T);
			sugi.call(GA);
			msec += System::usedRealTime(T);
			cr += sugi.numberOfCrossings();
		}
	}

	cout.precision(3);
	cout << "\r    " << std::fixed << 0.001*double(msec) << " s, " <<
		cr << " crossings (avg. " <<
		double(msec) / totalGraphs <<	" ms, " <<
		double(cr)/totalGraphs << " crossings per graph)" << endl;


#ifdef OGDF_LP_SOLVER
	cout << "-> OptimalRanking, Median, OptimalHierarchy, disconnected graphs... " << endl;
#else
	cout << "-> OptimalRanking, Median, disconnected graphs... " << endl;
#endif
	srand(4711);

	sugi.runs(40);
	sugi.transpose(true);

	msec = 0; cr = 0;
	for(int n = min_n; n <= max_n; n += step_n)
	{
		cout << "\r" << n << flush;
		for(int i = 0; i < numGraphs; ++i)
		{
			createDisconnectedGraph(G, n/7, 1.4, 2.6, 7);

			System::usedRealTime(T);
			sugi.call(GA);
			msec += System::usedRealTime(T);
			cr += sugi.numberOfCrossings();
		}
	}

	cout.precision(3);
	cout << "\r    " << std::fixed << 0.001*double(msec) << " s, " <<
		cr << " crossings (avg. " <<
		double(msec) / totalGraphs <<	" ms, " <<
		double(cr)/totalGraphs << " crossings per graph)" << endl;


	cout << "\n-> Runtime test planar biconnected graphs... " << endl;
	srand(4711);

	sugi.runs(10);
	sugi.transpose(false);
	sugi.setCrossMin(new BarycenterHeuristic);
	sugi.setLayout(new FastHierarchyLayout);

	LongestPathRanking *lpr = new LongestPathRanking;
	lpr->setSubgraph(new DfsAcyclicSubgraph);
	sugi.setRanking(lpr);

	const int nr_min = 64;
	const int nr_max = 4096;

	for(int n = nr_min; n <= nr_max; n *= 2)
	{
		msec = 0; cr = 0;
		for(int i = 0; i < numGraphs; ++i)
		{
			planarBiconnectedGraph(G, n, 3*n/2);

			System::usedRealTime(T);
			sugi.call(GA);
			msec += System::usedRealTime(T);
			cr += sugi.numberOfCrossings();
		}
		cout.precision(3);
		cout << n << ": \t" << std::fixed <<
			0.001*double(msec)/numGraphs << " s  \t" <<
			double(cr)/numGraphs << " crossings" << endl;
	}

	return true;
}
