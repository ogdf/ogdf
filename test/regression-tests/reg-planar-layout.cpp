//*********************************************************
//  Regression test for planar layout algorithms
//
//  Tested classes:
//    - PlanarStraightLayout
//    - PlanarDrawLayout
//    - MixedModelLayout
//
//  Author: Carsten Gutwenger
//*********************************************************

#include <ogdf/basic/graph_generators.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/planarlayout/PlanarStraightLayout.h>
#include <ogdf/planarlayout/PlanarDrawLayout.h>
#include <ogdf/planarlayout/MixedModelLayout.h>
#include <ogdf/planarity/EmbedderMaxFaceLayers.h>
#include <ogdf/planarity/SimpleEmbedder.h>
#include <ogdf/planarlayout/TriconnectedShellingOrder.h>
#include <ogdf/planarlayout/BiconnectedShellingOrder.h>

using namespace ogdf;

bool regPlanarLayout()
{
	const int numGraphs = 20;
	const int min_n = 50;
	const int max_n = 500;
	const int step_n = 50;

	int totalGraphs = numGraphs * (1+(max_n-min_n)/step_n) * 3 * 3;

	Graph G;
	GridLayout gl;

	PlanarStraightLayout psl;
	PlanarDrawLayout     pdl;
	MixedModelLayout     mml;

	cout << "-> Defaults, biconnected graphs... " << endl;
	srand(4711);

	__int64 T;
	__int64 msec = 0;

	for(int n = min_n; n <= max_n; n += step_n)
	{
		cout << "\r" << n << flush;
		for(int i = 0; i < numGraphs; ++i)
		{
			int m = 3*n/2;
			planarBiconnectedGraph(G, n, m);

			System::usedRealTime(T);
			psl.callGrid(G,gl);
			pdl.callGrid(G,gl);
			mml.callGrid(G,gl);
			msec += System::usedRealTime(T);

			m = 2*n;
			planarBiconnectedGraph(G, n, m);

			System::usedRealTime(T);
			psl.callGrid(G,gl);
			pdl.callGrid(G,gl);
			mml.callGrid(G,gl);
			msec += System::usedRealTime(T);

			m = 5*n/2;
			planarBiconnectedGraph(G, n, m);

			System::usedRealTime(T);
			psl.callGrid(G,gl);
			pdl.callGrid(G,gl);
			mml.callGrid(G,gl);
			msec += System::usedRealTime(T);
		}
	}

	cout.precision(3);
	cout << "\r    " << std::fixed << 0.001*double(msec) << " seconds (avg. " <<
		double(msec) / totalGraphs <<
		" ms per graph)" << endl;


	cout << "-> EmbedderMaxFace, connected graphs... " << endl;
	srand(4711);

	mml.setEmbedder(new EmbedderMaxFaceLayers);
	psl.setEmbedder(new EmbedderMaxFaceLayers);
	pdl.setEmbedder(new EmbedderMaxFaceLayers);

	msec = 0;
	for(int n = min_n; n <= max_n; n += step_n)
	{
		cout << "\r" << n << flush;
		for(int i = 0; i < numGraphs; ++i)
		{
			int b = n/25;

			planarCNBGraph(G, 25, 38, b);

			System::usedRealTime(T);
			psl.callGrid(G,gl);
			pdl.callGrid(G,gl);
			mml.callGrid(G,gl);
			msec += System::usedRealTime(T);

			planarCNBGraph(G, 25, 50, b);

			System::usedRealTime(T);
			psl.callGrid(G,gl);
			pdl.callGrid(G,gl);
			mml.callGrid(G,gl);
			msec += System::usedRealTime(T);

			planarCNBGraph(G, 25, 63, b);

			System::usedRealTime(T);
			psl.callGrid(G,gl);
			pdl.callGrid(G,gl);
			mml.callGrid(G,gl);
			msec += System::usedRealTime(T);
		}
	}

	cout.precision(3);
	cout << "\r    " << std::fixed << 0.001*double(msec) << " seconds (avg. " <<
		double(msec) / totalGraphs <<
		" ms per graph)" << endl;


	cout << "-> EmbedderMaxFace, Tric-Order, triconnected graphs... " << endl;
	srand(4711);

	psl.setShellingOrder(new TriconnectedShellingOrder);
	pdl.setShellingOrder(new TriconnectedShellingOrder);
	mml.setShellingOrder(new TriconnectedShellingOrder);

	msec = 0;
	for(int n = min_n; n <= max_n; n += step_n)
	{
		cout << "\r" << n << flush;
		for(int i = 0; i < numGraphs; ++i)
		{
			int m = 3*n/2;
			planarTriconnectedGraph(G, n, m);

			System::usedRealTime(T);
			psl.callGrid(G,gl);
			pdl.callGrid(G,gl);
			mml.callGrid(G,gl);
			msec += System::usedRealTime(T);

			m = 2*n;
			planarTriconnectedGraph(G, n, m);

			System::usedRealTime(T);
			psl.callGrid(G,gl);
			pdl.callGrid(G,gl);
			mml.callGrid(G,gl);
			msec += System::usedRealTime(T);

			m = 5*n/2;
			planarTriconnectedGraph(G, n, m);

			System::usedRealTime(T);
			psl.callGrid(G,gl);
			pdl.callGrid(G,gl);
			mml.callGrid(G,gl);
			msec += System::usedRealTime(T);
		}
	}

	cout.precision(3);
	cout << "\r    " << std::fixed << 0.001*double(msec) << " seconds (avg. " <<
		double(msec) / totalGraphs <<
		" ms per graph)" << endl;


	cout << "-> Runtime test layout, fixed embedding (ps&pd)... " << endl;
	srand(4711);

	psl.setShellingOrder(new BiconnectedShellingOrder);
	pdl.setShellingOrder(new BiconnectedShellingOrder);
	mml.setShellingOrder(new BiconnectedShellingOrder);
	mml.setEmbedder(new SimpleEmbedder);

	GraphAttributes GA(G);

	const int nr_min = 64;
	const int nr_max = 8192;

	for(int n = nr_min; n <= nr_max; n *= 2)
	{
		msec = 0;
		for(int i = 0; i < numGraphs; ++i)
		{
			planarBiconnectedGraph(G, n, 2*n);

			planarEmbed(G);

			System::usedRealTime(T);
			psl.callFixEmbed(GA);
			pdl.callFixEmbed(GA);
			mml.call(GA);
			msec += System::usedRealTime(T);
		}
		cout.precision(3);
		cout << n << ": \t" << std::fixed << 0.001*double(msec) << " seconds (avg. " <<
			double(msec) / (3*numGraphs) << " ms per graph)" << endl;
	}

	return true;
}
