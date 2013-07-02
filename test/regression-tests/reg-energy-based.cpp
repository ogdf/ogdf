//*********************************************************
//  Regression test for energy-based layout algorithms
//
//  Tested classes:
//    - FMMMLayout
//    - SpringEmbedderFR
//    - GEMLayout
//    - DavidsonHarelLayout
//
//  Author: Carsten Gutwenger
//*********************************************************

#include <ogdf/energybased/FMMMLayout.h>
#include <ogdf/energybased/SpringEmbedderFR.h>
#include <ogdf/energybased/GEMLayout.h>
#include <ogdf/energybased/DavidsonHarelLayout.h>
#include <ogdf/basic/graph_generators.h>
#include <ogdf/basic/simple_graph_alg.h>


using namespace ogdf;


void randomLayout(GraphAttributes &GA)
{
	const Graph &G = GA.constGraph();
	int max_x = (int)(2.0f * sqrt((float)G.numberOfNodes()));
	int max_y = max_x;

	node v;
	forall_nodes(v,G) {
		GA.x(v) = randomNumber(0,max_x);
		GA.y(v) = randomNumber(0,max_y);
	}
}

bool regEnergyBased()
{
	const int numGraphs = 10;
	const int min_n = 25;
	const int max_n = 100;
	const int step_n = 25;

	int totalGraphs = numGraphs * (1+(max_n-min_n)/step_n) * 3;

	Graph G;
	GraphAttributes GA(G);

	FMMMLayout          fmmm;
	SpringEmbedderFR    frl;
	GEMLayout           gem;
	DavidsonHarelLayout dhl;

	fmmm.newInitialPlacement(true);

	cout << "-> Defaults, biconnected graphs... " << endl;
	srand(4711);

	__int64 T;
	__int64 msecFMMM, msecFR, msecGEM, msecDH;

	msecFMMM = msecFR = msecGEM = msecDH = 0;
	for(int n = min_n; n <= max_n; n += step_n)
	{
		cout << "\r" << n << flush;
		for(int i = 0; i < numGraphs; ++i)
		{
			int m = 3*n/2;
			randomBiconnectedGraph(G, n, m);
			makeSimpleUndirected(G);

			randomLayout(GA);
			System::usedRealTime(T);
			fmmm.call(GA);
			msecFMMM += System::usedRealTime(T);

			randomLayout(GA);
			System::usedRealTime(T);
			frl .call(GA);
			msecFR += System::usedRealTime(T);

			randomLayout(GA);
			System::usedRealTime(T);
			gem .call(GA);
			msecGEM += System::usedRealTime(T);

			randomLayout(GA);
			System::usedRealTime(T);
			dhl.call(GA);
			msecDH += System::usedRealTime(T);

			m = 2*n;
			randomBiconnectedGraph(G, n, m);
			makeSimpleUndirected(G);

			randomLayout(GA);
			System::usedRealTime(T);
			fmmm.call(GA);
			msecFMMM += System::usedRealTime(T);

			randomLayout(GA);
			System::usedRealTime(T);
			frl .call(GA);
			msecFR += System::usedRealTime(T);

			randomLayout(GA);
			System::usedRealTime(T);
			gem .call(GA);
			msecGEM += System::usedRealTime(T);

			randomLayout(GA);
			System::usedRealTime(T);
			dhl.call(GA);
			msecDH += System::usedRealTime(T);

			m = 5*n/2;
			randomBiconnectedGraph(G, n, m);
			makeSimpleUndirected(G);

			randomLayout(GA);
			System::usedRealTime(T);
			fmmm.call(GA);
			msecFMMM += System::usedRealTime(T);

			randomLayout(GA);
			System::usedRealTime(T);
			frl .call(GA);
			msecFR += System::usedRealTime(T);

			randomLayout(GA);
			System::usedRealTime(T);
			gem .call(GA);
			msecGEM += System::usedRealTime(T);

			randomLayout(GA);
			System::usedRealTime(T);
			dhl.call(GA);
			msecDH += System::usedRealTime(T);
		}
	}

	cout.precision(3);
	cout << "\r  FMMM: " << std::fixed << 0.001*double(msecFMMM) << " seconds (avg. " <<
		double(msecFMMM) / totalGraphs <<
		" ms per graph)" << endl;
	cout << "\r  FR:   " << std::fixed << 0.001*double(msecFR) << " seconds (avg. " <<
		double(msecFR) / totalGraphs <<
		" ms per graph)" << endl;
	cout << "\r  GEM:  " << std::fixed << 0.001*double(msecGEM) << " seconds (avg. " <<
		double(msecGEM) / totalGraphs <<
		" ms per graph)" << endl;
	cout << "\r  DH:   " << std::fixed << 0.001*double(msecDH)  << " seconds (avg. " <<
		double(msecDH) / totalGraphs <<
		" ms per graph)" << endl;


	cout << "\n-> High-quality, planar-connected graphs... " << endl;
	srand(4711);

	fmmm.qualityVersusSpeed(FMMMLayout::qvsGorgeousAndEfficient);
	fmmm.useHighLevelOptions(true);
	frl.iterations(1000);

	msecFMMM = msecFR = msecGEM = msecDH = 0;
	for(int n = min_n; n <= max_n; n += step_n)
	{
		int b = n/25;

		cout << "\r" << n << flush;
		for(int i = 0; i < numGraphs; ++i)
		{
			planarCNBGraph(G, 25, 38, b);

			randomLayout(GA);
			System::usedRealTime(T);
			fmmm.call(GA);
			msecFMMM += System::usedRealTime(T);

			randomLayout(GA);
			System::usedRealTime(T);
			frl .call(GA);
			msecFR += System::usedRealTime(T);

			randomLayout(GA);
			System::usedRealTime(T);
			gem .call(GA);
			msecGEM += System::usedRealTime(T);

			randomLayout(GA);
			System::usedRealTime(T);
			dhl.call(GA);
			msecDH += System::usedRealTime(T);

			planarCNBGraph(G, 25, 50, b);

			randomLayout(GA);
			System::usedRealTime(T);
			fmmm.call(GA);
			msecFMMM += System::usedRealTime(T);

			randomLayout(GA);
			System::usedRealTime(T);
			frl .call(GA);
			msecFR += System::usedRealTime(T);

			randomLayout(GA);
			System::usedRealTime(T);
			gem .call(GA);
			msecGEM += System::usedRealTime(T);

			randomLayout(GA);
			System::usedRealTime(T);
			dhl.call(GA);
			msecDH += System::usedRealTime(T);

			planarCNBGraph(G, 25, 63, b);

			randomLayout(GA);
			System::usedRealTime(T);
			fmmm.call(GA);
			msecFMMM += System::usedRealTime(T);

			randomLayout(GA);
			System::usedRealTime(T);
			frl .call(GA);
			msecFR += System::usedRealTime(T);

			randomLayout(GA);
			System::usedRealTime(T);
			gem .call(GA);
			msecGEM += System::usedRealTime(T);

			randomLayout(GA);
			System::usedRealTime(T);
			dhl.call(GA);
			msecDH += System::usedRealTime(T);
		}
	}

	cout.precision(3);
	cout << "\r  FMMM: " << std::fixed << 0.001*double(msecFMMM) << " seconds (avg. " <<
		double(msecFMMM) / totalGraphs <<
		" ms per graph)" << endl;
	cout << "\r  FR:   " << std::fixed << 0.001*double(msecFR) << " seconds (avg. " <<
		double(msecFR) / totalGraphs <<
		" ms per graph)" << endl;
	cout << "\r  GEM:  " << std::fixed << 0.001*double(msecGEM) << " seconds (avg. " <<
		double(msecGEM) / totalGraphs <<
		" ms per graph)" << endl;
	cout << "\r  DH:   " << std::fixed << 0.001*double(msecDH) << " seconds (avg. " <<
		double(msecDH) / totalGraphs <<
		" ms per graph)" << endl;


	cout << "\n-> Runtime test layout (FMMM - FR - GEM)... " << endl;
	srand(4711);

	fmmm.qualityVersusSpeed(FMMMLayout::qvsNiceAndIncredibleSpeed);

	const int nr_min = 64;
	const int nr_max = 4096;

	for(int n = nr_min; n <= nr_max; n *= 2)
	{
		msecFMMM = msecFR = msecGEM = msecDH = 0;
		for(int i = 0; i < numGraphs; ++i)
		{
			planarBiconnectedGraph(G, n, 3*n/2);

			randomLayout(GA);
			System::usedRealTime(T);
			fmmm.call(GA);
			msecFMMM += System::usedRealTime(T);

			randomLayout(GA);
			System::usedRealTime(T);
			frl .call(GA);
			msecFR += System::usedRealTime(T);

			randomLayout(GA);
			System::usedRealTime(T);
			gem .call(GA);
			msecGEM += System::usedRealTime(T);
		}
		cout.precision(3);
		cout << n << ": \t" << std::fixed <<
			0.001*double(msecFMMM)/numGraphs << " s  - \t" <<
			0.001*double(msecFR)  /numGraphs << " s  - \t" <<
			0.001*double(msecGEM) /numGraphs << " s" << endl;
	}

	return true;
}
