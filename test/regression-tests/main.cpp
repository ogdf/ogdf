#include <ogdf/basic/basic.h>

using namespace ogdf;

// forward declarations of regression tests
extern bool regPlanarLayout();
extern bool regPlanarizationLayout();
extern bool regEnergyBased();
extern bool regSugiyama();
extern bool regPlanarityTest();
extern bool regSteinerTree();
extern bool regLCA();

struct regTest {
	const char *what;
	const char *list;
	bool(*call)();
};

static struct regTest regTests[] = {
	{
		"planar layouts",
		"PlanarStraightLayout, PlanarDrawLayout, MixedModelLayout",
		regPlanarLayout
	},
	{
		"planarization layouts",
		"PlanarizationLayout, PlanarizationGridLayout",
		regPlanarizationLayout
	},
	{
		"energy-based layouts",
		"FMMMLayout, SpringEmbedderFR, GEMLayout, DavidsonHarelLayout",
		regEnergyBased
	},
	{
		"Sugiyama layout",
		"SugiyamaLayout, LongestPathRanking, OptimalRanking,\n"
		"       BarycenterHeuristic, MedianHeuristic,\n"
		"       FastHierarchyLayout, OptimalHierarchyLayout",
		regSugiyama,
	},
	{
		"planarity testing and embedding",
		"BoothLueker, BoyerMyrvold",
		regPlanarityTest
	},
	{
		"Lowest Common Ancestor testing",
		"LCA",
		regLCA
	},
	{ NULL, NULL, NULL }
};

int regressionMain(int argc, char *const argv[])
{
	cout << "**************************************************\n";
	cout << "*              OGDF Regression Test              *\n";
	cout << "**************************************************\n\n";
	cout << "System:         " << Configuration::whichSystem() << "\n";
	cout << "LP-Solver:      " << Configuration::whichLPSolver() << "\n";
	cout << "Abacus:         " << std::boolalpha << Configuration::haveAbacus() << "\n";
	cout << "Memory-Manager: " << Configuration::whichMemoryManager() << "\n";
	cout << endl;

	int successful = 0;
	int failed = 0;

	for (int i = 0; regTests[i].what; ++i) {
		// if one argument is given: this is the only regression test
		// that should be run (its index in the regTests array set above)
		if (argc == 2 && atoi(argv[1]) != i) {
			continue;
		}
		cout << "Regression test for " << regTests[i].what << ".\n";
		cout << "Tests: " << regTests[i].list << "\n\n";
		if (regTests[i].call()) {
			++successful;
			cout << "Regression test successful." << endl;
		} else {
			++failed;
			cout << "Regression test failed!" << endl;
		}
		cout << "\n--------------------------------------------------\n";
	}
	cout << "successful: " << successful << endl;
	cout << "failed:     " << failed << endl;

	return failed;
}
