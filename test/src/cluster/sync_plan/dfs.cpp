#include <ogdf/basic/Graph.h>

using namespace ogdf;

class DFS {
	Graph G;

	NodeArray<int> discovery;
	NodeArray<int> finish;
	NodeArray<node> predecessor;

	int time = 0;

public:
	explicit DFS(const Graph& g)
		: G(g), discovery(G, -1), finish(G, -1), predecessor(G, nullptr) { }

	void dfs_visit(node u) {
		time += 1;
		discovery[u] = time;

		for (adjEntry adj : u->adjEntries) {
			node v = adj->twinNode();
			if (adj->isSource() && discovery[v] < 0) {
				predecessor[v] = u;
				dfs_visit(v);
			}
		}

		time += 1;
		finish[u] = time;
	}

	void dfs() {
		for (node n : G.nodes) {
			if (discovery[n] < 0) {
				dfs_visit(n);
			}
		}
	}
};

int main() { }
