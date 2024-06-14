#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/simple_graph_alg.h>

#include <json.hpp>
#include <sstream>

#include "return.h"

using nlohmann::json;

using namespace std;
using namespace ogdf;
static std::string progname;

int main(int argc, char* argv[]) {
	progname = argv[0];
	for (int i = 1; i < argc; ++i) {
		Graph G;
		GraphAttributes GA(G, GraphAttributes::nodeLabel | GraphAttributes::edgeLabel);
		string file = argv[i];
		int read_ret = readCG(file, Logger::Level::Default, GA, G);
		if (read_ret != SUCCESS) {
			Logger::slout(Logger::Level::Alarm)
					<< "IO Error: Couldn't read graph " << file << "!" << std::endl;
			return read_ret;
		}
		PQPlanarity pq(&G);
		try {
			std::ifstream i(file + ".json");
			if (!i.good()) {
				Logger::slout(Logger::Level::Alarm)
						<< "IO Error: Couldn't read config json " << file << ".json!" << std::endl;
				return ERROR_IO;
			}
			json j;
			i >> j;
			PQPlanOptions::applyConfigJSON(G, GA, pq, j);
		} catch (json::parse_error& e) {
			Logger::slout(Logger::Level::Alarm) << "IO Error: Couldn't parse config json " << file
												<< ".json: " << e.what() << std::endl
												<< "exception id: " << e.id << std::endl
												<< "byte position of error: " << e.byte << std::endl;
			return ERROR_IO;
		}

		int max_deg = 0, tot_deg = 0;
		for (auto& p : pq.matchings.getPipes()) {
			tot_deg += p.degree();
			if (p.degree() > max_deg) {
				max_deg = p.degree();
			}
		}

		json j {
				{"file", file},
				{"nodes", G.numberOfNodes()},
				{"edges", G.numberOfEdges()},
				{"planar", isPlanar(G)},
				{"pipes", pq.matchings.getPipeCount()},
				{"pipe-max-deg", max_deg},
				{"pipe-total-deg", tot_deg},
				{"partitions", pq.partitions.partitionCount()},
				{"qvertices", pq.partitions.qVertexCount()},
				{"components", pq.getComponents().connectedCount()},
				{"bc-nodes", pq.getComponents().bcTree().numberOfNodes()},
				{"bc-egdes", pq.getComponents().bcTree().numberOfEdges()},
		};
		cout << j << endl;
	}
	return SUCCESS;
}
