#include "PQPlanarity.h"
#include "utils/Logging.h"

class UndoMakeWheel : public PQPlanarity::UndoOperation {
public:
	int centre_idx;

	UndoMakeWheel(node centre) : centre_idx(centre->index()) { }

	void undo(PQPlanarity& pq) override { pq.contractWheel(pq.nodeFromIndex(centre_idx)); }

	ostream& print(ostream& os) const override {
		return os << "UndoMakeWheel(" << centre_idx << ")";
	}
};

void PQPlanarity::makeWheel(node centre, bool update_cuts) {
	PQ_PROFILE_START("makeWheel")
	bool is_cut = update_cuts && components.isCutVertex(centre);
	log.lout(Logger::Level::High) << "MAKE WHEEL centre " << fmtPQNode(centre) << " (update cut "
								  << update_cuts << is_cut << ")" << std::endl;
	Logger::Indent _(&log);
	log.lout(Logger::Level::Medium) << printIncidentEdges(centre->adjEntries) << std::endl;
	NodeArray<SList<adjEntry>> block_neigh;
	if (is_cut) {
		block_neigh.init(components.bcTree());
	}

	node prev_ray = nullptr, first_ray = nullptr;
	node centre_bc = components.biconnectedComponent(centre);
	List<adjEntry> adjEntries;
	centre->allAdjEntries(adjEntries);
	for (adjEntry adj : adjEntries) {
		node neigh = adj->twinNode();
		edge the_edge = adj->theEdge();
		bool reverse = adj->isSource();
		if (reverse) {
			G->reverseEdge(the_edge);
		}
		edge new_edge = G->split(the_edge);
		OGDF_ASSERT(neigh == the_edge->source());
		OGDF_ASSERT(the_edge->target() == new_edge->source());
		OGDF_ASSERT(new_edge->target() == centre);
		node spike = new_edge->source();
		is_wheel[spike] = true;
		if (reverse) {
			G->reverseEdge(the_edge);
			G->reverseEdge(new_edge);
		}

		if (GA != nullptr) {
			std::ostringstream ss;
			ss << "WheelSpike " << spike->index() << " for " << centre->index() << "-"
			   << neigh->index();
			GA->label(spike) = ss.str();
		}
		if (is_cut) {
			node neigh_bc = components.biconnectedComponent(neigh);
			if (components.isCutComponent(neigh_bc)) {
				block_neigh[components.findCommonBiconComp(centre_bc, neigh_bc)].pushBack(adj);
			} else {
				block_neigh[neigh_bc].pushBack(adj);
			}
		}
		components.nodeInserted(spike, centre_bc);

		if (prev_ray != nullptr) {
			G->newEdge(prev_ray, spike);
		}
		prev_ray = spike;
		if (first_ray == nullptr) {
			first_ray = spike;
		}
	}
	if (first_ray != nullptr) {
		G->newEdge(prev_ray, first_ray);
	}
	is_wheel[centre] = true;
	log.lout(Logger::Level::Medium) << printIncidentEdges(centre->adjEntries) << std::endl;

	if (is_cut) {
		components.cutReplacedByWheel(centre, block_neigh); // TODO use union-find to merge blocks faster
	}

	pushUndoOperation(new UndoMakeWheel(centre));
	PQ_PROFILE_STOP("makeWheel")
}

void PQPlanarity::contractWheel(node centre) {
	PQ_PROFILE_START("contractWheel")
	OGDF_ASSERT(is_wheel[centre]);
	log.lout(Logger::Level::High)
			<< "CONTRACT WHEEL centre " << fmtPQNode(centre, false) << std::endl;
	Logger::Indent _(&log);
	log.lout(Logger::Level::Medium) << printIncidentEdges(centre->adjEntries) << std::endl;
	List<adjEntry> adjEntries;
	centre->allAdjEntries(adjEntries);
	for (adjEntry adj : adjEntries) {
		OGDF_ASSERT(is_wheel[adj->twinNode()]);
		if (!adj->isSource()) {
			G->reverseEdge(adj->theEdge());
		}
		node res = G->contract(adj->theEdge());
		OGDF_ASSERT(res == centre);
	}
	is_wheel[centre] = false;
	log.lout(Logger::Level::Medium) << printIncidentEdges(centre->adjEntries) << std::endl;
	PQ_PROFILE_STOP("contractWheel")
}
