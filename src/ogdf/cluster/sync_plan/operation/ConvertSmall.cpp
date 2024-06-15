#include <ogdf/cluster/sync_plan/PQPlanarity.h>
#include <ogdf/cluster/sync_plan/basic/GraphUtils.h>
#include <ogdf/cluster/sync_plan/utils/Logging.h>

class UndoConvertSmall : public PQPlanarity::UndoOperation {
public:
	int small_idx, twin_idx;
	int small_first_adj_idx, twin_last_adj_idx;
#ifdef OGDF_DEBUG
	FrozenPipeBij bij;
#endif

	UndoConvertSmall(node small, node twin) : small_idx(small->index()), twin_idx(twin->index()) {
		if (small->degree() > 0) {
			small_first_adj_idx = small->adjEntries.head()->theEdge()->index();
			twin_last_adj_idx = twin->adjEntries.tail()->theEdge()->index();
		} else {
			small_first_adj_idx = twin_last_adj_idx = -1;
		}
#ifdef OGDF_DEBUG
		getFrozenPipeBijection(small, twin, bij);
#endif
	}

	void undo(PQPlanarity& pq) override {
		PQ_PROFILE_START("undo-convertSmall")
		node small = pq.nodeFromIndex(small_idx);
		node twin = pq.nodeFromIndex(twin_idx);
		if (small_first_adj_idx >= 0) {
			moveAdjToFront(*pq.G, pq.edgeFromIndex(small_first_adj_idx)->getAdj(small));
			moveAdjToBack(*pq.G, pq.edgeFromIndex(twin_last_adj_idx)->getAdj(twin));
		}
		if (pq.partitions.isQVertex(small)) {
			OGDF_ASSERT(small->degree() > 2);
			pq.partitions.releaseQVertex(small);
			pq.partitions.releaseQVertex(twin);
		} else {
			OGDF_ASSERT(small->degree() <= 2);
		}
		pq.matchings.matchNodes(small, twin);

#ifdef OGDF_DEBUG
		pq.verifyPipeBijection(small, twin, bij);
#endif
		PQ_PROFILE_STOP("undo-convertSmall")
	}

	ostream& print(ostream& os) const override {
		return os << "UndoConvertSmall(small=" << small_idx << ", twin=" << twin_idx
#ifdef OGDF_DEBUG
				  << ", bij=" << printFrozenBijection(bij)
#endif
				  << ")";
	}
};

PQPlanarity::Result PQPlanarity::convertSmall(node small) {
	if (small->degree() > 4 || partitions.isQVertex(small)) {
		return NOT_APPLICABLE;
	} else if (matchings.isMatchedPVertex(small)) {
		PQ_PROFILE_START("convertSmall")
		log.lout(Logger::Level::High) << "CONVERT SMALL degree " << small->degree() << std::endl;
		log.lout(Logger::Level::Minor) << matchings.printBijection(small) << std::endl;
		node twin = matchings.removeMatching(small);
		if (small->degree() > 2) {
			// this forces the edge bijection to stay valid
			int partition = partitions.makeQVertex(small);
			partitions.makeQVertex(twin, partition);
		}
		pushUndoOperation(new UndoConvertSmall(small, twin));
		if (small->degree() > 2) {
			if (components.isCutVertex(small)) {
				makeWheel(small);
			}
			if (components.isCutVertex(twin)) {
				makeWheel(twin);
			}
		}
		formatNode(twin);
		PQ_PROFILE_STOP("convertSmall")
	}
	formatNode(small);
	return SUCCESS;
}
