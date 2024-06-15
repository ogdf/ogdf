#include <ogdf/cluster/sync_plan/PQPlanarity.h>
#include <ogdf/cluster/sync_plan/PipeOrder.h>

bool PipeQueueByDegreePreferContract::comparePipes(const Pipe* x, const Pipe* y) const {
	if (x->heap_data != y->heap_data) {
		return x->heap_data < y->heap_data;
	}

	if (invert_degree) {
		return x->degree() < y->degree();
	} else {
		return x->degree() > y->degree();
	}
}

bool PipeQueueByDegreePreferContract::isQueue1(Pipe* p) const {
	if (p->pipe_priority >= 0) {
		return true;
	}
	bool ret = p->degree() <= 3 || PQ->canContract(p);
	if (invert_contract) {
		ret = !ret;
	}
	return ret;
}
