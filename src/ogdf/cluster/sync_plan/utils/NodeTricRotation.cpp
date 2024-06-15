#include <ogdf/cluster/sync_plan/PQPlanarity.h>
#include <ogdf/cluster/sync_plan/utils/Logging.h>
#include <ogdf/cluster/sync_plan/utils/NodeTricRotation.h>

using namespace pc_tree;
using namespace spqr_utils;

Logger SimpleSPQRTree::log;

#define logm log.lout(Logger::Level::Medium)
#define logd log.lout(Logger::Level::Minor)

ostream& operator<<(ostream& os, Triconnectivity::CompType t) {
	switch (t) {
	case ogdf::Triconnectivity::CompType::triconnected:
		return os << "R";
	case ogdf::Triconnectivity::CompType::polygon:
		return os << "S";
	case ogdf::Triconnectivity::CompType::bond:
		return os << "P";
	default:
		return os << "?";
	}
}

void SimpleSPQRTree::init() {
	if (!skel_array.empty()) {
		return;
	}
	OGDF_ASSERT(!GC.empty());

	// remove all parallel edges, Triconnectivity doesn't seem to like those
	par_replacement.init(GC);
	SListPure<edge> edges;
	EdgeArray<int> minIndex(GC), maxIndex(GC);
	parallelFreeSortUndirected(GC, edges, minIndex, maxIndex);
	edge main = nullptr;
	int edgecnt = GC.numberOfEdges();
	for (edge e : edges) {
		if (main != nullptr && minIndex[e] == minIndex[main] && maxIndex[e] == maxIndex[main]) {
			OGDF_ASSERT(main != nullptr);
			par_replacement[main].pushBack(GC.original(e));
			GC.delEdge(e);
		} else {
			main = e;
		}
	}
	logm << "Removing parallels reduced edges from " << edgecnt << " to " << GC.numberOfEdges()
		 << endl;
	edgecnt = GC.numberOfEdges();

	if (GC.numberOfNodes() == 2) {
		OGDF_ASSERT(GC.numberOfEdges() == 1);
		skels.init(GC, nullptr);
		twins.init(GC, nullptr);
		skel_array.resize(1, nullptr);
		skel_array[0] = new OverlappingGraphCopy(GC_skels);
		OverlappingGraphCopy& skel = *skel_array[0];
		skel.newNode(GC.firstNode());
		skel.newNode(GC.lastNode());
		skel.newEdge(GC.firstEdge());
		skels[GC.firstEdge()] = &skel;
		return;
	}

	// compute Triconnectivity information, will add virtual edges to GC
	int virtual_edges = GC.numberOfEdges(), real_edges = 0;
	Triconnectivity tri(GC, true);
	// OGDF_ASSERT(tri.checkComp());
	virtual_edges = GC.numberOfEdges() - virtual_edges;

	// derive structural information: pairing of twin edges and embeddings of rigids
	skels.init(GC, nullptr);
	twins.init(GC, nullptr);
	skel_array.resize(tri.m_numComp, nullptr);
	logm << "Got " << tri.m_numComp << " components, " << virtual_edges << " virtual edges" << endl;
	for (int i = 0; i < tri.m_numComp; ++i) {
		Comp& comp = tri.m_component[i];
		skel_array[i] = new OverlappingGraphCopy(GC_skels);
		OverlappingGraphCopy& skel = *skel_array[i];
		auto& l = logm << "Component " << i << " type " << comp.m_type << " got skel " << &skel
					   << " with " << comp.m_edges.size() << " edges";
		if (log.is_lout(Logger::Level::Minor)) {
			for (edge e : comp.m_edges) {
				l << " " << e << " (" << e->index() << ")";
			}
		}
		l << endl;
		if (comp.m_edges.empty()) {
			continue; // S- or P-node was merged with a neighbor during Triconnectivity construction
		}
		Logger::Indent _(log);
		for (edge e : comp.m_edges) {
			if (skels[e] == nullptr) {
				skels[e] = &skel;
				real_edges++;
			} else {
				logd << "Edge " << e << " (" << e->index() << ") is virtual" << endl;
				OGDF_ASSERT(GC.original(e) == nullptr); // only for virtual edges
				OGDF_ASSERT(twins[e] == nullptr);
				twins[e] = &skel;
				virtual_edges--;
				real_edges--;
			}

			if (skel.copy(e->source()) == nullptr) {
				skel.newNode(e->source());
			}
			if (skel.copy(e->target()) == nullptr) {
				skel.newNode(e->target());
			}
			skel.newEdge(e);
		}
		OGDF_ASSERT(comp.m_edges.size() == skel.numberOfEdges());
		if (comp.m_type == Triconnectivity::CompType::triconnected) {
			OGDF_ASSERT(isTriconnected(skel));
			OGDF_ASSERT(isRNode(skel));
			planar = planarEmbedPlanarGraph(skel);
			if (!planar) {
				break;
			}
		} else if (comp.m_type == Triconnectivity::CompType::bond) {
			OGDF_ASSERT(isPNode(skel));
		} else {
			OGDF_ASSERT(comp.m_type == Triconnectivity::CompType::polygon);
			OGDF_ASSERT(isSNode(skel));
		}
	}
	OGDF_ASSERT(!planar || virtual_edges == 0);
	OGDF_ASSERT(!planar || real_edges == edgecnt);
}

OverlappingGraphCopy* SimpleSPQRTree::getNonSSkel(node GC_n) const {
	OGDF_ASSERT(GC_n->degree() > 2);
	OverlappingGraphCopy* skel = skels[GC_n->firstAdj()->theEdge()];
	if (isSNode(*skel)) {
		node skel_n = skel->copy(GC_n);
		OGDF_ASSERT(skel_n->degree() == 2);
		OverlappingGraphCopy* twinF = getTwinSkel(skel, skel_n->firstAdj()->theEdge());
		OverlappingGraphCopy* twinL = getTwinSkel(skel, skel_n->lastAdj()->theEdge());
		OGDF_ASSERT((twinF == nullptr) != (twinL == nullptr));
		skel = twinF == nullptr ? twinL : twinF;
	}
	OGDF_ASSERT(!isSNode(*skel));
	return skel;
}

OverlappingGraphCopy* SimpleSPQRTree::getTwinSkel(OverlappingGraphCopy* skel, edge skel_e) const {
	return getTwinSkel_GC(skel, skel->original(skel_e));
}

OverlappingGraphCopy* SimpleSPQRTree::getTwinSkel_GC(OverlappingGraphCopy* skel, edge GC_e) const {
	if (twins[GC_e] == nullptr) {
		OGDF_ASSERT(GC.original(GC_e) != nullptr);
		return nullptr;
	} else {
		OGDF_ASSERT(GC.original(GC_e) == nullptr);
	}
	OverlappingGraphCopy* twin;
	if (skels[GC_e] == skel) {
		twin = twins[GC_e];
	} else {
		OGDF_ASSERT(twins[GC_e] == skel);
		twin = skels[GC_e];
	}
	OGDF_ASSERT(twin != skel);
	return twin;
}

#define log spqr.log

NodeSSPQRRotation::NodeSSPQRRotation(const SimpleSPQRTree& spqr, node n) : spqr(spqr) {
	m_G = &spqr.GC.original();
	m_n = n;
	incidentEdgeForLeaf.init(*this, nullptr);
	graphNodeForInnerNode.init(*this, nullptr);
	OGDF_ASSERT(m_n->graphOf() == m_G);
	OGDF_ASSERT(m_n->degree() >= 3);

	node GC_n = spqr.GC.copy(n);
	OverlappingGraphCopy* skel =
			GC_n->degree() > 2 ? spqr.getNonSSkel(GC_n) : spqr.skels[GC_n->firstAdj()->theEdge()];
	node skel_n = skel->copy(GC_n);
	logm << "Generating NodeSSPQRRotation for " << "G_n " << n->index() << " °" << n->degree()
		 << ", " << "GC_n " << GC_n->index() << " °" << GC_n->degree() << ", " << "skel_n "
		 << skel_n->index() << " °" << skel_n->degree() << endl;
	logd << "First skel is " << skel << " (" << skel->numberOfNodes() << "n, "
		 << skel->numberOfEdges() << "e) " << "of adj " << GC_n->firstAdj() << " ("
		 << GC_n->firstAdj()->theEdge()->index() << ")" << endl;
	Logger::Indent _(log);


	if (GC_n->degree() == 1) {
		OGDF_ASSERT(skel_n->degree() == 1);
		process(skel_n->firstAdj(), *skel, nullptr);
	} else if (GC_n->degree() == 2) {
		OGDF_ASSERT(isSNode(*skel));
		OGDF_ASSERT(skel_n->degree() == 2);
		PCNode* n1 = process(skel_n->firstAdj(), *skel, nullptr);
		PCNode* n2 = process(skel_n->lastAdj(), *skel, nullptr);
		if (n1->isLeaf()) {
			swap(n1, n2);
		}
		logm << "Degree-2 S-Node Shortcut with Nodes " << n1 << " and " << n2 << "!" << endl;
		OGDF_ASSERT(!n1->isLeaf()); // otherwise m_n would be deg 2
		n1->appendChild(n2);
		setRoot(n1);
	} else {
		PCNode* root;
		if (isPNode(*skel)) {
			root = newNode(PCNodeType::PNode);
			graphNodeForInnerNode[root] =
					spqr.GC.original(skel->original(skel_n->firstAdj()->twinNode()));
		} else {
			root = newNode(PCNodeType::CNode);
		}
		for (adjEntry adj : skel_n->adjEntries) {
			process(adj, *skel, root);
		}
	}
	OGDF_ASSERT(getLeafCount() == n->degree());
	OGDF_ASSERT(checkValid());
}

PCNode* NodeSSPQRRotation::process(adjEntry skel_adj, OverlappingGraphCopy& skel, PCNode* parent) {
	adjEntry GC_adj = getAdjInOrig(&skel, skel_adj);
	OGDF_ASSERT(spqr.GC.original(GC_adj->theNode()) == m_n);
	edge G_edge = spqr.GC.original(GC_adj->theEdge());
	OverlappingGraphCopy* twin = spqr.getTwinSkel_GC(&skel, GC_adj->theEdge());
	if (G_edge != nullptr) {
		OGDF_ASSERT(twin == nullptr);
		adjEntry G_adj = G_edge->getAdj(spqr.GC.original(GC_adj->theNode()));
		logd << "Real edge " << "skel_adj " << skel_adj << " (" << skel_adj->theEdge()->index()
			 << "), " << "GC_adj " << GC_adj << " (" << GC_adj->theEdge()->index() << "), "
			 << "G_adj " << G_adj << " (" << G_adj->theEdge()->index() << ") " << "with "
			 << spqr.par_replacement[GC_adj].size() << " parallels" << endl;
		OGDF_ASSERT(spqr.twins[GC_adj] == nullptr);
		if (!spqr.par_replacement[GC_adj].empty()) {
			if (!isPNode(skel) || skel.numberOfEdges() == 1) { // no P-node
				logd << "\tCreating parent for parallels, as current skel is no P-node" << endl;
				parent = newNode(pc_tree::PCNodeType::PNode, parent);
				graphNodeForInnerNode[parent] = spqr.GC.original(GC_adj->twinNode());
			}
			auto& l = logd << "\tParallel edges:";
			for (edge e : spqr.par_replacement[GC_adj]) {
				l << " " << e << " (" << e->index() << ")";
				OGDF_ASSERT(e->isIncident(m_n));
				incidentEdgeForLeaf[newNode(pc_tree::PCNodeType::Leaf, parent)] = e;
			}
			l << endl;
		}
		pc_tree::PCNode* l = newNode(pc_tree::PCNodeType::Leaf, parent);
		incidentEdgeForLeaf[l] = G_adj->theEdge();
		if (parent == nullptr) {
			return l;
		} else {
			return parent;
		}

	} else {
		OGDF_ASSERT(twin != nullptr);
		adjEntry twin_adj = getAdjInSkel(twin, GC_adj);

		logd << (isSNode(*twin) ? "Inlining " : "") << "Virtual edge " << "skel_adj " << skel_adj
			 << " (" << skel_adj->theEdge()->index() << "), " << "GC_adj " << GC_adj << " ("
			 << GC_adj->theEdge()->index() << ") " << "corresponds to skel " << twin << " ("
			 << twin->numberOfNodes() << "n, " << twin->numberOfEdges() << "e)" << endl;
		if (isSNode(*twin)) { // S-node
			OGDF_ASSERT(twin_adj->theNode()->degree() == 2);
			return process(twin_adj->cyclicSucc(), *twin, parent);
		} else { // P- or R-node
			Logger::Indent _(log);
			PCNode* n;
			if (isPNode(*twin)) {
				n = newNode(PCNodeType::PNode, parent);
				graphNodeForInnerNode[n] = spqr.GC.original(GC_adj->twinNode());
			} else {
				n = newNode(PCNodeType::CNode, parent);
			}
			adjEntry adj = twin_adj->cyclicSucc();
			while (adj != twin_adj) {
				process(adj, *twin, n);
				adj = adj->cyclicSucc();
			}
			OGDF_ASSERT(n->getDegree() >= twin_adj->theNode()->degree());
			return n;
		}
	}
}

void NodeSSPQRRotation::mapPartnerEdges() {
	if (!isTrivial() || knowsPartnerEdges()) {
		return;
	}
	bundleEdgesForLeaf.init(*this);

	node partner = getTrivialPartnerPole();
	node GC_n = spqr.GC.copy(m_n);
	OverlappingGraphCopy* skel;
	if (GC_n->degree() > 2) {
		skel = spqr.getNonSSkel(GC_n);
		OGDF_ASSERT(isPNode(*skel));
	} else {
		skel = spqr.skels[GC_n->firstAdj()->theEdge()];
		OGDF_ASSERT(isSNode(*skel) || skel->numberOfEdges() == 1);
	}
	node skel_n = skel->copy(GC_n);
	int i = 0, deg = 0;
	for (adjEntry adj : skel_n->adjEntries) {
		PCNode* leaf = getLeaves()[i];
		edge GC_e = skel->original(adj->theEdge());
		edge orig_e = spqr.GC.original(GC_e);
		if (orig_e != nullptr) {
			if (GC_n->degree() == 2 && spqr.par_replacement[GC_e].empty()) {
				// very special case: we're in an S-node just due to the removal of parallel edges
				// the incident edge of the S-node that does not have parallels actually corresponds to
				// all edges of the partner pole that are not part of the parallel and in the current block
				for (adjEntry adj : partner->adjEntries) {
					if (adj->twinNode() == m_n) {
						continue;
					}
					edge GC_e = spqr.GC.copy(adj->theEdge());
					if (GC_e == nullptr) {
						continue;
					}
					for (edge e : spqr.par_replacement[GC_e]) {
						bundleEdgesForLeaf[leaf].pushBack(e);
					}
					bundleEdgesForLeaf[leaf].pushBack(adj->theEdge());
				}
				OGDF_ASSERT(getIncidentEdgeForLeaf(leaf) == orig_e);
			} else {
				for (edge e : spqr.par_replacement[GC_e]) {
					OGDF_ASSERT(getIncidentEdgeForLeaf(leaf) == e);
					bundleEdgesForLeaf[leaf].pushBack(e);
					deg++;
					i++;
					leaf = getLeaves()[i];
				}
				OGDF_ASSERT(getIncidentEdgeForLeaf(leaf) == orig_e);
				bundleEdgesForLeaf[leaf].pushBack(orig_e);
			}
		} else {
#ifdef OGDF_DEBUG
			{
				// check that the leaf and skel_n->adjEntry orders match
				List<edge> out;
				getIncidentRealEdgesInSubtree(adj, *skel, out);
				OGDF_ASSERT(out.size() == 1);
				OGDF_ASSERT(incidentEdgeForLeaf[leaf] == out.front());
			}
#endif
			getIncidentRealEdgesInSubtree(adj->twin(), *skel, bundleEdgesForLeaf[leaf]);
		}
		for (edge e : bundleEdgesForLeaf[leaf]) {
			OGDF_ASSERT(e->graphOf() == m_G);
			OGDF_ASSERT(e->isIncident(partner));
		}
		deg += bundleEdgesForLeaf[leaf].size();
		i++;
	}
#ifdef OGDF_DEBUG
	OGDF_ASSERT(i == getLeafCount());
	int exp_deg = 0; // get all edges of partner in this block
	for (adjEntry adj : partner->adjEntries) {
		edge GC_e = spqr.GC.copy(adj->theEdge());
		if (GC_e != nullptr) {
			exp_deg += 1 + spqr.par_replacement[GC_e].size();
		}
	}
	OGDF_ASSERT(exp_deg <= partner->degree());
	OGDF_ASSERT(deg == exp_deg);
#endif
}

void NodeSSPQRRotation::getIncidentRealEdgesInSubtree(adjEntry skel_adj, OverlappingGraphCopy& skel,
		List<edge>& out) {
	adjEntry GC_adj = getAdjInOrig(&skel, skel_adj);
	edge G_edge = spqr.GC.original(GC_adj->theEdge());
	OverlappingGraphCopy* twin = spqr.getTwinSkel_GC(&skel, GC_adj->theEdge());
	if (G_edge != nullptr) {
		for (edge e : spqr.par_replacement[GC_adj]) {
			out.pushBack(e);
		}
		out.pushBack(G_edge);
	} else {
		OGDF_ASSERT(twin != nullptr);
		adjEntry twin_adj = getAdjInSkel(twin, GC_adj);
		if (isSNode(*twin)) {
			OGDF_ASSERT(twin_adj->theNode()->degree() == 2);
			return getIncidentRealEdgesInSubtree(twin_adj->cyclicSucc(), *twin, out);
		} else {
			adjEntry adj = twin_adj->cyclicSucc();
			while (adj != twin_adj) {
				getIncidentRealEdgesInSubtree(adj, *twin, out);
				adj = adj->cyclicSucc();
			}
		}
	}
}
