#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/cluster/sync_plan/utils/NodeSPQRRotation.h>

Logger NodeSPQRRotation::logger;

#define logd logger.lout(Logger::Level::Minor)
#define log logger.lout(Logger::Level::Medium)

ostream& operator<<(ostream& os, DynamicSPQRForest::TNodeType t) {
	switch (t) {
	case DynamicSPQRForest::TNodeType::RComp:
		return os << "R";
	case DynamicSPQRForest::TNodeType::SComp:
		return os << "S";
	case DynamicSPQRForest::TNodeType::PComp:
		return os << "P";
	default:
		return os << "?";
	}
}

RigidEmbedding::RigidEmbedding(Graph& G) : spqr(G, true), rigids(spqr.spqrTree(), nullptr) {
	for (node bc : spqr.bcTree().nodes) {
		if (spqr.typeOfBNode(bc) != BCTree::BNodeType::BComp) {
			continue;
		}
		if (spqr.numberOfEdges(bc) <= 3 || spqr.numberOfNodes(bc) <= 2) {
			continue;
		}
		spqr.createSPQR(bc);
		node root = spqr.spqrroot(bc);
		OGDF_ASSERT(root != nullptr);
	}

	for (node n : spqr.spqrTree().nodes) {
		if (spqr.typeOfTNode(n) != DynamicSPQRForest::TNodeType::RComp) {
			continue;
		}
		GraphCopySimple* skel = rigids[n] = new GraphCopySimple();
		skel->createEmpty(spqr.auxiliaryGraph());
		for (edge e : spqr.hEdgesSPQR(n)) {
			if (skel->copy(e->source()) == nullptr) {
				skel->newNode(e->source());
			}
			if (skel->copy(e->target()) == nullptr) {
				skel->newNode(e->target());
			}
			skel->newEdge(e);
		}
		OGDF_ASSERT(spqr.hEdgesSPQR(n).size() == skel->numberOfEdges());
		OGDF_ASSERT(isTriconnected(*skel));
		if (!planarEmbed(*skel)) {
			throw GraphNotPlanarException();
		}
	}
}

node NodeSPQRRotation::findSPQRApex(node n, bool clear) {
	// find subtree with edges incident to max_deg
	SList<node> todo;

	int reals = 0;
	for (adjEntry adj : n->adjEntries) {
		bool is_virt = spqr.twinEdge(adj->theEdge()) != nullptr;
		bool is_real = spqr.original(adj->theEdge()) != nullptr;
		node t = spqr.spqrproper(adj->theEdge());
		OGDF_ASSERT(t != nullptr);
		logd << "Processing adj " << adj << " of r"
			 << (is_real ? spqr.original(adj->theEdge())->index() : -1) << " v"
			 << (is_virt ? spqr.twinEdge(adj->theEdge())->index() : -1) << " edge "
			 << adj->theEdge()->index() << " in SPQR " << spqr.typeOfTNode(t) << "-node "
			 << t->index() << endl;
		OGDF_ASSERT(is_virt != is_real);
		if (is_virt) {
			continue;
		}
		reals++;

		if (clear) {
			if (highest_with_edges[t] != nullptr) {
				todo.pushBack(t);
			}
		} else {
			if (highest_with_edges[t] == nullptr) {
				highest_with_edges[t] = t;
				todo.pushBack(t);
			}
			edges[t].pushBack(adj);
		}
	}
	{
		node o = spqr.original(n);
		if (spqr.typeOfGNode(o) == BCTree::GNodeType::Normal) {
			OGDF_ASSERT(reals == o->degree());
		} else {
			OGDF_ASSERT(reals <= o->degree());
		}
	}

	// find highest node in interesting subtree and list its children
	node highest = nullptr;
	while (!todo.empty()) {
		node next = todo.popFrontRet();
		logd << spqr.typeOfTNode(next) << "-node " << next->index() << ", " << edges[next].size()
			 << " real edges and " << children[next].size() << " (<= °" << next->degree()
			 << ") children, pred is " << highest_with_edges[next]->index() << endl;
		if (todo.empty() && highest == nullptr) {
			logd << "Queue size 1!" << endl;
			highest = next;
			break;
		}

		node par = spqr.spqrParent(next);
		if (par == nullptr) {
			OGDF_ASSERT(highest == nullptr);
			logd << "Parent is null, setting highest" << endl;
			highest = next;
			continue;
		}

		if (clear) {
			logd << "Parent " << spqr.typeOfTNode(par) << "-node " << par->index() << " is "
				 << (highest_with_edges[par] != nullptr ? "un" : "already ") << "processed" << endl;
			edges[next].clear();
			children[next].clear();
			highest_with_edges[next] = nullptr;
			if (highest_with_edges[par] != nullptr) {
				todo.pushBack(par);
			}
		} else {
			logd << "Parent " << spqr.typeOfTNode(par) << "-node " << par->index() << " is "
				 << (highest_with_edges[par] == nullptr ? "un" : "already ") << "processed" << endl;
			children[par].pushBack(next);
			if (highest_with_edges[par] == nullptr) {
				highest_with_edges[par] = highest_with_edges[next];
				todo.pushBack(par);
			}
		}
	}

	if (clear) {
		edges[highest].clear();
		children[highest].clear();
		highest_with_edges[highest] = nullptr;
	}
	return highest_with_edges[highest];
}

pc_tree::PCNode* NodeSPQRRotation::addLeaf(pc_tree::PCNode* n, adjEntry adj) {
	OGDF_ASSERT(spqr.original(adj->theEdge()) != nullptr);
	pc_tree::PCNode* l = newNode(pc_tree::PCNodeType::Leaf, n);
	incidentEdgeForLeaf[l] = adj->theEdge();
	logd << "Created leaf " << l << " for adj " << adj->index() << " " << adj << endl;
	return l;
}

pc_tree::PCNode* NodeSPQRRotation::makePCNode(node t, node t_parent, pc_tree::PCNode* parent) {
	OGDF_ASSERT(highest_with_edges[t] != nullptr); // check that arrays haven't been cleared yet
	pc_tree::PCNode* n;
	log << spqr.typeOfTNode(t) << "-node " << t->index() << " containing "
		<< spqr.hEdgesSPQR(t).size() << " edges. Target node has " << edges[t].size()
		<< " real edges and " << children[t].size() << " (<= °" << t->degree() << ") children."
		<< endl;
	auto& l = logd << " ";
	for (edge e : spqr.hEdgesSPQR(t)) {
		l << " " << (spqr.twinEdge(e) == nullptr ? "r" : "v") << e->index() << " " << e;
	}
	l << endl;
	Logger::Indent _(logger);

	if (spqr.typeOfTNode(t) == DynamicSPQRForest::TNodeType::SComp) {
		if (parent == nullptr) {
			// if an S-node is the root, we create a (deg-2) P-node which we will remove later on
			OGDF_ASSERT((edges[t].size() + children[t].size()) == 2);
			OGDF_ASSERT(children[t].size() > 0);
			n = makePCNode(children[t].front(), t, nullptr);
			for (adjEntry adj : edges[t]) {
				addLeaf(n, adj);
			}
			for (node ct : children[t]) {
				if (ct == children[t].front()) {
					continue;
				}
				OGDF_ASSERT(ct != t_parent);
				makePCNode(ct, t, n);
			}
			log << "Root S-node replaced by first child: " << n << endl;
			return n;
		} else {
			OGDF_ASSERT((edges[t].size() + children[t].size()) == 1);
			n = parent;
		}
	} else if (spqr.typeOfTNode(t) == DynamicSPQRForest::TNodeType::RComp) {
		n = newNode(pc_tree::PCNodeType::CNode, parent);
	} else {
		OGDF_ASSERT(spqr.typeOfTNode(t) == DynamicSPQRForest::TNodeType::PComp);
		n = newNode(pc_tree::PCNodeType::PNode, parent);
		graphNodeForInnerNode[n] = spqr.hEdgesSPQR(t).front()->opposite(m_n);
		OGDF_ASSERT(graphNodeForInnerNode[n]->degree() >= (edges[t].size() + children[t].size()));
	}
	if (n == parent) {
		log << "Using parent node " << n << endl;
	} else {
		log << "Created node " << n << endl;
	}
	node gn = graphNodeForInnerNode[n];
	if (gn != nullptr) {
		node ggn = spqr.original(gn);
		log << "H-Graph node for P-node is node " << gn->index() << " of degree " << gn->degree()
			<< " actual G-Graph node will be " << ggn->index() << " of degree " << ggn->degree()
			<< endl;
	}

	if (spqr.typeOfTNode(t) == DynamicSPQRForest::TNodeType::RComp) {
		const GraphCopySimple& skel = *rigids[t];
		node alloc = skel.copy(m_n);
		for (adjEntry adj : alloc->adjEntries) {
			edge e = skel.original(adj->theEdge());
			edge et = spqr.twinEdge(e);
			if (et != nullptr) {
				node nt = spqr.spqrproper(et);
				if (nt == t_parent) {
					// order would be [c_n ... c_k] parent [c_k-1 ... c_0]
					// change to      [c_n ... c_k] [c_0 ... c_k-1] parent
					n->flip();
				} else {
					makePCNode(nt, t, n);
				}
			} else {
				addLeaf(n, e->getAdj(m_n));
			}
		}

		logd << "Children identified by tree were:" << endl;
		Logger::Indent _(logger);
		for (adjEntry adj : edges[t]) {
			logd << "Adj " << adj->index() << " " << adj << endl;
		}
		for (node ct : children[t]) {
			logd << "Child " << spqr.typeOfTNode(ct) << "-node " << ct->index() << endl;
		}
		OGDF_ASSERT(alloc->degree() == n->getDegree());

	} else {
		for (adjEntry adj : edges[t]) {
			addLeaf(n, adj);
		}
		for (node ct : children[t]) {
			OGDF_ASSERT(ct != t_parent);
			makePCNode(ct, t, n);
		}
	}

	if (n != parent) {
		log << "Result: " << n << endl;
		if (parent == nullptr) {
			OGDF_ASSERT(n->getDegree() == (edges[t].size() + children[t].size()));
		} else {
			OGDF_ASSERT(n->getDegree() == (edges[t].size() + children[t].size() + 1));
		}
	}
	return n;
}

void NodeSPQRRotation::mapPartnerEdges() {
	if (!isTrivial() || knowsPartnerEdges()) {
		return;
	}
	bundleEdgesForLeaf.init(*this);
	EdgeArray<PCNode*> mapping(*getGraph());
	generateLeafForIncidentEdgeMapping(mapping);
	NodeArray<PCNode*> snodes(spqr.spqrTree(), nullptr); // TODO reuse mappings?
	node pole = getTrivialPartnerPole();
	OGDF_ASSERT(pole != nullptr);
	OGDF_ASSERT(pole->graphOf() == m_G);

	node pnode = apex;
	if (spqr.typeOfTNode(pnode) == DynamicSPQRForest::TNodeType::SComp) {
		// if an S-node is the root, skip to its first child
		OGDF_ASSERT((edges[pnode].size() + children[pnode].size()) == 2);
		OGDF_ASSERT(children[pnode].size() > 0);
		pnode = children[pnode].front();
	}
	OGDF_ASSERT(spqr.typeOfTNode(pnode) == DynamicSPQRForest::TNodeType::PComp);
	OGDF_ASSERT(spqr.hEdgesSPQR(pnode).size() == getLeafCount());

	int reals = 0, reals2 = 0, bundle = 0;
	for (edge e : spqr.hEdgesSPQR(pnode)) {
		bool is_virt = spqr.twinEdge(e) != nullptr;
		bool is_real = spqr.original(e) != nullptr;
		OGDF_ASSERT(is_virt != is_real);
		if (is_real) {
			PCNode* l = mapping[e];
			OGDF_ASSERT(l != nullptr);
			OGDF_ASSERT(bundleEdgesForLeaf[l].empty());
			bundleEdgesForLeaf[l].pushBack(e);
			reals++;
			continue;
		}
		node snode = spqr.spqrNodeOf(spqr.twinEdge(e));
		OGDF_ASSERT(spqr.typeOfTNode(snode) == DynamicSPQRForest::TNodeType::SComp);
		OGDF_ASSERT(edges[snode].size() == 1);
		PCNode* l = mapping[edges[snode].front()];
		OGDF_ASSERT(l != nullptr);
		OGDF_ASSERT(bundleEdgesForLeaf[l].empty());
		OGDF_ASSERT(snodes[snode] == nullptr);
		snodes[snode] = l;
	}

	for (adjEntry adj : pole->adjEntries) {
		if (spqr.twinEdge(adj->theEdge()) != nullptr) {
			continue;
		}
		if (adj->twinNode() == m_n) {
			reals2++;
			continue;
		}
		// only process real edges that are not incident to first pole
		node tn = spqr.spqrproper(adj->theEdge());
		OGDF_ASSERT(tn != nullptr);
		while (snodes[tn] == nullptr) {
			if (spqr.spqrParent(tn) == nullptr) {
				// if we reached the root, the path continues up to the parent of the parallel
				tn = spqr.spqrParent(pnode);
				break;
			}
			tn = spqr.spqrParent(tn);
		}
		OGDF_ASSERT(tn != nullptr);
		OGDF_ASSERT(spqr.typeOfTNode(tn) == DynamicSPQRForest::TNodeType::SComp);
		OGDF_ASSERT(snodes[tn] != nullptr);
		bundleEdgesForLeaf[snodes[tn]].pushBack(adj->theEdge());
		bundle++;
	}
	OGDF_ASSERT(reals == reals2);
	if (spqr.typeOfGNode(spqr.original(pole)) == BCTree::GNodeType::Normal) {
		OGDF_ASSERT(reals + bundle == spqr.original(pole)->degree());
	} else {
		OGDF_ASSERT(reals + bundle <= spqr.original(pole)->degree());
	}
}
