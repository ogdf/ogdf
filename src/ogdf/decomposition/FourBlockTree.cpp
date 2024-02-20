#include <ogdf/decomposition/FourBlockTree.h>

template<typename _C, typename _L, typename _T = typename _C::value_type>
std::vector<_T> countingSort(const _C& input, _L key) {
    std::vector<size_t> counts;
    std::vector<std::vector<_T>> buckets;
    for (const auto& x : input) {
        const size_t k = key(x);
        if (k >= counts.size()) {
            counts.resize(k + 1, 0);
        }
        ++counts[k];
    }
    for (size_t i = 1; i < counts.size(); ++i) {
        counts[i] += counts[i - 1];
    }
    std::vector<_T> res(input.size());
    for (auto it = input.rbegin(); it != input.rend(); ++it) {
        res[--counts[key(*it)]] = *it;
    }
    return res;
}

void ogdf::FourBlockTreeBuilder::populateIndices() {
    for (const auto v : m_g.nodes) {
        size_t i = 0;
        for (const auto a : v->adjEntries) {
            m_indices[a] = i++;
        }
    }
}

void ogdf::FourBlockTreeBuilder::populateSepTriangles() {
    /* find all triangles [Chiba, Nishizeki] */

    /* sort vertices by degree using CountingSort */
    const auto byDescDegree =
            countingSort(m_g.nodes, [&](ogdf::node v) { return m_g.numberOfNodes() - v->degree(); });

    ogdf::NodeArray<bool> marked(m_g, false);
    ogdf::NodeArray<bool> deleted(m_g, false);
    ogdf::NodeArray<ogdf::adjEntry> a_vw(m_g, nullptr); // a_vw[w]->theEdge() = {v,w} for all w adj to v

    /* for all vertices ordered by non-ascending degree */
    for (const auto v : byDescDegree) {
        /* mark all neighbors of v */
        for (const auto a_vw_ : v->adjEntries) {
            const auto w = a_vw_->twinNode();
            marked[w] = true;
            a_vw[w] = a_vw_;
        }

        for (const auto a_vu : v->adjEntries) // for each marked vertex u
        {
            const auto u = a_vu->twinNode();
            if (deleted[u]) {
                continue;
            }
            for (const auto a_uw : u->adjEntries) // for each vertex w adj to u
            {
                const auto w = a_uw->twinNode();
                if (deleted[w]) {
                    continue;
                }
                if (marked[w]) {
                    /* {u,v,w} is a triangle */

                    const auto ix_vu = m_indices[a_vu];
                    const auto ix_vw = m_indices[a_vw[w]];
                    if ((v->degree() + ix_vu - ix_vw + 1) % v->degree() > 2) { // <=> |ix_vu - ix_vw| > 1
                        /* {u,v,w} is separating */
                        m_sepTriangles.push_back(
                                {a_vu->theEdge(), a_vw[w]->theEdge(), a_uw->theEdge()});
                    }
                }
            }

            marked[u] = false;
        }

        deleted[v] = true;
    }
}

void ogdf::FourBlockTreeBuilder::firstDfs() {
    ogdf::NodeArray<bool> visited(m_g, false);
    ogdf::EdgeArray<bool> traversed(m_g, false);
    std::vector<ogdf::adjEntry> stack = {m_externalFace};
    visited[m_root] = traversed[m_externalFace->theEdge()] = true;
    m_parentEdge[m_root] = m_externalFace;

    /* if stack.back() = nullptr, we are done with the most recently visited node on level stack.size() */
    /* otherwise we are always in the subtree rooted at stack.back()->theNode(), traversing stack.back()->theEdge() */
    /* stack[m_depth[v]]->theNode() = v forall v on currently active path from root */

    while (!stack.empty()) {
        const auto adj = stack.back();
        if (adj) {
            const auto w = adj->twinNode();
            const auto e = adj->theEdge();
            if (!visited[w]) {
                /* e is tree edge */
                m_depth[w] = stack.size(); // = m_depth[v] + 1
                m_parentEdge[w] = adj->twin();
                visited[w] = true;
                traversed[e] = true;
                m_isTreeEdge[e] = true;
                if (!adj->isSource()) {
                    m_g.reverseEdge(e);
                }
                stack.push_back(w->firstAdj());
            } else if (!traversed[e]) {
                /* e is tree edge */
                traversed[e] = true;
                m_isTreeEdge[e] = false;
                m_lowpoint[e] = w;
                if (!adj->isSource()) {
                    m_g.reverseEdge(e);
                }

                /* indices at w */
                const auto ix_e = m_indices[adj->twin()];
                const auto ix_p = m_indices[m_parentEdge[w]];
                const auto ix_d = m_indices[stack[m_depth[w]]];

                const size_t ang_e_d = (w->degree() + ix_d - ix_e) % w->degree();
                const size_t ang_e_p = (w->degree() + ix_p - ix_e) % w->degree();
                if (ang_e_d < ang_e_p) {
                    /* in clockwise order around w: e, active tree edge = stack[m_depth[w]], m_parentEdge[w] */
                    m_returnSide[e] = RIGHT;
                    m_angularDistance[e] = ang_e_d;
                } else {
                    m_returnSide[e] = LEFT;
                    m_angularDistance[e] = w->degree() - ang_e_d;
                }

                /* go to next edge of v */
                stack.back() = adj->succ();
            } else {
                /* go to next edge of v */
                stack.back() = adj->succ();
            }
        } else {
            stack.pop_back();

            if (!stack.empty()) {
                /* w is done now, m_parentEdge[w] inherits properties of outgoing edges of w */
                const auto w = stack.back()->twinNode();
                const auto p = stack.back()->theEdge(); // = m_parentEdge[w]->theEdge()

                m_lowpoint[p] = w;
                m_returnSide[p] = 0;
                m_angularDistance[p] = 0;
                for (const auto a : w->adjEntries) {
                    /* only consider outgoing edges of w */
                    if (!a->isSource()) {
                        continue;
                    }
                    const auto e = a->theEdge();

                    if (m_depth[m_lowpoint[e]] < m_depth[m_lowpoint[p]]) {
                        m_lowpoint[p] = m_lowpoint[e];
                        m_returnSide[p] = m_returnSide[e];
                        m_angularDistance[p] = m_angularDistance[e];
                    } else if (m_depth[m_lowpoint[e]] == m_depth[m_lowpoint[p]]) {
                        m_returnSide[p] |= m_returnSide[e];
                        m_angularDistance[p] = std::max(m_angularDistance[e], m_angularDistance[e]);
                    }
                }

                /* go to next edge of v */
                stack.back() = stack.back()->succ();
            }
        }
    }
}

void ogdf::FourBlockTreeBuilder::orderTriangles() {
    /* sort edges e using Counting-/RadixSort on */
    /* 1. e->source() */
    /* 2. -m_depth[m_lowpoint[e]] */
    /* 3. m_angularDistance[e] */
    const auto byAngularDistance =
            countingSort(m_g.edges, [&](ogdf::edge e) { return m_angularDistance[e]; });
    const auto byDepthLowpoint = countingSort(byAngularDistance,
            [&](ogdf::edge e) { return m_g.numberOfNodes() - m_depth[m_lowpoint[e]]; });
    ogdf::NodeArray<std::vector<ogdf::edge>> edge_order(m_g);
    for (const auto e : byDepthLowpoint) {
        edge_order[e->source()].push_back(e);
    }

    std::vector<unsigned char> untraversedEdges(m_sepTriangles.size(),
            3); // number of edges for each triangle that have not yet been traversed
    ogdf::EdgeArray<std::vector<size_t>> triangleIndices(
            m_g); // indices into untraversedEdges and m_sepTriangles
    for (size_t i = 0; i < m_sepTriangles.size(); ++i) {
        const auto& t = m_sepTriangles[i];
        for (const auto e : t) {
            triangleIndices[e].push_back(i);
        }
    }

    /* find order of triangles */
    /* indices into m_sepTriangles are also indices into these two arrays */
    std::vector<size_t> internalAngle(m_sepTriangles.size(), 0);
    std::vector<size_t> foundWithEdge(m_sepTriangles.size(), 0);

    /* second DFS */
    std::vector<std::pair<ogdf::node, std::vector<ogdf::edge>::const_iterator>> stack = {
            {m_root, edge_order[m_root].cbegin()}};
    size_t edge_id = 0;
    while (!stack.empty()) {
        auto& [v, it] = stack.back();
        if (it == edge_order[v].cend()) {
            stack.pop_back();
        } else {
            ++edge_id;
            const auto e = *it++;
            const auto w = e->target();
            if (m_isTreeEdge[e]) {
                stack.emplace_back(w, edge_order[w].cbegin());
            }

            std::vector<size_t> completeTriangles; // elements are indices into m_sepTriangles
            for (const auto i : triangleIndices[e]) {
                if (--untraversedEdges[i] == 0) {
                    completeTriangles.push_back(i);
                    foundWithEdge[i] = edge_id;

                    /* make e first edge in triangle */
                    auto& t = m_sepTriangles[i];
                    if (e == t[1]) {
                        std::swap(t[0], t[1]);
                    } else if (e == t[2]) {
                        std::swap(t[0], t[2]);
                    }

                    /* "direct" t clockwise */
                    if (t[1]->isIncident(w) ^ (m_returnSide[e] == LEFT)) {
                        std::swap(t[1], t[2]);
                    }
                }
            }
            if (completeTriangles.size() > 1) {
                /* e must be a back edge, and all triangles are on the same side of e */
                const bool trianglesRight = m_returnSide[e] == LEFT;

                const auto ix_e = m_indices[e->adjTarget()];

                for (const auto i : completeTriangles) {
                    const auto& t = m_sepTriangles[i];
                    const auto a = [&]() {
                        /* find adjEntry a on t with a->theNode() == w and a->twinNode() != v */
                        for (const auto side : t) {
                            for (const auto adj : {side->adjSource(), side->adjTarget()}) {
                                if (adj->theNode() == w && adj->twinNode() != v) {
                                    return adj;
                                }
                            }
                        }
                        /* unreachable */
                        std::terminate();
                    }();
                    const auto ang_a_e = (w->degree() + ix_e - m_indices[a]) % w->degree();
                    if (trianglesRight) {
                        internalAngle[i] = ang_a_e;
                    } else {
                        internalAngle[i] = w->degree() - ang_a_e;
                    }
                }
            }
        }
    }

    std::vector<size_t> tmp;
    tmp.reserve(m_sepTriangles.size());
    for (size_t i = 0; i < m_sepTriangles.size(); ++i) {
        tmp.push_back(i);
    }
    const auto byInternalAngle = countingSort(tmp, [&](size_t i) { return internalAngle[i]; });
    const auto byFoundWithEdge =
            countingSort(byInternalAngle, [&](size_t i) { return foundWithEdge[i]; });
    std::vector<triangle_t> res;
    for (const size_t i : byFoundWithEdge) {
        res.push_back(m_sepTriangles[i]);
    }

    m_sepTriangles = std::move(res);
}

ogdf::FourBlockTree ogdf::FourBlockTreeBuilder::buildTree() {
    std::vector<std::unique_ptr<ogdf::FourBlockTree>> blocks;

    /* used to set FourBlockTree.parent */
    ogdf::AdjEntryArray<ogdf::FourBlockTree*> parentInv(m_g, nullptr);

    /* isInner[v] => v is part of an inner block that is about to be copied from m_g */
    /* initialized once to avoid linear cost for each block */
    ogdf::NodeArray<bool> isInner(m_g, false);

    /* the ID the edge will have in b */
    /* initialized once to avoid linear cost for each block */
    ogdf::EdgeArray<size_t> edgeIdsForB(m_g, 0);

    for (const auto& t : m_sepTriangles) {
        /* u,v,w in clockwise order around t */
        const auto [uv, vw, wu] = t;
        const auto v = uv->commonNode(vw);
        const auto w = vw->commonNode(wu);
        const auto u = wu->commonNode(uv);

        /* create new block b */
        blocks.push_back(std::make_unique<ogdf::FourBlockTree>());
        auto& b = *blocks.back();

        /* split m_g at t */
        const auto u_ = m_g.newNode();
        const auto v_ = m_g.newNode();
        const auto w_ = m_g.newNode();
        const auto vw_ = m_g.newEdge(v_, w_);
        const auto wu_ = m_g.newEdge(w_, u_);
        const auto uv_ = m_g.newEdge(u_, v_);
        m_vertexIds[u_] = m_vertexIds[u];
        m_vertexIds[v_] = m_vertexIds[v];
        m_vertexIds[w_] = m_vertexIds[w];

        /* move edges incident to v */
        auto adj = vw->getAdj(v)->cyclicSucc();
        while (adj != uv->getAdj(v)) {
            const auto tmp = adj->cyclicSucc();
            if (adj->isSource()) {
                m_g.moveSource(adj->theEdge(), uv_->adjTarget(), ogdf::Direction::before);
            } else {
                m_g.moveTarget(adj->theEdge(), uv_->adjTarget(), ogdf::Direction::before);
            }
            adj = tmp;
        }

        /* move edges incident to w */
        adj = wu->getAdj(w)->cyclicSucc();
        while (adj != vw->getAdj(w)) {
            const auto tmp = adj->cyclicSucc();
            if (adj->isSource()) {
                m_g.moveSource(adj->theEdge(), vw_->adjTarget(), ogdf::Direction::before);
            } else {
                m_g.moveTarget(adj->theEdge(), vw_->adjTarget(), ogdf::Direction::before);
            }
            adj = tmp;
        }

        /* move edges incident to u */
        adj = uv->getAdj(u)->cyclicSucc();
        while (adj != wu->getAdj(u)) {
            const auto tmp = adj->cyclicSucc();
            if (adj->isSource()) {
                m_g.moveSource(adj->theEdge(), wu_->adjTarget(), ogdf::Direction::before);
            } else {
                m_g.moveTarget(adj->theEdge(), wu_->adjTarget(), ogdf::Direction::before);
            }
            adj = tmp;
        }

        /* set parent pointers to b */
        for (const auto a : {uv->getAdj(v), vw->getAdj(w), wu->getAdj(u)}) {
            auto& c = parentInv[a];
            if (c) {
                c->parent = &b;
            }
            c = &b;
        }

        /* list all vertices to be copied to b */
        std::vector<ogdf::node> innerVertices = {u_, v_, w_}; // vertices in inner block
        isInner[u_] = isInner[v_] = isInner[w_] = true;
        size_t b_m = 0; // number of edges in inner block
        size_t nextEdgeIdForB = 0;
        for (size_t i = 0; i < innerVertices.size(); ++i) {
            for (const auto a : innerVertices[i]->adjEntries) {
                /* count edges */
                if (a->isSource()) {
                    ++b_m;
                    edgeIdsForB[a] = nextEdgeIdForB++;
                }

                /* add neighbor x if not already added */
                const auto x = a->twinNode();
                if (!isInner[x]) {
                    isInner[x] = true;
                    innerVertices.push_back(x);
                }
            }
        }

        /* copy innerVertices with their edges to b */
        b.vertexIds.init(b.g);
        b.parent = nullptr;
        b.parentFace = nullptr;
        const auto bDummy = b.g.newNode();
        std::vector<ogdf::edge> bEdges;
        for (size_t i = 0; i < b_m; ++i) {
            bEdges.push_back(b.g.newEdge(bDummy, bDummy));
        }
        for (const auto x : innerVertices) {
            /* copy vertex */
            const auto b_v = b.g.newNode();
            b.vertexIds[b_v] = m_vertexIds[x];

            for (const auto a : x->adjEntries) {
                /* copy adjEntry */
                const auto b_e = bEdges[edgeIdsForB[a]];
                if (a->isSource()) {
                    b.g.moveSource(b_e, b_v);
                } else {
                    b.g.moveTarget(b_e, b_v);
                }

                /* set c.parent and c.parentFace for children c of b */
                const auto c = parentInv[a];
                if (c) {
                    const auto b_a = a->isSource() ? b_e->adjSource() : b_e->adjTarget();
                    c->parent = &b;
                    c->parentFace = b_a;
                }
            }

            if (x == v_) {
                b.externalFace = b_v->firstAdj();
            }
        }
        b.g.delNode(bDummy);
    }

    /* root of 4-block tree */
    ogdf::FourBlockTree res;
    {
        /* list all vertices to be copied to res */
        std::vector<ogdf::node> innerVertices = {m_root}; // vertices in outermost block
        isInner[m_root] = true;
        size_t b_m = 0; // number of edges in outermost block
        size_t nextEdgeIdForB = 0;
        for (size_t i = 0; i < innerVertices.size(); ++i) {
            for (const auto a : innerVertices[i]->adjEntries) {
                /* count edges */
                if (a->isSource()) {
                    ++b_m;
                    edgeIdsForB[a] = nextEdgeIdForB++;
                }

                /* add neighbor w if not already added */
                const auto w = a->twinNode();
                if (!isInner[w]) {
                    isInner[w] = true;
                    innerVertices.push_back(w);
                }
            }
        }

        /* copy innerVertices with their edges to res */
        res.vertexIds.init(res.g);
        res.parent = nullptr;
        res.parentFace = nullptr;
        const auto bDummy = res.g.newNode();
        std::vector<ogdf::edge> bEdges;
        for (size_t i = 0; i < b_m; ++i) {
            bEdges.push_back(res.g.newEdge(bDummy, bDummy));
        }
        for (const auto v : innerVertices) {
            /* copy vertex */
            const auto b_v = res.g.newNode();
            res.vertexIds[b_v] = m_vertexIds[v];

            for (const auto a : v->adjEntries) {
                /* copy adjEntry */
                const auto b_e = bEdges[edgeIdsForB[a]];
                if (a->isSource()) {
                    res.g.moveSource(b_e, b_v);
                } else {
                    res.g.moveTarget(b_e, b_v);
                }

                /* set res.externalFace */
                if (a == m_externalFace) {
                    res.externalFace = a->isSource() ? b_e->adjSource() : b_e->adjTarget();
                }

                /* set c.parent and c.parentFace for children c of res */
                const auto c = parentInv[a];
                if (c) {
                    const auto b_a = a->isSource() ? b_e->adjSource() : b_e->adjTarget();
                    c->parent = &res;
                    c->parentFace = b_a;
                }
            }
        }
        res.g.delNode(bDummy);
    }

    /* move inner blocks to their parents */
    for (auto& b : blocks) {
        auto& p = *b->parent;
        p.children.push_back(std::move(b));
    }

    return res;
}

ogdf::FourBlockTreeBuilder::FourBlockTreeBuilder(ogdf::Graph& g, ogdf::NodeArray<size_t>& vertexIds,
        ogdf::adjEntry externalFace)
    : m_g(g)
    , m_vertexIds(vertexIds)
    , m_externalFace(externalFace)
    , m_root(externalFace->theNode())
    , m_indices(g, 0)
    , m_depth(g, 0)
    , m_parentEdge(g, nullptr)
    , m_isTreeEdge(g, false)
    , m_lowpoint(g, nullptr)
    , m_returnSide(g, 0)
    , m_angularDistance(g, 0) { }

ogdf::FourBlockTree ogdf::FourBlockTreeBuilder::call() {
    populateIndices();
    populateSepTriangles();
    firstDfs();
    orderTriangles();
    return buildTree();
}
