/** \file
 * \brief Implementation of the random cluster-based k-spanner
 * algorithm of Baswana and Sen 2007.
 *
 * \author Finn Stutzenstein
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.md in the OGDF root directory for details.
 *
 * \par
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * Version 2 or 3 as published by the Free Software Foundation;
 * see the file LICENSE.txt included in the packaging of this file
 * for details.
 *
 * \par
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * \par
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, see
 * http://www.gnu.org/copyleft/gpl.html
 */

#pragma once

#include <ogdf/basic/NodeSet.h>
#include <ogdf/graphalg/SpannerIteratedWrapper.h>
#include <ogdf/graphalg/SpannerModule.h>

namespace ogdf {

/**
 * Randomized multiplicative spanner calculation by forming clusters.
 *
 * S. Baswana und S. Sen. A simple and linear time randomized algorithm for computing
 * sparse spanners in weighted graphs. Random Structures & Algorithms 30.4 (2007), S. 532-563.
 * doi: https://doi.org/10.1002/rsa.20130.
 *
 * Conditions for the graph:
 * - undirected
 * - weighted
 *
 * The stretch \f$s\f$ must satisfy: \f$s\in\{1,3,5,7,\ldots\}\f$.
 *
 * The preconditions can be checked with SpannerBaswanaSen::preconditionsOk.
 *
 * Calculates a \f$(2k-1)\f$-spanner with size \f$\mathcal{O}(kn^{1+1/k})\f$. There are
 * no guarantees for the lightness. The \e expected runtime is \f$\mathcal{O}(km)\f$.
 *
 * @ingroup ga-spanner
 */
template<typename TWeight>
class SpannerBaswanaSen : public SpannerModule<TWeight> {
	//! The copy of GA.constGraph() to work on. Edges will be removed in each iteration.
	GraphCopySimple m_G;
	//! the current cluster for each iteration in phase 1 and the final cluster from phase 1 which is used in phase 2
	NodeArray<node> m_cluster;

	EpsilonTest m_eps;

public:
	//! @copydoc ogdf::SpannerModule::preconditionsOk
	virtual bool preconditionsOk(const GraphAttributes& GA, double stretch,
			std::string& error) override {
		if (GA.directed()) {
			error = "The graph must be undirected";
			return false;
		}
		double integralPart;
		if (std::modf(stretch, &integralPart) != 0.0) {
			error = "The stretch is required to be an integer, not " + to_string(m_stretch);
			return false;
		}
		int intStretch = static_cast<int>(stretch);
		if (intStretch < 1) {
			error = "The stretch must be >= 1.0";
			return false;
		}
		if (intStretch % 2 == 0) {
			error = "The stretch must be odd";
			return false;
		}
		return true;
	}

private:
	virtual void init(const GraphAttributes& GA, double stretch, GraphCopySimple& spanner,
			EdgeArray<bool>& inSpanner) override {
		SpannerModule<TWeight>::init(GA, stretch, spanner, inSpanner);
		m_G.clear();
		m_G.init(GA.constGraph());
		m_cluster.init(m_G);
	}

	virtual typename SpannerModule<TWeight>::ReturnType execute() override {
		phase1();
		phase2();
		return SpannerModule<TWeight>::ReturnType::Feasible;
	}

	/**
	 * Adds \p e from \e m_G to the spanner and sets inSpanner
	 */
	inline void addToSpanner(edge e) {
		m_spanner->newEdge(m_G.original(e));
		(*m_inSpanner)[m_G.original(e)] = true;
	}

	/**
	 * \returns the weights of an edge \p e from m_G
	 */
	inline TWeight weight(edge e) { return getWeight(*m_GA, m_G.original(e)); }

	/**
	 * Phase 1: Forming clusters
	 */
	void phase1() {
		// init
		NodeSet<true> clusterCenters(m_G);
		for (node v : m_G.nodes) {
			m_cluster[v] = v; // At the beginning each node is it's own cluster
			clusterCenters.insert(v);
		}

		// A and B will be explained below. They are created here, so they can be reused.
		NodeArray<edge> A(m_G, nullptr);
		NodeArray<bool> B(m_G, false);

		// stretch = 2k-1
		// the stretch must be integer and uneven.
		// => k=(stretch+1)/2
		int k = (static_cast<int>(m_stretch) + 1) / 2;
		const double probability = pow(m_G.numberOfNodes(), -1.0 / k);
		for (int iteration = 1; iteration <= k - 1; iteration++) {
			assertTimeLeft();

			NodeArray<node> oldCluster = m_cluster; // oldCluster is the cluster of iteration i-1

			// 1: Sample new cluster centers
			NodeSet<true> sampledClusterCenters(m_G);
			for (node oldCenter : clusterCenters.nodes()) {
				if (randomDouble(0.0, 1.0) <= probability) {
					sampledClusterCenters.insert(oldCenter);
				}
			}

			// 2 & 3: Nearest neighbors and growing clusters
			for (node v : m_G.nodes) {
				if (sampledClusterCenters.isMember(oldCluster[v])) {
					continue;
				}

				// v is not a member of a sampled cluster, so we have to search for a nearest
				// neighbor
				// - v will be added to the cluster of the NN
				// - A[x] saves the least weight edge from v to the cluster member of the cluster x
				//   (or nullptr if there is no such edge.)

				// both "min" values are with respect to sampled clusters, not for all adjacent
				// clusters
				TWeight minWeight = std::numeric_limits<TWeight>::max();
				edge minEdge = nullptr;

				for (adjEntry a : v->adjEntries) {
					node center = oldCluster[a->twinNode()];
					edge e = a->theEdge();

					// maybe add v to a new cluster
					if (sampledClusterCenters.isMember(center)) {
						// twin is an element of a sampled cluster, so it is a (possible nearest)
						// neighbor
						if (m_eps.less(weight(e), minWeight)) {
							m_cluster[v] = center; // this line will automatically add v to the new
									// sampled cluster since the last assignment is
									// done on the min edge.
							minWeight = weight(e);
							minEdge = e;
						}
					}

					// fill A[x]
					if (A[center] == nullptr || m_eps.less(weight(e), weight(A[center]))) {
						A[center] = e;
					}
				}

				// B[<clustercenter>] stores whether all edges to a cluster should be deleted.

				// add sampled clusters (or all clusters) to v
				for (node center : m_G.nodes) {
					// center is only a possible center. It is a center, if A[center] != nullptr.
					// Some possibilities:
					// - minEdge == nullptr (case (a) in paper): v is not adjacent to any sampled cluster.
					//   So, add every least weight edge to the spanner.
					// - minEdge == A[center] (case (b), first part, in paper): v is adjacent to a sampled
					//   cluster and we have the min edge here: Add it!
					// - not minEdge, but less weight (case (b), second part, in paper): v is adjacent to a
					//   sampled cluster. The edge to the cluster has strictly less weight than the min edge.
					if (A[center] != nullptr
							&& (minEdge == nullptr || minEdge == A[center]
									|| m_eps.less(weight(A[center]), minWeight))) {
						addToSpanner(A[center]);
						B[center] = true;
					}
					A[center] = nullptr;
				}

#ifdef OGDF_DEBUG
				for (node _v : m_G.nodes) {
					OGDF_ASSERT(A[_v] == nullptr);
				}
#endif

				// Finally, delete all edges from v indicated by B
				ListPure<adjEntry> adjs; // copy, so we can iterate and delete.
				v->allAdjEntries(adjs);
				for (adjEntry a : adjs) {
					node center = oldCluster[a->twinNode()];
					if (B[center]) {
						m_G.delEdge(a->theEdge());
					}
					B[center] = false;
				}

#ifdef OGDF_DEBUG
				for (node _v : m_G.nodes) {
					OGDF_ASSERT(!B[_v]);
				}
#endif
			}

			// 4: Removing intra-cluster edges
			for (edge _e : m_GA->constGraph().edges) {
				edge e = m_G.copy(_e);
				if (e == nullptr) {
					continue;
				}

				if (m_cluster[e->source()] == m_cluster[e->target()]) {
					// inter-cluster connection.
					m_G.delEdge(e);
				}
			}

			clusterCenters = std::move(sampledClusterCenters);
		}
	}

	/**
	 * Phase 2: Vertex-Cluster Joining
	 */
	void phase2() {
		NodeArray<edge> A(m_G, nullptr);
		for (node v : m_G.nodes) {
			assertTimeLeft();

			for (adjEntry a : v->adjEntries) {
				node center = m_cluster[a->twinNode()];
				edge e = a->theEdge();

				// fill A[x]
				if (A[center] == nullptr || weight(e) < weight(A[center])) {
					A[center] = e;
				}
			}

			for (node center : m_G.nodes) {
				if (A[center] != nullptr) {
					addToSpanner(A[center]);
					A[center] = nullptr;
				}
			}

			// Delete all edges from v in the original graph.
			adjEntry a;
			while ((a = v->firstAdj()) != nullptr) {
				m_G.delEdge(a->theEdge());
			}

#ifdef OGDF_DEBUG
			for (node _v : m_G.nodes) {
				OGDF_ASSERT(A[_v] == nullptr);
			}
#endif
		}
	}

	using SpannerModule<TWeight>::getWeight;
	using SpannerModule<TWeight>::assertTimeLeft;
	using SpannerModule<TWeight>::m_GA;
	using SpannerModule<TWeight>::m_stretch;
	using SpannerModule<TWeight>::m_spanner;
	using SpannerModule<TWeight>::m_inSpanner;
};

/**
 * Use the ogdf::SpannerIteratedWrapper to execute the ogdf::SpannerBaswanaSen algorithm
 * up to 1000 times.
 *
 * @ingroup ga-spanner
 */
template<typename TWeight>
class SpannerBaswanaSenIterated : public SpannerIteratedWrapper<TWeight> {
public:
	SpannerBaswanaSenIterated()
		: SpannerIteratedWrapper<TWeight>(new SpannerBaswanaSen<TWeight>(), 1000) {};
};

}
