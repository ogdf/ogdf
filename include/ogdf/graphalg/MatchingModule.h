/** \file
 * \brief Declaration of interface for matching algorithms
 *
 * \author Joshua Sangmeister
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

#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/Logger.h>
#include <ogdf/basic/Module.h>
#include <ogdf/basic/NodeArray.h>
#include <ogdf/graphalg/matching_blossom/utils.h>

#include <tuple>
#include <unordered_set>

namespace ogdf {

using namespace matching_blossom;

template<typename TWeight>
class MatchingModule : public Module, public Logger {
public:
	virtual ~MatchingModule() { }

	/**
	 * @brief Computes a minimum weight perfect matching in \p G.
	 *
	 * @param G the graph to compute the matching in
	 * @param weights the weights of the edges of G
	 * @param matching output parameter to store the matching in
	 * @returns whether or not a perfect matching was found
	 */
	bool minimumWeightPerfectMatching(const Graph& G, const EdgeArray<TWeight>& weights,
			std::unordered_set<edge>& matching) {
		return doCall(G, weights, matching);
	}

	/**
	 * @brief Computes a minimum weight perfect matching in \p GA.
	 *
	 * @param GA The graph to compute the matching in, together with its edge weights
	 * @param matching output parameter to store the matching in
	 * @returns whether or not a perfect matching was found
	 */
	bool minimumWeightPerfectMatching(const GraphAttributes& GA, std::unordered_set<edge>& matching) {
		return doCall(GA, matching);
	}

	/**
	 * @brief Computes a minimum weight perfect matching in \p G.
	 *
	 * @param G the graph to compute the matching in
	 * @param weights the weights of the edges of G
	 * @returns a tuple consisting of a boolean indicating whether or not a perfect matching was
	 * 			found and the weight of the matching, if one was found
	 */
	std::tuple<bool, TWeight> minimumWeightPerfectMatching(const Graph& G,
			const EdgeArray<TWeight>& weights) {
		std::unordered_set<edge> matching;
		bool result = minimumWeightPerfectMatching(G, weights, matching);
		TWeight weight = matchingWeight(matching, weights);
		return std::make_tuple(result, weight);
	}

	/**
	 * @brief Computes a minimum weight perfect matching in \p G.
	 *
	 * @param GA The graph to compute the matching in, together with its edge weights
	 * @returns a tuple consisting of a boolean indicating whether or not a perfect matching was
	 * 			found and the weight of the matching, if one was found
	 */
	std::tuple<bool, TWeight> minimumWeightPerfectMatching(const GraphAttributes& GA) {
		std::unordered_set<edge> matching;
		bool result = minimumWeightPerfectMatching(GA, matching);
		TWeight weight = matchingWeight(matching, GA);
		return std::make_tuple(result, weight);
	}

	/**
	 * @brief Computes a maximum weight perfect matching in \p G.
	 *
	 * @param G the graph to compute the matching in
	 * @param weights the weights of the edges of G
	 * @param matching output parameter to store the matching in
	 * @returns whether or not a perfect matching was found
	 */
	bool maximumWeightPerfectMatching(const Graph& G, const EdgeArray<TWeight>& weights,
			std::unordered_set<edge>& matching) {
		EdgeArray<TWeight> invertedWeights(G);
		copyWeights(G, weights, invertedWeights, true);
		return doCall(G, invertedWeights, matching);
	}

	/**
	 * @brief Computes a maximum weight perfect matching in \p GA.
	 *
	 * @param GA The graph to compute the matching in, together with its edge weights
	 * @param matching output parameter to store the matching in
	 * @returns whether or not a perfect matching was found
	 */
	bool maximumWeightPerfectMatching(const GraphAttributes& GA, std::unordered_set<edge>& matching) {
		const Graph& G = GA.constGraph();
		EdgeArray<TWeight> weights(G);
		copyWeights(G, GA, weights, true);
		return doCall(G, weights, matching);
	}

	/**
	 * @brief Computes a maximum weight perfect matching in \p G.
	 *
	 * @param G the graph to compute the matching in
	 * @param weights the weights of the edges of G
	 * @returns a tuple consisting of a boolean indicating whether or not a perfect matching was
	 * 			found and the weight of the matching, if one was found
	 */
	std::tuple<bool, TWeight> maximumWeightPerfectMatching(const Graph& G,
			const EdgeArray<TWeight>& weights) {
		std::unordered_set<edge> matching;
		bool result = maximumWeightPerfectMatching(G, weights, matching);
		TWeight weight = matchingWeight(matching, weights);
		return std::make_tuple(result, weight);
	}

	/**
	 * @brief Computes a maximum weight perfect matching in \p G.
	 *
	 * @param GA The graph to compute the matching in, together with its edge weights
	 * @returns a tuple consisting of a boolean indicating whether or not a perfect matching was
	 * 			found and the weight of the matching, if one was found
	 */
	std::tuple<bool, TWeight> maximumWeightPerfectMatching(const GraphAttributes& GA) {
		std::unordered_set<edge> matching;
		bool result = maximumWeightPerfectMatching(GA, matching);
		TWeight weight = matchingWeight(matching, GA);
		return std::make_tuple(result, weight);
	}

	/**
	 * @brief Computes a maximum weight matching in \p G.
	 *
	 * @param G the graph to compute the matching in
	 * @param weights the weights of the edges of G
	 * @param matching output parameter to store the matching in
	 */
	void maximumWeightMatching(const Graph& G, const EdgeArray<TWeight>& weights,
			std::unordered_set<edge>& matching) {
		doMaximumWeightMatching(G, weights, matching);
	}

	/**
	 * @brief Computes a maximum weight matching in \p GA.
	 *
	 * @param GA The graph to compute the matching in, together with its edge weights
	 * @param matching output parameter to store the matching in
	 */
	void maximumWeightMatching(const GraphAttributes& GA, std::unordered_set<edge>& matching) {
		doMaximumWeightMatching(GA.constGraph(), GA, matching);
	}

	/**
	 * @brief Computes a maximum weight matching in \p G.
	 *
	 * @param G the graph to compute the matching in
	 * @param weights the weights of the edges of G
	 * @returns the weight of the maximum matching
	 */
	TWeight maximumWeightMatching(const Graph& G, const EdgeArray<TWeight>& weights) {
		std::unordered_set<edge> matching;
		maximumWeightMatching(G, weights, matching);
		return matchingWeight(matching, weights);
	}

	/**
	 * @brief Computes a maximum weight matching in \p G.
	 *
	 * @param GA The graph to compute the matching in, together with its edge weights
	 * @returns the weight of the maximum matching
	 */
	TWeight maximumWeightMatching(const GraphAttributes& GA) {
		std::unordered_set<edge> matching;
		maximumWeightMatching(GA, matching);
		return matchingWeight(matching, GA);
	}

	//! Calculate the weight of the \p matching with the weights given in \p weights.
	TWeight matchingWeight(const std::unordered_set<edge>& matching,
			const EdgeArray<TWeight>& weights) {
		return getMatchingWeight(matching, weights);
	}

	//! Calculate the weight of the \p matching with the weights given in \p GA.
	TWeight matchingWeight(const std::unordered_set<edge>& matching, const GraphAttributes& GA) {
		return getMatchingWeight(matching, GA);
	}

protected:
	//! Main call function for the algorithm with an EdgeArray as weight container.
	virtual bool doCall(const Graph& G, const EdgeArray<TWeight>& weights,
			std::unordered_set<edge>& matching) = 0;

	//! Main call function for the algorithm with GraphAttrubtes as weight container.
	virtual bool doCall(const GraphAttributes& GA, std::unordered_set<edge>& matching) = 0;

private:
	template<class WeightContainer>
	void copyWeights(const Graph& G, const WeightContainer& weights, EdgeArray<TWeight>& copy,
			bool invert = false) {
		for (edge e : G.edges) {
			copy[e] = (invert ? -1 : 1) * getWeight<TWeight>(e, weights);
		}
	}

	template<class WeightContainer>
	TWeight getMatchingWeight(const std::unordered_set<edge>& matching,
			const WeightContainer& weights) {
		TWeight value = 0;
		for (edge e : matching) {
			value += getWeight<TWeight>(e, weights);
		}
		return value;
	}

	template<class WeightContainer>
	void doMaximumWeightMatching(const Graph& G, const WeightContainer& weights,
			std::unordered_set<edge>& matching) {
		// Duplicate the graph and connect all nodes with their copies with 0-weight edges
		GraphCopySimple graphCopy(G);
		EdgeArray<TWeight> weightsCopy(graphCopy, 0);
		NodeArray<node> nodeMap(graphCopy);
		for (node origNode : G.nodes) {
			node v = graphCopy.copy(origNode);
			node dup = graphCopy.newNode();
			nodeMap[v] = dup;
			graphCopy.newEdge(v, dup);
		}
		for (edge origEdge : G.edges) {
			edge e = graphCopy.copy(origEdge);
			edge newEdge = graphCopy.newEdge(nodeMap[e->source()], nodeMap[e->target()]);
			weightsCopy[newEdge] = weightsCopy[e] = -getWeight<TWeight>(origEdge, weights);
		}

		// Calculate the MinimumWeightPerfectMatching on the new graph
		std::unordered_set<edge> matchingCopy;
#ifdef OGDF_DEBUG
		bool result =
#endif
				doCall(graphCopy, weightsCopy, matchingCopy);
		OGDF_ASSERT(result);

		// Convert the matching to the original graph
		for (edge e : matchingCopy) {
			if (edge origEdge = graphCopy.original(e)) {
				matching.insert(origEdge);
			}
		}
	};
};

}
