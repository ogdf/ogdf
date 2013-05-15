/*
 * $Revision: 3418 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-18 11:06:44 +0200 (Do, 18. Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of the pivot MDS.
 *
 * \author Mark Ortmann, University of Konstanz
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.txt in the root directory of the OGDF installation for details.
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
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/

#include <ogdf/energybased/PivotMDS.h>


namespace ogdf {

const double PivotMDS::EPSILON = 1 - 1e-10;
const double PivotMDS::FACTOR = -0.5;


void PivotMDS::call(GraphAttributes& GA)
{
	if (DIMENSION_COUNT > 2) {
		OGDF_ASSERT(GA.attributes() & GraphAttributes::threeD);
	}
	if (!isConnected(GA.constGraph())) {
		OGDF_THROW_PARAM(PreconditionViolatedException,pvcConnected);
		return;
	}
	if (m_hasEdgeCostsAttribute
			&& !(GA.attributes() & GraphAttributes::edgeDoubleWeight)) {
				OGDF_THROW(PreconditionViolatedException);
		return;
	}
	pivotMDSLayout(GA);
}


void PivotMDS::centerPivotmatrix(Array<Array<double> >& pivotMatrix)
{
	int numberOfPivots = pivotMatrix.size();
	// this is ensured since the graph size is at least 2!
	int nodeCount = pivotMatrix[0].size();

	double normalizationFactor = 0;
	double rowColNormalizer;
	Array<double> colNormalization(numberOfPivots);

	for (int i = 0; i < numberOfPivots; i++) {
		rowColNormalizer = 0;
		for (int j = 0; j < nodeCount; j++) {
			rowColNormalizer += pivotMatrix[i][j] * pivotMatrix[i][j];
		}
		normalizationFactor += rowColNormalizer;
		colNormalization[i] = rowColNormalizer / nodeCount;
	}
	normalizationFactor = normalizationFactor / (nodeCount * numberOfPivots);
	for (int i = 0; i < nodeCount; i++) {
		rowColNormalizer = 0;
		for (int j = 0; j < numberOfPivots; j++) {
			double square = pivotMatrix[j][i] * pivotMatrix[j][i];
			pivotMatrix[j][i] = square + normalizationFactor
					- colNormalization[j];
			rowColNormalizer += square;
		}
		rowColNormalizer /= numberOfPivots;
		for (int j = 0; j < numberOfPivots; j++) {
			pivotMatrix[j][i] = FACTOR * (pivotMatrix[j][i] - rowColNormalizer);
		}
	}
}


void PivotMDS::pivotMDSLayout(GraphAttributes& GA)
{
	const Graph& G = GA.constGraph();
	if (G.numberOfNodes() <= 1) {
		// make it exception save
		node v;
		forall_nodes(v,G)
		{
			GA.x(v) = 0.0;
			GA.y(v) = 0.0;
			if (DIMENSION_COUNT > 2)
				GA.z(v) = 0.0;
		}
		return;
	}
	// check whether the graph is a path or not
	const node head = getRootedPath(G);
	if (head != 0) {
		doPathLayout(GA, head);
	} else {
		Array<Array<double> > pivDistMatrix;
		// compute the pivot matrix
		getPivotDistanceMatrix(GA, pivDistMatrix);
		// center the pivot matrix
		centerPivotmatrix(pivDistMatrix);
		// init the coordinate matrix
		Array<Array<double> > coord(DIMENSION_COUNT);
		for (int i = 0; i < coord.size(); i++) {
			coord[i].init(G.numberOfNodes());
		}
		// init the eigen values array
		Array<double> eVals(DIMENSION_COUNT);
		singularValueDecomposition(pivDistMatrix, coord, eVals);
		// compute the correct aspect ratio
		for (int i = 0; i < coord.size(); i++) {
			eVals[i] = sqrt(eVals[i]);
			for (int j = 0; j < G.numberOfNodes(); j++) {
				coord[i][j] *= eVals[i];
			}
		}
		// set the new positions to the graph
		int i = 0;
		node v;
		forall_nodes(v,G)
		{
			GA.x(v) = coord[0][i];
			GA.y(v) = coord[1][i];
			if (DIMENSION_COUNT > 2){
				GA.z(v) = coord[2][i];//cout << coord[2][i] << "\n";
			}
			++i;
		}
	}
}


void PivotMDS::doPathLayout(GraphAttributes& GA, const node& v)
{
	double xPos = 0;
	node prev = v;
	node cur = v;
	edge e;
	// since the given node is the beginning of the path just
	// use bfs and increment the x coordinate by the average
	// edge costs.
	do {
		GA.x(cur) = xPos;
		GA.y(cur) = 0;
		node adj;
		forall_adj_edges(e,cur) {
			adj = e->opposite(cur);
			if (!(adj == prev) || adj == cur) {
				prev = cur;
				cur = adj;
				if(m_hasEdgeCostsAttribute) {
					xPos+=GA.doubleWeight(e);
				} else {
					xPos += m_edgeCosts;
				}
				break;
			}
			prev = cur;
		}
	} while (prev != cur);
}


void PivotMDS::eigenValueDecomposition(
	Array<Array<double> >& K,
	Array<Array<double> >& eVecs,
	Array<double>& eValues)
{
	randomize(eVecs);
	const int p = K.size();
	double r = 0;
	for (int i = 0; i < DIMENSION_COUNT; i++) {
		eValues[i] = normalize(eVecs[i]);
	}
	while (r < EPSILON) {
		if (isnan(r) || isinf(r)) {
			// Throw arithmetic exception (Shouldn't occur
			// for DIMEMSION_COUNT = 2
			OGDF_THROW(AlgorithmFailureException);
			return;
		}
		// remember prev values
		Array<Array<double> > tmpOld(DIMENSION_COUNT);
		for (int i = 0; i < DIMENSION_COUNT; i++) {
			tmpOld[i].init(p);
			for (int j = 0; j < p; j++) {
				tmpOld[i][j] = eVecs[i][j];
				eVecs[i][j] = 0;
			}
		}
		// multiply matrices
		for (int i = 0; i < DIMENSION_COUNT; i++) {
			for (int j = 0; j < p; j++) {
				for (int k = 0; k < p; k++) {
					eVecs[i][k] += K[j][k] * tmpOld[i][j];
				}
			}
		}
		// orthogonalize
		for (int i = 0; i < DIMENSION_COUNT; i++) {
			for (int j = 0; j < i; j++) {
				double fac = prod(eVecs[j], eVecs[i])
						/ prod(eVecs[j], eVecs[j]);
				for (int k = 0; k < p; k++) {
					eVecs[i][k] -= fac * eVecs[j][k];
				}
			}
		}
		// normalize
		for (int i = 0; i < DIMENSION_COUNT; i++) {
			eValues[i] = normalize(eVecs[i]);
		}
		r = 1;
		for (int i = 0; i < DIMENSION_COUNT; i++) {
			// get absolute value (abs only defined for int)
			double tmp = prod(eVecs[i], tmpOld[i]);
			if (tmp < 0) {
				tmp *= -1;
			}
			r = min(r, tmp);
		}
	}
}


void PivotMDS::getPivotDistanceMatrix(
	const GraphAttributes& GA,
	Array<Array<double> >& pivDistMatrix)
{
	const Graph& G = GA.constGraph();
	// lower the number of pivots if necessary
	int numberOfPivots = min(G.numberOfNodes(), m_numberOfPivots);
	// number of pivots times n matrix used to store the graph distances
	pivDistMatrix.init(numberOfPivots);
	for (int i = 0; i < numberOfPivots; i++) {
		pivDistMatrix[i].init(G.numberOfNodes());
	}
	// edges costs array
	EdgeArray<double> edgeCosts;
	bool hasEdgeCosts = false;
	// already checked whether this attribute exists or not (see call method)
	if (m_hasEdgeCostsAttribute) {
		edgeCosts.init(G);
		edge e;
		forall_edges(e, G)
		{
			edgeCosts[e] = GA.doubleWeight(e);
		}
		hasEdgeCosts = true;
	}
	// used for min-max strategy
	NodeArray<double> minDistances(G, std::numeric_limits<double>::infinity());
	NodeArray<double> shortestPathSingleSource(G);
	// the current pivot node
	node pivNode = G.firstNode();
	for (int i = 0; i < numberOfPivots; i++) {
		// get the shortest path from the currently processed pivot node to
		// all other nodes in the graph
		shortestPathSingleSource.fill(std::numeric_limits<double>::infinity());
		if (hasEdgeCosts) {
			dijkstra_SPSS(pivNode, G, shortestPathSingleSource, edgeCosts);
		} else {
			bfs_SPSS(pivNode, G, shortestPathSingleSource, m_edgeCosts);
		}
		copySPSS(pivDistMatrix[i], shortestPathSingleSource);
		// update the pivot and the minDistances array ... to ensure the
		// correctness set minDistance of the pivot node to zero
		minDistances[pivNode] = 0;
		node v;
		forall_nodes(v,G)
		{
			minDistances[v] = min(minDistances[v], shortestPathSingleSource[v]);
			if (minDistances[v] > minDistances[pivNode]) {
				pivNode = v;
			}
		}
	}
}


void PivotMDS::copySPSS(Array<double>& copyTo, NodeArray<double>& copyFrom)
{
	for (int i = 0; i < copyTo.size(); i++) {
		copyTo[i] = copyFrom[i];
	}
}


node PivotMDS::getRootedPath(const Graph& G)
{
	node head = 0;
	int degree;
	edge e;
	node v;
	node adj;
	NodeArray<bool> visited(G, false);
	SListPure<node> neighbors;
	// in every path there are two nodes with degree 1 and
	// each node has at most degree 2
	forall_nodes(v,G)
	{
		degree = 0;
		visited[v] = true;
		neighbors.pushBack(v);
		forall_adj_edges(e,v) {
			adj = e->opposite(v);
			if (!visited[adj])
			{
				neighbors.pushBack(adj);
				visited[adj]=true;
				++degree;
			}
		}
		if (degree > 2) {
			neighbors.clear();
			return 0;
		}
		if (degree == 1) {
			head = v;
		}
		for(SListConstIterator<node> it = neighbors.begin(); it.valid(); ++it) {
			visited[*it] = false;
		}
		neighbors.clear();
	}
	return head;
}


double PivotMDS::normalize(Array<double>& x)
{
	double norm = sqrt(prod(x, x));
	if (norm != 0) {
		for (int i = 0; i < x.size(); i++) {
			x[i] /= norm;
		}
	}
	return norm;
}


double PivotMDS::prod(const Array<double>& x, const Array<double>& y)
{
	double result = 0;
	for (int i = 0; i < x.size(); i++) {
		result += x[i] * y[i];
	}
	return result;
}


void PivotMDS::randomize(Array<Array<double> >& matrix)
{
	srand(SEED);
	for (int i = 0; i < matrix.size(); i++) {
		for (int j = 0; j < matrix[i].size(); j++) {
			matrix[i][j] = ((double) rand()) / RAND_MAX;
		}
	}
}


void PivotMDS::selfProduct(const Array<Array<double> >& d, Array<Array<double> >& result)
{
	double sum;
	for (int i = 0; i < d.size(); i++) {
		for (int j = 0; j <= i; j++) {
			sum = 0;
			for (int k = 0; k < d[0].size(); k++) {
				sum += d[i][k] * d[j][k];
			}
			result[i][j] = sum;
			result[j][i] = sum;
		}
	}
}


void PivotMDS::singularValueDecomposition(
	Array<Array<double> >& pivDistMatrix,
	Array<Array<double> >& eVecs,
	Array<double>& eVals)
{
	const int l = pivDistMatrix.size();
	const int n = pivDistMatrix[0].size();
	Array<Array<double> > K(l);
	for (int i = 0; i < l; i++) {
		K[i].init(l);
	}
	// calc C^TC
	selfProduct(pivDistMatrix, K);

	Array<Array<double> > tmp(DIMENSION_COUNT);
	for (int i = 0; i < DIMENSION_COUNT; i++) {
		tmp[i].init(l);
	}

	eigenValueDecomposition(K, tmp, eVals);

	// C^Tx
	for (int i = 0; i < DIMENSION_COUNT; i++) {
		eVals[i] = sqrt(eVals[i]);
		for (int j = 0; j < n; j++) { // node j
			eVecs[i][j] = 0;
			for (int k = 0; k < l; k++) { // pivot k
				eVecs[i][j] += pivDistMatrix[k][j] * tmp[i][k];
			}
		}
	}
	for (int i = 0; i < DIMENSION_COUNT; i++) {
		normalize(eVecs[i]);
	}
}

} /* namespace ogdf */
