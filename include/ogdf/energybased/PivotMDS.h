/*
 * $Revision: 3415 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-17 15:40:21 +0200 (Mi, 17. Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration of the pivot MDS. By setting the number of pivots to
 * infinity this algorithm behaves just like classical MDS.
 * See Brandes and Pich: Eigensolver methods for progressive multidi-
 * mensional scaling of large data.
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

#ifndef OGDF_PIVOT_MDS_H
#define OGDF_PIVOT_MDS_H


#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/graphalg/ShortestPathAlgorithms.h>
#include <ogdf/module/LayoutModule.h>


namespace ogdf {

#ifndef isnan
template<typename T>
inline bool isnan(T value)
{
	return value != value;
}
#endif


#ifndef isinf
// requires #include <limits>
template<typename T>
inline bool isinf(T value)
{
return std::numeric_limits<T>::has_infinity &&
value == std::numeric_limits<T>::infinity();
}
#endif


class OGDF_EXPORT PivotMDS : public LayoutModule {
public:
	PivotMDS() : m_numberOfPivots(250), m_edgeCosts(100), m_hasEdgeCostsAttribute(false) { }

	virtual ~PivotMDS() { }

	//! Sets the number of pivots. If the new value is smaller or equal 0
	//! the default value (250) is used.
	void setNumberOfPivots(int numberOfPivots) {
		m_numberOfPivots = (numberOfPivots < DIMENSION_COUNT) ? DIMENSION_COUNT : numberOfPivots;
	}

	//! Sets the desired distance between adjacent nodes. If the new value is smaller or equal
	//! 0 the default value (100) is used.
	void setEdgeCosts(double edgeCosts){
		m_edgeCosts = edgeCosts;
	}

	//! Calls the layout algorithm for graph attributes \a GA.
	void call(GraphAttributes& GA);

	void useEdgeCostsAttribute(bool useEdgeCostsAttribute) {
		m_hasEdgeCostsAttribute = useEdgeCostsAttribute;
	}

private:

	//! The dimension count determines the number of evecs that
	//! will be computed. Nevertheless PivotMDS only takes the first two
	//! with the highest eigenwert into account.
	const static int DIMENSION_COUNT = 3;

	//! Convergence factor used for power iteration.
	const static double EPSILON;

	//! Factor used to center the pivot matrix.
	const static double FACTOR;

	//! Seed of the random number generator.
	const static unsigned int SEED = 0;

	//! The number of pivots.
	int m_numberOfPivots;

	//! The costs to traverse an edge.
	double m_edgeCosts;

	//! Tells whether the pivot mds is based on uniform edge costs or a
	//! edge costs attribute
	bool m_hasEdgeCostsAttribute;

	//! Centers the pivot matrix.
	void centerPivotmatrix(Array<Array<double> >& pivotMatrix);

	//! Computes the pivot mds layout of the given connected graph of \a GA.
	void pivotMDSLayout(GraphAttributes& GA);

	void copySPSS(Array<double>& copyTo, NodeArray<double>& copyFrom);

	//! Computes the layout of a path.
	void doPathLayout(GraphAttributes& GA, const node& v);

	//! Computes the eigen value decomposition based on power iteration.
	void eigenValueDecomposition(
		Array<Array<double> >& K,
		Array<Array<double> >& eVecs,
		Array<double>& eValues);

	//! Computes the pivot distance matrix based on the maxmin strategy
	void getPivotDistanceMatrix(const GraphAttributes& GA, Array<Array<double> >& pivDistMatrix);

	//! Checks whether the given graph is a path or not.
	node getRootedPath(const Graph& G);

	//! Normalizes the vector \a x.
	double normalize(Array<double>& x);

	//! Computes the product of two vectors \a x and \a y.
	double prod(const Array<double>& x, const Array<double>& y);

	//! Fills the given \a matrix with random doubles d 0 <= d <= 1.
	void randomize(Array<Array<double> >& matrix);

	//! Computes the self product of \a d.
	void selfProduct(const Array<Array<double> >&d, Array<Array<double> >& result);

	//! Computes the singular value decomposition of matrix \a K.
	void singularValueDecomposition(
		Array<Array<double> >& K,
		Array<Array<double> >& eVecs,
		Array<double>& eVals);
};


} // end namespace ogdf

#endif
