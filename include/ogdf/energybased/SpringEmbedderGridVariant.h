/** \file
 * \brief Declaration of Spring-Embedder Grid Variant algorithm.
 *
 * \author Carsten Gutwenger
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

#pragma once

#include <ogdf/module/LayoutModule.h>
#include <ogdf/basic/Array2D.h>
#include <ogdf/energybased/SpringForceModel.h>
#include <ogdf/basic/GraphCopyAttributes.h>

namespace ogdf {

//! The spring-embedder layout algorithm with force approximation using hte grid variant approach.
/**
 * @ingroup gd-energy
 *
 * The implementation used in SpringEmbedderGridVariant is based on
 * the following publication:
 *
 * Thomas M. J. Fruchterman, Edward M. Reingold: <i>%Graph Drawing by Force-directed
 * Placement</i>. Software - Practice and Experience 21(11), pp. 1129-1164, 1991.
 *
 * <H3>Optional parameters</H3>
 * Fruchterman/Reingold layout provides the following optional parameters.
 *
 * <table>
 *   <tr>
 *     <th><i>Option</i><th><i>Type</i><th><i>Default</i><th><i>Description</i>
 *   </tr><tr>
 *     <td><i>iterations</i><td>int<td>400
 *     <td>The number of iterations performed in the optimization.
 *   </tr><tr>
 *     <td><i>noise</i><td>bool<td>true
 *     <td>If set to true, (small) random perturbations are performed.
 *   </tr><tr>
 *     <td><i>minDistCC</i><td>double<td>20.0
 *     <td>The minimum distance between connected components.
 *   </tr><tr>
 *     <td><i>pageRatio</i><td>double<td>1.0
 *     <td>The page ratio.
 *   </tr><tr>
 *     <td><i>scaling</i><td> #Scaling <td> Scaling::scaleFunction
 *     <td>The scaling method for scaling the inital layout.
 *   </tr><tr>
 *     <td><i>scaleFunctionFactor</i><td>double<td>8.0
 *     <td>The scale function factor (used if scaling = scaleFunction).
 *   </tr><tr>
 *     <td><i>userBoundingBox</i><td>rectangle<td>(0.0,100.0,0.0,100.0)
 *     <td>The user bounding box for scaling (used if scaling = scUserBoundingBox).
 *   </tr>
 * </table>
 */
class OGDF_EXPORT SpringEmbedderGridVariant : public LayoutModule
{
public:
	//! The scaling method used by the algorithm.
	enum class Scaling {
		input,             //!< bounding box of input is used.
		userBoundingBox,   //!< bounding box set by userBoundingBox() is used.
		scaleFunction,     //!< automatic scaling is used with parameter set by scaleFunctionFactor() (larger factor, larger b-box).
		useIdealEdgeLength //!< use the given ideal edge length to scale the layout suitably.
	};


	//! Creates an instance of Fruchterman/Reingold layout.
	SpringEmbedderGridVariant();

	// destructor
	~SpringEmbedderGridVariant() { }


	//! Calls the layout algorithm for graph attributes \a GA.
	virtual void call(GraphAttributes &GA) override;
	//void call_mt(GraphAttributes &GA);


	//! Returns the currently used force model.
	SpringForceModel forceModel() const {
		return m_forceModel;
	}

	//! Sets the used force model to \a fm.
	void forceModel(SpringForceModel fm) {
		m_forceModel = fm;
	}

	//! Returns the currently used force model for the improvement step.
	SpringForceModel forceModelImprove() const {
		return m_forceModelImprove;
	}

	//! Sets the used force model for the improvement step to \a fm.
	void forceModelImprove(SpringForceModel fm) {
		m_forceModelImprove = fm;
	}

	//! Returns the currently used <i>average convergence factor</i>.
	/**
	 * This factor is used for detecting convergence of the energy system.
	 * With respect to the average displacement of a node in a single step, we assume
	 * to have convergence if it is at most \a avgConvergenceFactor * \a idealEdgeLength.
	 */
	double avgConvergenceFactor() const {
		return m_avgConvergenceFactor;
	}

	//! Sets the <i>average convergence factor</i> to \a f.
	void avgConvergenceFactor(double f) {
		if(f >= 0)
			m_avgConvergenceFactor = f;
	}

	//! Returns the currently used <i>maximum</i> convergence factor.
	/**
	 * This factor is used for detecting convergence of the energy system.
	 * With respect to the maximum displacement of a node in a single step, we assume
	 * to have convergence if it is at most \a maxConvergenceFactor * \a idealEdgeLength.
	 */
	double maxConvergenceFactor() const {
		return m_maxConvergenceFactor;
	}

	//! Sets the <i>maximum</i> convergence factor to \a f.
	void maxConvergenceFactor(double f) {
		if(f >= 0)
			m_maxConvergenceFactor = f;
	}

	//! Returns the current setting of iterations.
	/**
	 * This setting limits the number of optimization rounds. If convergence (with respect to node displacement)
	 * is detected, the optimization process is immediately finished.
	 */
	int iterations() const {
		return m_iterations;
	}

	//! Sets the number of iterations to \a i.
	void iterations(int i) {
		if (i >= 0)
			m_iterations = i;
	}

	//! Returns the current setting of iterations for the improvement phase.
	int iterationsImprove() const {
		return m_iterationsImprove;
	}

	//! Sets the number of iterations for the improvement phase to \a i.
	void iterationsImprove(int i) {
		if (i >= 0)
			m_iterationsImprove = i;
	}

	//! Returns the current setting of ideal edge length.
	double idealEdgeLength() const {
		return m_idealEdgeLength;
	}

	double coolDownFactor() const {
		return m_coolDownFactor;
	}

	double forceLimitStep() const {
		return m_forceLimitStep;
	}

	//! Sets the ideal edge length to \a len.
	/**
	 * Edge lengths are measured between the centers of the two nodes, i.e.,
	 * node sizes are not taken into account.
	 */
	void idealEdgeLength(double len) {
		m_idealEdgeLength = len;
	}

	//! Returns the current setting of noise.
	bool noise() const {
		return m_noise;
	}

	//! Sets the parameter noise to \a on.
	void noise(bool on) {
		m_noise = on;
	}

	//! Returns the minimum distance between connected components.
	double minDistCC() const { return m_minDistCC; }

	//! Sets the minimum distance between connected components to \a x.
	void minDistCC(double x) { m_minDistCC = x; }

	//! Returns the page ratio.
	double pageRatio() { return m_pageRatio; }

	//! Sets the page ration to \a x.
	void pageRatio(double x) { m_pageRatio = x; }

	//! Returns the current scaling method.
	Scaling scaling() const {
		return m_scaling;
	}

	//! Sets the method for scaling the inital layout to \a sc.
	void scaling(Scaling sc) {
		m_scaling = sc;
	}

	//! Returns the current scale function factor.
	double scaleFunctionFactor() const {
		return m_scaleFactor;
	}

	//! Sets the scale function factor to \a f.
	void scaleFunctionFactor(double f) {
		m_scaleFactor = f;
	}

	//! Sets the user bounding box (used if scaling method is scUserBoundingBox).
	void userBoundingBox(double xmin, double ymin, double xmax, double ymax) {
		m_bbXmin = xmin;
		m_bbYmin = ymin;
		m_bbXmax = xmax;
		m_bbYmax = ymax;
	}

	//! Returns the maximal number of used threads.
	unsigned int maxThreads() const { return m_maxThreads; }

	//! Sets the maximal number of used threads to \a n.
	void maxThreads(unsigned int n) { m_maxThreads = n; }

private:
	struct NodeInfo
	{
		DPoint m_pos;

		int m_adjBegin;
		int m_adjStop;

		int m_gridX;
		int m_gridY;

		ListIterator<int> m_lit;
	};


	class ForceModelBase;
	class ForceModelFR;
	class ForceModelFRModAttr;
	class ForceModelFRModRep;
	class ForceModelEades;
	class ForceModelHachul;
	class ForceModelGronemann;

	class Master;
	class Worker;

	int    m_iterations;         //!< The number of iterations.
	int    m_iterationsImprove;  //!< The number of iterations for the improvement phase.
	double m_idealEdgeLength;    //!< The ideal edge length.
	double m_coolDownFactor;
	double m_forceLimitStep;

	double m_xleft;       //!< Bounding box (minimal x-coordinate).
	double m_xright;      //!< Bounding box (maximal x-coordinate).
	double m_ysmall;      //!< Bounding box (minimal y-coordinate).
	double m_ybig;        //!< Bounding box (maximal y-coordinate).

	SpringForceModel m_forceModel; //! The used force model.
	SpringForceModel m_forceModelImprove; //! The used force model for the improvement phase.
	bool m_noise;            //!< Perform random perturbations?

	Scaling m_scaling;    //!< The scaling method.
	double m_scaleFactor; //!< The factor used if scaling type is scScaleFunction.

	double m_bbXmin; //!< User bounding box (minimal x-coordinate).
	double m_bbYmin; //!< User bounding box (maximal x-coordinate).
	double m_bbXmax; //!< User bounding box (minimal y-coordinate).
	double m_bbYmax; //!< User bounding box (maximal y-coordinate).

	double m_minDistCC; //!< The minimal distance between connected components.
	double m_pageRatio; //!< The page ratio.

	double m_avgConvergenceFactor; //!< convergence if avg. displacement is at most this factor times ideal edge length
	double m_maxConvergenceFactor; //!< convergence if max. displacement is at most this factor times ideal edge length

	unsigned int m_maxThreads;	//!< The maximal number of used threads.
};


} // end namespace ogdf
