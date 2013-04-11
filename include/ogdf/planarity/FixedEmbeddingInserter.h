/*
 * $Revision: 3368 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-04 20:07:31 +0200 (Do, 04. Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration of class FixedEmbeddingInserter.
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

#ifdef _MSC_VER
#pragma once
#endif

#ifndef OGDF_FIXED_EMB_INSERTER_H
#define OGDF_FIXED_EMB_INSERTER_H


#include <ogdf/module/EdgeInsertionModule.h>
#include <ogdf/planarity/RemoveReinsertType.h>


namespace ogdf {


	//! Edge insertion module that inserts each edge optimally into a fixed embedding.
	class OGDF_EXPORT FixedEmbeddingInserter : public EdgeInsertionModule
	{
	public:
		//! Creates an instance of fixed embedding edge inserter with default settings.
		FixedEmbeddingInserter();

		//! Creates an instance of fixed embedding edge inserter with the same settings as \a inserter.
		FixedEmbeddingInserter(const FixedEmbeddingInserter &inserter);

		//! Destructor.
		~FixedEmbeddingInserter() { }

		//! Returns a new instance of the fixed embedding inserter with the same option settings.
		EdgeInsertionModule *clone() const;

		//! Assignment operator. Copies option settings only.
		FixedEmbeddingInserter &operator=(const FixedEmbeddingInserter &inserter);


		/**
		 *  @name Optional parameters
		 *  @{
		 */

		//! Sets the remove-reinsert postprocessing method.
		void removeReinsert(RemoveReinsertType rrOption) {
			m_rrOption = rrOption;
		}

		//! Returns the current setting of the remove-reinsert postprocessing method.
		RemoveReinsertType removeReinsert() const {
			return m_rrOption;
		}


		//! Sets the option <i>percentMostCrossed</i> to \a percent.
		/**
		 * This option determines the portion of most crossed edges used if the remove-reinsert
		 * method is set to #rrMostCrossed. This portion is number of edges * percentMostCrossed() / 100.
		 */
		void percentMostCrossed(double percent) {
			m_percentMostCrossed = percent;
		}

		//! Returns the current setting of option percentMostCrossed.
		double percentMostCrossed() const {
			return m_percentMostCrossed;
		}

		//! Sets the option <i>keepEmbedding</i> to \a keep.
		/**
		 * This option determines if the planar embedding of the planarized representation \a PG passed to the call-method
		 * is preserved, or if always a new embedding is computed. If <i>keepEmbedding</i> is set to true,
		 * \a PG must always be planarly embedded.
		 */
		void keepEmbedding(bool keep) {
			m_keepEmbedding = keep;
		}

		//! Returns the current setting of option <i>keepEmbedding</i>.
		bool keepEmbeding() const {
			return m_keepEmbedding;
		}

		/** @}
		 *  @name Further information
		 *  @{
		 */

		//! Returns the number of runs performed by the remove-reinsert method after the algorithm has been called.
		int runsPostprocessing() const {
			return m_runsPostprocessing;
		}

		//! @}

	private:
		//! Implements the algorithm call.
		virtual ReturnType doCall(
			PlanRepLight              &pr,
			const Array<edge>         &origEdges,
			const EdgeArray<int>      *costOrig,
			const EdgeArray<bool>     *pForbiddenOrig,
			const EdgeArray<__uint32> *pEdgeSubGraphs);

		RemoveReinsertType m_rrOption; //!< The remove-reinsert method.
		double m_percentMostCrossed;   //!< The portion of most crossed edges considered.
		bool m_keepEmbedding;

		int m_runsPostprocessing; //!< Runs of remove-reinsert method.
	};

} // end namespace ogdf

#endif
