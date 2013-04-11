/*
 * $Revision: 3368 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-04 20:07:31 +0200 (Do, 04. Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Definition of RemoveReinsertType (used for postprocessing
 *        in edge insertion algorithms).
 *
 * \author Karsten Klein, Carsten Gutwenger
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

#ifndef OGDF_REMOVE_REINSERT_TYPE_H
#define OGDF_REMOVE_REINSERT_TYPE_H


namespace ogdf {

	//! The postprocessing method for edge insertion algorithms.
	enum RemoveReinsertType {
		rrNone,        //!< No postprocessing.
		rrInserted,    //!< Postprocessing only with the edges that have to be inserted.
		rrMostCrossed, //!< Postprocessing with the edges involved in the most crossings.
		rrAll,         //!< Postproceesing with all edges.
		rrIncremental, //!< Full postprocessing after each edge insertion.
		rrIncInserted  //!< Postprocessing for (so far) inserted edges after each edge insertion.
	};

}

#endif
