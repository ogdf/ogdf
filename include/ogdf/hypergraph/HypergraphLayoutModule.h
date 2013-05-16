/*
 * $Revision: 3505 $
 *
 * last checkin:
 *   $Author: beyer $
 *   $Date: 2013-05-16 14:49:47 +0200 (Do, 16. Mai 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration of an interface for hypergraph layout
 *        algorithms. Any hypergraph layout must follow this prescription.
 *
 * \author Ondrej Moris
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

#ifndef OGDF_HYPERGRAPH_LAYOUT_MODULE_H
#define OGDF_HYPERGRAPH_LAYOUT_MODULE_H

#include <ogdf/hypergraph/HypergraphAttributes.h>

namespace ogdf {

/**
 * \brief Interface of hypergraph layout algorithms.
 *
 */
class OGDF_EXPORT HypergraphLayoutModule
{
public:

	//! Initializes a layout module.
	HypergraphLayoutModule()
	{
	}

	virtual ~HypergraphLayoutModule()
	{
	}

	/**
	 * \brief Computes a layout of hypergraph given by \a HA.
	 *
	 * This method is the actual algorithm call and must be implemented by
	 * derived classes.
	 * @param HA is the input hypergraph attributes class.
	 */
	virtual void call(HypergraphAttributes &HA) = 0;

	/**
	 * \brief Computes a layout of a hypergraph given by \a HA.
	 *
	 * @param HA is the input hypergraph attributes class.
	 */
	void operator()(HypergraphAttributes &HA)
	{
		call(HA);
	}

	OGDF_MALLOC_NEW_DELETE;
};

} // end namespace ogdf

#endif
