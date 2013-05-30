/*
 * $Revision: 3505 $
 *
 * last checkin:
 *   $Author: beyer $
 *   $Date: 2013-05-16 14:49:47 +0200 (Do, 16. Mai 2013) $
 ***************************************************************/

/** \file
 * \brief Declares HypergraphAttributes storing specific attributes
 *        related to hypergraph layout drawings.
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

#ifndef OGDF_HYPERGRAPH_ATTRIBUTES_H
#define OGDF_HYPERGRAPH_ATTRIBUTES_H

#include <ogdf/hypergraph/HypergraphArray.h>
#include <ogdf/hypergraph/EdgeStandardRep.h>
#include <ogdf/basic/GraphAttributes.h>

namespace ogdf {

/**
 * \brief Stores additional attributes of a hypergraph.
 *
 * Because of different representation standards of hypergraphs (edge or
 * subset) there are two different classes for hypergraph attributes. The
 * main reason is that some edge standard atrributes are pointless in t
 * subset standard representation and vice-versa. Common attributes, currently
 * various hypernode attributes only, are pushed to the superclass
 * HypergraphAttributes.
 *
 * Similarly to GraphAttributes, attributes are simply stored in hypernode
 * or hyperedge arrays.
 */
class OGDF_EXPORT HypergraphAttributes
{
protected:

	//!< Only points to an existing hypergraph.
	const Hypergraph * m_hypergraph;

	//!< Label of a hypernode.
	HypernodeArray<string> m_label;

	//!< Coordinate x of a hypernod.e
	HypernodeArray<double> m_x;

	//!< Coordinate y of a hypernode.
	HypernodeArray<double> m_y;

	//!< Width of a hypernode bounding box.
	HypernodeArray<double> m_width;

	//!< Height of a hypernodes bounding box.
	HypernodeArray<double> m_height;

	//!< Shape of a hypernode.
	HypernodeArray<int>    m_shape;

public:

	//! Initializes new instance of class HypergraphAttributes.
	HypergraphAttributes()
	  : m_hypergraph(0)
	{
	}

	//! Initializes new instance of class HypergraphAttributes.
	HypergraphAttributes(const Hypergraph &H)
	  : m_hypergraph(&H)
	{
		m_x.init(H, 0.0);
		m_y.init(H, 0.0);
		m_width .init(H,10.0);
		m_height.init(H,10.0);
		m_label.init(H);
		m_shape.init(H, shEllipse);
	}

	//! Desctructor.
	virtual ~HypergraphAttributes()
	{
	}

	const Hypergraph & constHypergraph() const
	{
		return *m_hypergraph;
	}

	//! Returns the x-coordinate of hypernode \a v.
	const double &x(hypernode v)
	{
		return m_x[v];
	}

	//! Sets the x-coordinate of hypernode \a v.
	void setX(hypernode v, double pX)
	{
		m_x[v] = pX;
	}

	//! Returns the y-coordinate of hypernode \a v.
	const double &y(hypernode v) {
		return m_y[v];
	}

	//! Sets the y-coordinate of hypernode \a v.
	void setY(hypernode v, double pY)
	{
		m_y[v] = pY;
	}

	//! Returns the width of the bounding box of hypernode \a v.
	const double &width(hypernode v)
	{
		return m_width[v];
	}

	//! Sets the the width of hypernode \a v.
	void setWidth(hypernode v, int pWidth)
	{
		m_width[v] = pWidth;
	}

	//! Returns the height of the bounding box of hypernode \a v.
	const double &height(hypernode v)
	{
		return m_height[v];
	}

	//! Sets the the height of hypernode \a v.
	void setHeight(hypernode v, int pHeight)
	{
		m_height[v] = pHeight;
	}

	//! Returns the shape of hypernode \a v.
	int shape(hypernode v)
	{
		return m_shape[v];
	}

	//! Returns the label of hypernode \ v.
	string &label(hypernode v)
	{
		return m_label[v];
	}

};

/**
 * \brief Stores additional attributes of edge standard representation of
 *        a hypergraph.
 *
 * Since edge standard representation is in fact just an ordinary graph,
 * all its attributes are stored in wrapped instance of GraphAttributes
 * class. Some nodes in the representation are yet specific since they
 * represent hypernodes, attributes of these are driven explicitly from
 * this class.
 *
 * Superclass is declared as pure virtual to make dynamic casting possible.
 */
class OGDF_EXPORT HypergraphAttributesES : virtual public HypergraphAttributes
{
private:

	//! Wrapped graph atrributes reference.
	GraphAttributes *m_repGA;

	//! Edge standard representation reference.
	EdgeStandardRep *m_repG;

	//! The type of of edge standard representation.
	EdgeStandardType::Type m_type;

public:

	//! Initializes new instance of class HypergraphAttributes.
	HypergraphAttributesES()
	  : m_repGA(0), m_repG(0)
	{
	}

	//! Initializes new instance of class HypergraphAttributes.
	HypergraphAttributesES(const Hypergraph &pH, EdgeStandardType::Type pType = EdgeStandardType::star)
	  : HypergraphAttributes(pH)
	{
		m_repG = new EdgeStandardRep(pH, pType);

		m_repGA = new  GraphAttributes(m_repG->constGraph(),
				GraphAttributes::nodeGraphics |
				GraphAttributes::edgeGraphics);

		m_repGA->setDirected(false);
	}

	//! Destructor.
	virtual ~HypergraphAttributesES()
	{
		delete m_repGA;
		delete m_repG;
	}

	EdgeStandardType::Type type() const
	{
		return m_type;
	}

	const Graph & repGraph() const
	{
		return m_repG->constGraph();
	}

	const GraphAttributes & repGA() const
	{
		return *m_repGA;
	}

	//! Returns the type of representation node \a v.
	HypernodeElement::Type type(hypernode v)
	{
		return v->type();
	}

	//! Returns the type of representation node \a v.
	HypernodeElement::Type type(node v)
	{
		if (m_repG->hypernodeMap(v))
			return m_repG->hypernodeMap(v)->type();
		else
			return HypernodeElement::dummy;
	}

	//! Returns the x-coordinate of representation node \a v.
	const double &x(node v)
	{
		return m_repGA->x(v);
	}

	//! Sets the x-coordinate of a representation node \a v.
	void setX(node v, double pX)
	{
		if (m_repG->hypernodeMap(v))
			setX(m_repG->hypernodeMap(v), pX);
		else
			m_repGA->x(v) = pX;
	}

	//! Sets the x-coordinate of hypernode \a v.
	void setX(hypernode v, double pX)
	{
		m_x[v] = pX;
		m_repGA->x(m_repG->nodeMap(v)) = pX;
	}

	//! Returns the y-coordinate of a representation node \a v.
	const double &y(node v)
	{
		return m_repGA->y(v);
	}

	//! Sets the x-coordinate of hypernode \a v.
	void setY(hypernode v, double pY)
	{
		m_y[v] = pY;
		m_repGA->y(m_repG->nodeMap(v)) = pY;
	}

	//! Sets the y-coordinate of a representation node \a v.
	void setY(node v, double pY)
	{
		if (m_repG->hypernodeMap(v))
			setY(m_repG->hypernodeMap(v), pY);
		else
			m_repGA->y(v) = pY;
	}

	//! Returns the width of a representation node \a v.
	const double &width(node v)
	{
		return m_repGA->width(v);
	}

	//! Sets the the width of hypernode \a v.
	void setWidth(hypernode v, double pWidth)
	{
		m_width[v] = pWidth;
		m_repGA->width(m_repG->nodeMap(v)) = pWidth;
	}

	//! Sets the the width of a representation node \a v.
	void setWidth(node v, double pWidth)
	{
		if (m_repG->hypernodeMap(v))
			setWidth(m_repG->hypernodeMap(v), pWidth);
		else
			m_repGA->width(v) = pWidth;
	}

	//! Returns the height of a representation node \a v.
	const double &height(node v)
	{
		return m_repGA->height(v);
	}

	//! Sets the the height of hypernode \a v.
	void setHeight(hypernode v, double pHeight)
	{
		m_height[v] = pHeight;
		m_repGA->height(m_repG->nodeMap(v)) = pHeight;
	}

	//! Sets the the height of a representation node \a v.
	void setHeight(node v, double pHeight)
	{
		if (m_repG->hypernodeMap(v))
			setHeight(m_repG->hypernodeMap(v), pHeight);
		else
			m_repGA->height(v) = pHeight;
	}

	//! Returns the list of bend points of edge \a e.
	DPolyline &bends(edge e)
	{
		return m_repGA->bends(e);
	}

	////! Writes the hypergraph (edge standard representation) into SVG format.
	//void writeSVG(ostream &os, int fontSize, const string &fontColor) const
	//{
	//    m_repGA->writeSVG(os, fontSize, fontColor);
	//}
	//
	////! Writes the hypergraph (edge standard representation) into SVG format.
	//void writeSVG(const char *fileName, int fontSize, const string &fontColor) const
	//{
	//    m_repGA->writeSVG(fileName, fontSize, fontColor);
	//}

	void clearAllBends()
	{
		m_repGA->clearAllBends();
	}

	//! Removes unnecessary bend points in orthogonal segements.
	void removeUnnecessaryBendsHV()
	{
		m_repGA->removeUnnecessaryBendsHV();
	}

	//! Returns the bounding box of the hypergraph.
	const DRect boundingBox() const
	{
		return m_repGA->boundingBox();
	}

};

} // end namespace ogdf

#endif
