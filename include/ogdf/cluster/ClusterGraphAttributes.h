/*
 * $Revision: 3427 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-22 10:55:42 +0200 (Mo, 22. Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Declares ClusterGraphAttributes, an extension of class
 * GraphAttributes,  to store clustergraph layout informations
 * like cluster cage positions and sizes that can be accessed
 * over the cluster/cluster ID
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

#ifndef OGDF_CLUSTER_GRAPH_ATTRIBUTES_H
#define OGDF_CLUSTER_GRAPH_ATTRIBUTES_H


#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/cluster/ClusterArray.h>


namespace ogdf {

	//! Stores information associated with a cluster.
	struct OGDF_EXPORT ClusterInfo
	{
		double m_x, m_y;	//!< position of lower left corner
		double m_w, m_h;	//!< width and height

		string m_label;		//!< cluster label

		Stroke m_stroke;	//!< stroke (style of boundary)
		Fill   m_fill;		//!< fill (style of interior)
	};


	//! Stores additional attributes of a clustered graph (like layout information).
	class OGDF_EXPORT ClusterGraphAttributes : public GraphAttributes
	{
		const ClusterGraph* m_pClusterGraph;//!< Only points to existing graphs.

	private:
		ClusterArray<ClusterInfo> m_clusterInfo;     //!< cluster attributes
		ClusterArray<string>      m_clusterTemplate; //!< Name of cluster template.

	public:
		// don't hide these inherited methods by overloading
		using GraphAttributes::x;
		using GraphAttributes::y;
		using GraphAttributes::width;
		using GraphAttributes::height;

		using GraphAttributes::label;

		using GraphAttributes::strokeType;
		using GraphAttributes::setStrokeType;
		using GraphAttributes::strokeColor;
		using GraphAttributes::strokeWidth;

		using GraphAttributes::fillPattern;
		using GraphAttributes::setFillPattern;
		using GraphAttributes::fillColor;
		using GraphAttributes::fillBgColor;

		/**
		 * @name Construction and management of attributes
		 */
		//@{

		//! Constructs cluster graph attributes for no associated graph.
		ClusterGraphAttributes() : GraphAttributes(), m_pClusterGraph(0) { }

		//! Constructs cluster graph attributes for cluster graph \a cg with attributes \a initAttributes.
		/**
		 * \remark All attributes in ClusterElement are always available.
		 */
		ClusterGraphAttributes(ClusterGraph& cg, long initAttributes = 0);

		virtual ~ClusterGraphAttributes() { }

		//! Initializes the cluster graph attributes for cluster graph \a cg with attributes \a initAttributes.
		virtual void init(ClusterGraph &cg, long initAttributes = 0);

		//! Forbidden initialization, use init(ClusterGraph &cg, long initAttributes) instead!
		virtual void init(const Graph &, long) {
			OGDF_THROW(Exception); // We need a cluster graph for initialization
		}

		//! Initializes the attributes according to \a initAttributes.
		virtual void initAtt(long initAttributes = 0) {
			GraphAttributes::initAttributes(initAttributes);
		}

		//! Returns the associated cluster graph.
		const ClusterGraph& constClusterGraph() const { return *m_pClusterGraph; }


		//@}
		/**
		 * @name Cluster attributes
		 */
		//@{

		//! Returns the x-position of cluster \a c's cage (lower left corner).
		double x(cluster c) const { return m_clusterInfo[c].m_x; }

		//! Returns the x-position of cluster \a c's cage (lower left corner).
		double& x(cluster c) { return m_clusterInfo[c].m_x; }

		//! Returns the y-position of cluster \a c's cage (lower left corner).
		double y(cluster c) const { return m_clusterInfo[c].m_y; }

		//! Returns the y-position of cluster \a c's cage (lower left corner).
		double& y(cluster c) { return m_clusterInfo[c].m_y; }

		//! Returns the width of cluster \a c.
		double width(cluster c) const { return m_clusterInfo[c].m_w; }

		//! Returns the width of cluster \a c.
		double& width(cluster c) { return m_clusterInfo[c].m_w; }

		//! Returns the height of cluster \a c.
		double height(cluster c) const { return m_clusterInfo[c].m_h; }

		//! Returns the height of cluster \a c.
		double& height(cluster c) { return m_clusterInfo[c].m_h; }

		//! Returns the stroke type of cluster \a c.
		StrokeType strokeType(cluster c) const {
			return m_clusterInfo[c].m_stroke.m_type;
		}

		//! Sets the stroke type of cluster \a c to \a st.
		void setStrokeType(cluster c, StrokeType st) {
			m_clusterInfo[c].m_stroke.m_type = st;
		}

		//! Returns the stroke color of cluster \a c.
		const Color &strokeColor(cluster c) const {
			return m_clusterInfo[c].m_stroke.m_color;
		}

		//! Returns the stroke color of cluster \a c.
		Color &strokeColor(cluster c) {
			return m_clusterInfo[c].m_stroke.m_color;
		}

		//! Returns the stroke width of cluster \a c.
		const float &strokeWidth(cluster c) const {
			return m_clusterInfo[c].m_stroke.m_width;
		}

		//! Returns the stroke width of cluster \a c.
		float &strokeWidth(cluster c) {
			return m_clusterInfo[c].m_stroke.m_width;
		}

		//! Returns the fill pattern of cluster \a c.
		FillPattern fillPattern(cluster c) const {
			return m_clusterInfo[c].m_fill.m_pattern;
		}

		//! Sets the fill pattern of cluster \a c to \a fp.
		void setFillPattern(cluster c, FillPattern fp) {
			m_clusterInfo[c].m_fill.m_pattern = fp;
		}

		//! Returns the fill color of cluster \a c.
		const Color &fillColor(cluster c) const {
			return m_clusterInfo[c].m_fill.m_color;
		}

		//! Returns the fill color of cluster \a c.
		Color &fillColor(cluster c) {
			return m_clusterInfo[c].m_fill.m_color;
		}

		//! Returns the background color of fill patterns for cluster \a c.
		const Color &fillBgColor(cluster c) const {
			return m_clusterInfo[c].m_fill.m_bgColor;
		}

		//! Returns the background color of fill patterns for cluster \a c.
		Color &fillBgColor(cluster c) {
			return m_clusterInfo[c].m_fill.m_bgColor;
		}

		//! Returns the label of cluster \a c.
		const string &label(cluster c) const {
			return m_clusterInfo[c].m_label;
		}

		//! Returns the label of cluster \a c.
		string &label(cluster c) {
			return m_clusterInfo[c].m_label;
		}

		//! Returns the template of cluster \a c.
		const string &templateCluster(cluster c) const { return m_clusterTemplate[c]; }

		//! Returns the template of cluster \a c.
		string &templateCluster(cluster c) { return m_clusterTemplate[c]; }

		//! Returns the cluster info structure of cluster \a c.
		const ClusterInfo& clusterInfo(cluster c) const { return m_clusterInfo[c]; }

		//! Returns the cluster info structure of cluster \a c.
		ClusterInfo& clusterInfo(cluster c) { return m_clusterInfo[c]; }

		//@}
		/**
		 * @name Utility functions
		 */
		//@{

		//! Returns the bounding box of the layout.
		const DRect boundingBox() const;

		//! Returns the maximum cluster index used.
		int maxClusterID() const { return m_pClusterGraph->clusterIdCount()-1; }

		//! Updates positions of cluster boundaries wrt to children and child clusters
		void updateClusterPositions(double boundaryDist = 1.0);

		//! Returns the parent cluster of node \a v.
		cluster clusterOf(node v) { return m_pClusterGraph->clusterOf(v); }

		//@}

	};


} // end namespace ogdf


#endif
