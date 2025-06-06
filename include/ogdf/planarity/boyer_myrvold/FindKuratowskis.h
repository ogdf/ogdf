/** \file
 * \brief Declaration of the class FindKuratowskis
 *
 * \author Jens Schmidt
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

#include <ogdf/basic/Array.h>
#include <ogdf/basic/ArrayBuffer.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/basic.h>

namespace ogdf {
class BoyerMyrvoldPlanar;
enum class BoyerMyrvoldEdgeType;
template<class E>
class ListPure;

/**
 * %List of externally active nodes strictly between x and y for minortypes \a B and \a E
 *
 * In case of extracting without bundles, all external paths and lists of their start-
 * and endnodes are added.
 */
struct ExternE {
	node theNode;
	SListPure<int> startnodes;
	SListPure<node> endnodes;
	SListPure<SListPure<edge>> externalPaths;
};

//! Saves information about a pertinent node w between two stopping vertices.
/** In particular links to appropriate highest-XY-Path and z-nodes are maintained
 */
struct WInfo {
	node w;

	//!  All possible base minortypes on w
	enum class MinorType {
		A = 0x0001, // minor A
		B = 0x0002, // minor B
		C = 0x0004, // minor C
		D = 0x0008, // minor D
		E = 0x0010 // minor E
	};
	int minorType;

	ArrayBuffer<adjEntry>* highestXYPath;
	ArrayBuffer<adjEntry>* zPath;

	bool pxAboveStopX;
	bool pyAboveStopY;

	SListPure<SListPure<edge>> pertinentPaths;

	SListIterator<ExternE> externEStart;
	SListIterator<ExternE> externEEnd;
	node firstExternEAfterW;
};

inline int operator&(int lhs, WInfo::MinorType rhs) { return lhs & static_cast<int>(rhs); }

inline int operator|=(int& lhs, WInfo::MinorType rhs) {
	lhs |= static_cast<int>(rhs);
	return lhs;
}

//! A Kuratowski Structure is a special graph structure containing severals subdivisions
class OGDF_EXPORT KuratowskiStructure {
	friend class FindKuratowskis;
	friend class ExtractKuratowskis;

public:
	//! Constructor
	KuratowskiStructure() { }

	//! Destructor
	~KuratowskiStructure() { }

	//! Copy constructor
	KuratowskiStructure(const KuratowskiStructure& orig) { copy(orig); }

	//! Assignment
	KuratowskiStructure& operator=(const KuratowskiStructure& orig) {
		copy(orig);
		return *this;
	}

	//! Reset all data members
	void clear();

	//! The current node to embed
	node V;
	//! DFI of the current node to embed
	int V_DFI;

	//! The root of the bicomp containing \p stopX and \p stopY
	node R;
	//! Real node of virtual node #R.
	/** This is redundant, but virtual node will be deleted later on
	 */
	node RReal;
	//! First stopping node
	node stopX;
	//! Second stopping node
	node stopY;

protected:
	//! Holds information about all pertinent nodes \a w of the bicomp containing \a V
	/** Those were not embedded because of the two stopping nodes. In addition,
	 * links to the highest-XY-path and the z-nodes of w and the minortype is saved.
	 */
	SListPure<WInfo> wNodes;

	//! The whole highestFacePath of the bicomp containing \a V
	/** The construct the highestFacePath, delete all edges of \a V except the two
	 * edges on the external face. The highestFacePath is the path starting at the
	 * first external edge along the unique face back to \a V.
	 */
	ArrayBuffer<adjEntry> highestFacePath;

	//! The appropriate paths of the highestFacePath for each wNode
	SListPure<ArrayBuffer<adjEntry>> highestXYPaths;

	//! External face path of bicomp containing \a V in direction CCW
	SListPure<adjEntry> externalFacePath;

	//! A list of all edges in all externally active paths (bundles only)
	SListPure<edge> externalSubgraph;

	//! A list of all edges in pertinent paths (bundles only)
	SListPure<edge> pertinentSubgraph;

	//! A path from the \a zNode in minortype \a D to node \a V for each highest XY-Path
	/** zNodes are cut-vertices not contained in the external face path
	 */
	SListPure<ArrayBuffer<adjEntry>> zPaths;

	//! List of externally active nodes strictly between x and y for minortypes \a B and \a E
	SListPure<ExternE> externE;

	//! List of all virtual startnodes of paths starting at #stopX (only without bundles)
	SListPure<int> stopXStartnodes;
	//! List of all virtual startnodes of paths starting at #stopY (only without bundles)
	SListPure<int> stopYStartnodes;
	//! List of all endnodes of paths starting at #stopX (only without bundles)
	SListPure<node> stopXEndnodes;
	//! List of all endnodes of paths starting at #stopY (only without bundles)
	SListPure<node> stopYEndnodes;

	//! Copies class
	void copy(const KuratowskiStructure& orig);
	//! Used in copy constructor
	void copyPointer(const KuratowskiStructure& orig, SListPure<WInfo>& list);
};

//! This class collects information about Kuratowski Subdivisions which is used for extraction later.
/** \pre Graph has to be simple.
 */
class FindKuratowskis {
public:
	//! Constructor
	explicit FindKuratowskis(BoyerMyrvoldPlanar* bm);

	//! Destructor
	~FindKuratowskis() { }

	//! Adds Kuratowski Structure on current node with root \p root and stopping nodes \p stopx and \p stopy
	void addKuratowskiStructure(const node currentNode, const node root, const node stopx,
			const node stopy);

	//! Get-method for the list of all KuratowskiStructures
	inline SListPure<KuratowskiStructure>& getAllKuratowskis() { return allKuratowskis; }

	//! Constant get-method for the list of all KuratowskiStructures
	inline const SListPure<KuratowskiStructure>& getAllKuratowskis() const {
		return allKuratowskis;
	}

	FindKuratowskis& operator=(const FindKuratowskis&) = delete;

protected:
	//! Link to class BoyerMyrvoldPlanar
	BoyerMyrvoldPlanar* pBM;

	//! Input graph
	Graph& m_g;

	//! The embedding grade
	const int& m_embeddingGrade;

	//! If true, bundles are extracted, otherwise single paths?
	const bool m_bundles;

	//! Links appropriate WInfo to node
	NodeArray<WInfo*> m_getWInfo;

	//! List of all Kuratowski Structures
	SListPure<KuratowskiStructure> allKuratowskis;

	//! Current Kuratowski Structure
	KuratowskiStructure k;

	//! Value used as marker for visited nodes etc.
	/** Used during computation of the external face path and the highest x-y-path
	 */
	int m_nodeMarker;
	//! Array maintaining visited bits on each node
	NodeArray<int> m_wasHere;

	//! Link to non-virtual vertex of a virtual Vertex.
	/** A virtual vertex has negative DFI of the DFS-Child of related non-virtual Vertex
	 */
	const NodeArray<node>& m_realVertex;

	//! The one and only DFI-NodeArray
	const NodeArray<int>& m_dfi;

	//! Returns appropriate node from given DFI
	const Array<node>& m_nodeFromDFI;

	//! Links to opposite adjacency entries on external face in clockwise resp. ccw order
	/** m_link[0]=CCW, m_link[1]=CW
	 */
	const NodeArray<adjEntry> (&m_link)[2];

	//! The adjEntry which goes from DFS-parent to current vertex
	const NodeArray<adjEntry>& m_adjParent;

	//! The DFI of the least ancestor node over all backedges
	/** If no backedge exists, the least ancestor is the DFI of that node itself
	 */
	const NodeArray<int>& m_leastAncestor;

	//! Contains the type of each edge
	EdgeArray<BoyerMyrvoldEdgeType>& m_edgeType;

	//! The lowpoint of each node
	NodeArray<int>& m_lowPoint;

	//! The highest DFI in a subtree with node as root
	const NodeArray<int>& m_highestSubtreeDFI;

	//! A list to all separated DFS-children of node
	/** The list is sorted by lowpoint values (in linear time)
	 */
	const NodeArray<ListPure<node>>& m_separatedDFSChildList;

	//! Identifies the rootnode of the child bicomp the given backedge points to
	const EdgeArray<node>& m_pointsToRoot;

	/**
	 * Stores for each (virtual) bicomp root how many backedges to its bicomp
	 * still have to be embedded. The value is set during the walkup, and it is
	 * used and decreased while embedding backedges during the walkdown.
	 */
	NodeArray<int>& m_numUnembeddedBackedgesInBicomp;

	//! Holds information, if node is the source of a backedge.
	/** This information refers to the adjEntries on the targetnode
	 * and is used during the walkdown
	 */
	NodeArray<SListPure<adjEntry>>& m_backedgeFlags;

	//! List of virtual bicomp roots, that are pertinent to the current embedded node
	NodeArray<SListPure<node>>& m_pertinentRoots;

	//! Finds root node of the bicomp containing the stopping node \p stopX
	node findRoot(node stopX) const;

	//! Extracts the highestFace Path of the bicomp containing both stopping nodes
	void extractHighestFacePath(ArrayBuffer<adjEntry>& highestFacePath, int marker);

	//! Extracts externalFacePath in direction CCW and splits highestFacePath in highestXYPaths
	void extractExternalFacePath(SListPure<adjEntry>& externalFacePath,
			const ArrayBuffer<adjEntry>& highestFacePath, int marker, int highMarker);

	//! Assign pertinent nodes to the different minortypes and extracts inner externalPaths
	void splitInMinorTypes(const SListPure<adjEntry>& externalFacePath, int marker);

	//! Extracts external subgraph from node \p stop to ancestors of node with DFI \p root (without bundles)
	void extractExternalSubgraph(const node stop, int root, SListPure<int>& externalStartnodes,
			SListPure<node>& externalEndnodes);
	//! Extracts external subgraph from node \p stop to ancestors of node with DFI \p root (with bundles)
	void extractExternalSubgraphBundles(const node stop, int root,
			SListPure<edge>& externalSubgraph, int nodeMarker);

	//! Extracts pertinent subgraph from all wNodes to \a V (without bundles)
#if 1
	void extractPertinentSubgraph(SListPure<WInfo>& W_All);
#else
	void extractPertinentSubgraph(SListPure<WInfo>& W_All, const node V);
#endif

	//! Extracts pertinent subgraph from all wNodes to \a V (with bundles)
	void extractPertinentSubgraphBundles(const SListPure<WInfo>& W_All, const node V,
			SListPure<edge>& pertinentSubgraph, int nodeMarker);
};

}
