/*
 * $Revision: 3505 $
 *
 * last checkin:
 *   $Author: beyer $
 *   $Date: 2013-05-16 14:49:47 +0200 (Do, 16. Mai 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration and a partial implementation of a Hypergraph class
 *        partly based on the original classes for handling hypergraphs
 *        written by Martin Gronemann.
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

#ifndef OGDF_HYPERGRAPH_H
#define OGDF_HYPERGRAPH_H

#include <ogdf/basic/List.h>
#include <ogdf/basic/GraphList.h>

//! Iteration over all adjacency list entries \a adj of a graph element \a ge.
#define forall_adj_elements(adj,ge) for((adj)=(v)->firstAdj();(adj);(adj)=(adj)->succ())

//! Iteration over all hypernodes \a v of hypergraph \a H.
#define forall_hypernodes(v,H) for((v)=(H).firstHypernode(); (v); (v)=(v)->succ())

//! Iteration over all hypernodes \a v of hypergraph \a H in reverse order.
#define forall_rev_hypernodes(v,H) for((v)=(H).lastHypernode(); (v); (v)=(v)->pred())

//! Iteration over all hyperedges \a e of hypergraph \a H.
#define forall_hyperedges(e,H) for((e)=(H).firstHyperedge(); (e); (e)=(e)->succ())

//! Iteration over all hyperedges \a e of hypergraph \a G in reverse order.
#define forall_rev_hyperedges(e,H) for((e)=(H).lastHyperedge(); (e); (e)=(e)->pred())

namespace ogdf {

class OGDF_EXPORT Hypergraph;
class OGDF_EXPORT HypernodeElement;
class OGDF_EXPORT HyperedgeElement;
class OGDF_EXPORT AdjHypergraphElement;

//!< The type of hypernodes.
typedef HypernodeElement *hypernode;

//!< The type of hyperedges.
typedef HyperedgeElement *hyperedge;

//!< The type of adjacency entries.
typedef AdjHypergraphElement *adjHypergraphEntry;

//! Class for adjacency list elements.
/**
 * Adjacency list elements represent the occurrence of hypernodes of a
 * hyperedge or hyperedges of a hypernode adjacency lists.
 */
class OGDF_EXPORT AdjHypergraphElement : private GraphElement
{
	friend class Hypergraph;
	friend class GraphListBase;
	friend class GraphList<AdjHypergraphElement>;

private:

	//!< The associated hyperedge or hypernode.
	GraphElement *m_element;

	//!< The corresponding adjacency entry.
	/**
	 * Note: For instance, if this AdjHypergraphElement is associated with
	 * a hypernode v, its element is an hyperedge e (incident with v) then
	 * the twin AdjHypergraphElement is the one associated with e and
	 * containing element v. This is a different from m_twin in Graph_d.h!
	 */
	adjHypergraphEntry m_twin;

	//!< The (unique) index of the adjacency entry.
	int m_index;

	//! Constructs an adjacency element for a given hyper{node,edge}.
	AdjHypergraphElement(GraphElement *pElement)
	  : m_element(pElement)
	{ }

	//! Constructs an adjacency entry for a given hyper{node,edge} and index.
	AdjHypergraphElement(GraphElement *pElement, int pIndex)
	  : m_element(pElement), m_index(pIndex)
	{ }

public:

	//! Returns the index of this adjacency element.
	int index() const
	{
		return m_index;
	}

	//! Returns the element associated with this adjacency entry.
	/**
	 * It is always a pointer to an instance of HypernodeElement or
	 * HyperedgeElement, you are expected to explicitly overtype it
	 * to the correct type.
	 */
	GraphElement * element() const
	{
		return m_element;
	}

	//! Returns the pointer to a twin adjacency list.
	adjHypergraphEntry twin() const
	{
		return m_twin;
	}

	//! Returns the successor in the adjacency list.
	adjHypergraphEntry succ() const
	{
		return (adjHypergraphEntry) m_next;
	}

	//! Returns the predecessor in the adjacency list.
	adjHypergraphEntry pred() const
	{
		return (adjHypergraphEntry) m_prev;
	}

	//! Returns the cyclic successor in the adjacency list.
	adjHypergraphEntry cyclicSucc() const;

	//! Returns the cyclic predecessor in the adjacency list.
	adjHypergraphEntry cyclicPred() const;

	OGDF_NEW_DELETE;

}; // class AdjHypegraphElement

//! Class for the representation of hyperedges.
class OGDF_EXPORT HyperedgeElement : private GraphElement
{
	friend class Hypergraph;
	friend class GraphListBase;
	friend class GraphList<HyperedgeElement>;

public:

	//! The type of hyperedges (for not there is only a single type).
	enum Type {
		normal = 0x100001
	};

private:

	//!< The adjacency list of the hyperedge.
	GraphList<AdjHypergraphElement> m_adjHypernodes;

	//! The (unique) index of the hyperedge.
	int m_index;

	//! The number of incidend hypernodes.
	int m_cardinality;

	//! The type of the hyperedge.
	Type m_type;

	//! The hypergraph containing the hyperedge (if any).
	Hypergraph *m_hypergraph;

	//! Constructs an hyperedge element between hypernodes.
	/**
	 * @param pIndex is the index of the hyperedge.
	 */
	HyperedgeElement(int pIndex)
	  : m_index(pIndex), m_cardinality(0), m_type(HyperedgeElement::normal), m_hypergraph(0) { }

	//! Constructs an hyperedge element between hypernodes.
	/**
	 * @param pIndex is the index of the hyperedge.
	 * @param pType is the type of the hyperedge.
	 */
	HyperedgeElement(int pIndex, Type pType)
	  : m_index(pIndex), m_cardinality(0), m_type(pType), m_hypergraph(0) { }

public:

	//! Returns the index of a hyperedge.
	int index() const
	{
		return m_index;
	}

	//! Returns the number of incident hypernodes.
	int cardinality() const
	{
		return m_cardinality;
	}

	//! Returns the type of the hyperedge.
	Type type() const
	{
		return m_type;
	}

	//! Returns the hypergraph containing the hyperedge.
	Hypergraph * hypergraph() const
	{
		return m_hypergraph;
	}

	//! Returns the first entry in the adjaceny list.
	adjHypergraphEntry firstAdj() const
	{
		return m_adjHypernodes.begin();
	}

	//! Returns the last entry in the adjacency list.
	adjHypergraphEntry lastAdj () const
	{
		return m_adjHypernodes.rbegin();
	}

	//! Returns the incident hypernodes of the hyperedge.
	GraphList<AdjHypergraphElement> incidentHypernodes() const
	{
		return m_adjHypernodes;
	}

	//! Returns a list with all incident hypernodes of the hyperedge.
	template<class NODELIST> void allHypernodes(NODELIST &hypernodes) const
	{
		hypernodes.clear();
		for (adjHypergraphEntry adj = firstAdj(); adj; adj = adj->succ())
			hypernodes.pushBack((hypernode) adj->element());
	}

	//! Returns true iff \a v is incident to the hyperedge.
	bool incident(hypernode v) const
	{
		for (adjHypergraphEntry adj = firstAdj(); adj; adj = adj->succ())
			if (((hypernode) adj->element()) == v)
				return true;
		return false;
	}

	//! Returns the successor in the list of all hyperedges.
	hyperedge succ() const
	{
		return (hyperedge) m_next;
	}

	//! Returns the predecessor in the list of all hyperedges.
	hyperedge pred() const
	{
		return (hyperedge) m_prev;
	}

	//! Equality operator.
	bool operator==(const hyperedge e) const
	{
		return e->index() == m_index && e->hypergraph() == m_hypergraph;
	}

	friend ostream & operator<<(ostream &os, ogdf::hyperedge e);

	OGDF_NEW_DELETE;

}; // class EdgeElement

//! Class for the representation of hypernodes.
class OGDF_EXPORT HypernodeElement : private GraphElement
{
	friend class Hypergraph;
	friend class GraphListBase;
	friend class GraphList<HypernodeElement>;

public:

	//! The type of hypernodes.
	enum Type {
		normal = 0x0000001, //! Default type.
		dummy  = 0x0000002, //! Temporary hypernode.
		OR     = 0x0000003, //! Electric circuit: OR gate.
		BUF    = 0x0000004, //! Electric circuit: Buffer gate (iscas85).
		AND    = 0x0000005, //! Electric circuit: AND gate.
		NOR    = 0x0000006, //! Electric circuit: NOR gate.
		NOT    = 0x0000007, //! Electric circuit: NOT gate.
		XOR    = 0x0000008, //! Electric circuit: XOR gate.
		DFF    = 0x0000009, //! Electric circuit: D-Flip-Flop gate (max500nodes).
		NAND   = 0x0000010, //! Electric circuit: NAND gate.
		INPUT  = 0x0000011, //! Electric circuit: Input.
		OUTPUT = 0x0000012  //! Electric circuit: Output.
	};

private:

	//!< The adjacency list of the hypernode.
	GraphList<AdjHypergraphElement> m_adjHyperedges;

	//!< The (unique) index of the hypernode.
	int m_index;

	//!< The number of incident hyperedges.
	int m_degree;

	//!< The type of the hypernode.
	Type m_type;

	//!< The hypergraph containing the hypernode (if any).
	Hypergraph * m_hypergraph;

	//!< Constructor.
	HypernodeElement(int pIndex)
	  : m_index(pIndex), m_degree(0), m_type(normal), m_hypergraph(0)
	{
	}

	//!< Constructor.
	HypernodeElement(int pIndex, Type pType)
	  : m_index(pIndex), m_degree(0), m_type(pType), m_hypergraph(0)
	{
	}

public:

	//! Returns the (unique) hypernode index.
	int index() const
	{
		return m_index;
	}

	//! Returns the hypernode degree.
	int degree() const
	{
		return m_degree;
	}

	//! Returns the hypergraph containing the hypernode.
	Hypergraph * hypergraph() const
	{
		return m_hypergraph;
	}

	//! Returns the type of hypernode.
	Type type() const
	{
		return m_type;
	}

	//! Sets the type of hypernode.
	void type(Type pType)
	{
		m_type = pType;
	}

	//! Returns the first entry in the adjaceny list.
	adjHypergraphEntry firstAdj() const
	{
		return m_adjHyperedges.begin();
	}

	//! Returns the last entry in the adjacency list.
	adjHypergraphEntry lastAdj() const
	{
		return m_adjHyperedges.rbegin();
	}

	//! Returns a list with all incident hyperedges of the hypernode.
	template<class NODELIST> void allHyperedges(NODELIST &hyperedges) const
	{
		hyperedges.clear();
		for (adjHypergraphEntry adj = firstAdj(); adj; adj = adj->succ())
			hyperedges.pushBack((hyperedge) adj->element());
	}

	//! Returns true iff \a v is adjacent to the hypernode.
	bool adjacent(hypernode v) const
	{
		for (adjHypergraphEntry adj = firstAdj(); adj; adj = adj->succ())
			if (((hyperedge) adj->element())->incident(v))
				return true;
		return false;
	}

	//! Returns the successor in the list of all hypernodes.
	hypernode succ() const
	{
		return (hypernode) m_next;
	}

	//! Returns the predecessor in the list of all hypernodes.
	hypernode pred() const
	{
		return (hypernode) m_prev;
	}

	//! Equality operator.
	bool operator==(const hypernode v) const
	{
		return v->index() == m_index && v->hypergraph() == m_hypergraph;
	}

	OGDF_NEW_DELETE;

}; // class HypernodeElement

class HypergraphArrayBase;
template<class T> class HypernodeArray;
template<class T> class HyperedgeArray;
class OGDF_EXPORT HypergraphObserver;

class OGDF_EXPORT Hypergraph
{
	//!< The list of all hypernodes.
	GraphList<HypernodeElement> m_hypernodes;

	//!< The list of all hyperedges.
	GraphList<HyperedgeElement> m_hyperedges;

	//!< The number of hypernodes in the hypergraph.
	int m_nHypernodes;

	//!< The number of hyperedges in the hypergraph.
	int m_nHyperedges;

	//!< The Index that will be assigned to the next created hypernode.
	int m_hypernodeIdCount;

	//!< The Index that will be assigned to the next created hyperedge.
	int m_hyperedgeIdCount;

	//!< The current table size of hypernode arrays within the hypergraph.
	int m_hypernodeArrayTableSize;

	//!< The current table size of hyperedge arrays within the hypergraph.
	int m_hyperedgeArrayTableSize;

	//!< The registered hypergraph arrays & observers.
	mutable ListPure<HypergraphArrayBase *> m_hypernodeArrays;
	mutable ListPure<HypergraphArrayBase *> m_hyperedgeArrays;
	mutable ListPure<HypergraphObserver *> m_observers;

public:

	//! Constructs an empty hypergraph.
	Hypergraph();

	//! Constructs a hypergraph that is a copy of \a H.
	Hypergraph(const Hypergraph &H);

	//! Destructor.
	~Hypergraph();

	//! Returns true iff the hypergraph is empty (ie. contains no hypernodes).
	bool empty() const
	{
		return m_nHypernodes == 0;
	}

	//! Returns the list of all hypernodes.
	GraphList<HypernodeElement> hypernodes() const
	{
		return m_hypernodes;
	}

	//! Returns the list of all hyperedges.
	GraphList<HyperedgeElement> hyperedges() const
	{
		return m_hyperedges;
	}

	//! Returns the number of hypernodes in the hypergraph.
	int numberOfHypernodes() const
	{
		return m_nHypernodes;
	}

	//! Returns the number of hyperedges in the hypergraph.
	int numberOfHyperedges() const
	{
		return m_nHyperedges;
	}

	//! Returns the largest used hypernode index.
	int maxHypernodeIndex() const
	{
		return m_hypernodeIdCount - 1;
	}

	//! Returns the largest used hyperedge index.
	int maxHyperedgeIndex() const
	{
		return m_hyperedgeIdCount - 1;
	}

	//! Returns the first hypernode in the list of all hypernodes.
	hypernode firstHypernode() const
	{
		return m_hypernodes.begin ();
	}

	//! Returns the last hypernode in the list of all hypernodes.
	hypernode lastHypernode () const
	{
		return m_hypernodes.rbegin();
	}

	//! Returns the first hyperedge in the list of all hyperedges.
	hyperedge firstHyperedge() const
	{
		return m_hyperedges.begin ();
	}

	//! Returns the last hyperedge in the list of all hyperedges.
	hyperedge lastHyperEdge () const
	{
		return m_hyperedges.rbegin();
	}

	//! Returns the table size of hypernode arrays with the hypergraph.
	int hypernodeArrayTableSize() const
	{
		return m_hypernodeArrayTableSize;
	}

	//! Returns the table size of hyperedge arrays within the hypergraph.
	int hyperedgeArrayTableSize() const
	{
		return m_hyperedgeArrayTableSize;
	}

	//! Creates a new hypernode and returns it.
	hypernode newHypernode();

	//! Creates a new hypernode with given \a index and returns it.
	hypernode newHypernode(int pIndex);

	//! Creates a new hypernode with given \a type and returns it.
	hypernode newHypernode(HypernodeElement::Type pType);

	//! Creates a new hypernode with given \a indexn and \a type and returns it.
	hypernode newHypernode(int pIndex, HypernodeElement::Type pType);

	//! Creates a new hyperedge btween \a hypernodes and returns it.
	/**
	 * @param hypernodes are the hypernodes of the newly created hyperedge.
	 * @return the newly created hyperedge.
	 */
	hyperedge newHyperedge(List<hypernode> &hypernodes);

	//! Creates a new hyperedge between \a hypernodes and returns it.
	/**
	 * @param pIndex is the unique index of the newly created hyperedge.
	 * @param hypernodes are the hypernodes of the newly created hyperedge.
	 * @return the newly created hyperedge.
	 */
	hyperedge newHyperedge(int pIndex, List<hypernode> &hypernodes);

	//! Removes hypernode \a v and all incident hyperedges from the hypergraph.
	/**
	 * @param v is the hypernode that will be deleted.
	 */
	void delHypernode(hypernode v);

	//! Removes hyperedge \a e from the hypergraph.
	/**
	 * @param e is the hyperegde that will be deleted.
	 */
	void delHyperedge(hyperedge e);

	//! Removes all hypernodes and all hyperedges from the hypergraph.
	void clear();

	//! Returns a randomly chosen hypergraph.
	hypernode randomHypernode() const;

	//! Returns a randomly chosen hyperedge.
	hyperedge randomHyperedge() const;

	//! Returns a list with all hypernodes of the hypergraph.
	template<class LIST> void allHypernodes(LIST &hypernodes) const
	{
		hypernodes.clear();
		for (hypernode v = m_hypernodes.begin(); v; v = v->succ())
			hypernodes.pushBack(v);
	}

	//! Returns a list with all hyperedges of the hypergraph.
	template<class LIST> void allHyperedges(LIST &hyperedges) const
	{
		hyperedges.clear();
		for (hyperedge e = m_hyperedges.begin(); e; e = e->succ())
			hyperedges.pushBack(e);
	}

	//! Reads hypergraph in bench format from the input stream.
	void readBenchHypergraph(istream &is);

	//! Reads hypergraph in bench format from the file.
	void readBenchHypergraph(const char *filename);

	//! Reads hypergraph in pla format from the input stream.
	void readPlaHypergraph(istream &is);

	//! Reads hypergraph in pla format from the file.
	void loadPlaHypergraph(const char *fileName);

	//! Checks the consistency of the data structure.
	bool consistency() const;

	//! Registers a node array.
	ListIterator<HypergraphArrayBase *>
		registerHypernodeArray(HypergraphArrayBase *pHypernodeArray) const;

	ListIterator<HypergraphArrayBase *>
		registerHyperedgeArray(HypergraphArrayBase *pHyperedgeArray) const;

	//! Registers a hypergraph observer (e.g. a EdgeStandardRep).
	ListIterator<HypergraphObserver *>
		registerObserver(HypergraphObserver *pObserver) const;

	//! Unregisters a hypernode array.
	void unregisterHypernodeArray(ListIterator<HypergraphArrayBase *> it) const;

	//! Unregisters an hyperedge array.
	void unregisterHyperedgeArray(ListIterator<HypergraphArrayBase *> it) const;

	//! Unregisters a hypergraph observer.
	void unregisterObserver(ListIterator<HypergraphObserver *> it) const;

	Hypergraph &operator=(const Hypergraph &H);

	friend ostream & operator<<(ostream &os, ogdf::Hypergraph &H);

	friend istream & operator>>(istream &is, ogdf::Hypergraph &H);

	OGDF_MALLOC_NEW_DELETE;

private:

	void initArrays();

	void initObservers();

	int nextEntry(char *buffer, int from, string stop);

	HypernodeElement::Type gateType(string gate);

}; // class Hypergraph

} // end of namespace ogdf

#endif /* HYPERGRAPH_H_ */
