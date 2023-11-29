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

#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/Observer.h>
#include <ogdf/basic/RegisteredArray.h>

//! Iteration over all adjacency list entries \p adj of a graph element \p ge.
#define forall_adj_elements(adj, ge) for ((adj) = (v)->firstAdj(); (adj); (adj) = (adj)->succ())

//! Iteration over all hypernodes \p v of hypergraph \p H.
#define forall_hypernodes(v, H) for ((v) = (H).firstHypernode(); (v); (v) = (v)->succ())

//! Iteration over all hypernodes \p v of hypergraph \p H in reverse order.
#define forall_rev_hypernodes(v, H) for ((v) = (H).lastHypernode(); (v); (v) = (v)->pred())

//! Iteration over all hyperedges \p e of hypergraph \p H.
#define forall_hyperedges(e, H) for ((e) = (H).firstHyperedge(); (e); (e) = (e)->succ())

//! Iteration over all hyperedges \p e of hypergraph \p H in reverse order.
#define forall_rev_hyperedges(e, H) for ((e) = (H).lastHyperedge(); (e); (e) = (e)->pred())

namespace ogdf {

class OGDF_EXPORT Hypergraph;
class OGDF_EXPORT HypernodeElement;
class OGDF_EXPORT HyperedgeElement;
class OGDF_EXPORT AdjHypergraphElement;

//! The type of hypernodes.
using hypernode = HypernodeElement*;

//! The type of hyperedges.
using hyperedge = HyperedgeElement*;

//! The type of adjacency entries.
using adjHypergraphEntry = AdjHypergraphElement*;

//! Class for adjacency list elements.
/**
 * Adjacency list elements represent the occurrence of hypernodes of a
 * hyperedge or hyperedges of a hypernode adjacency lists.
 */
class OGDF_EXPORT AdjHypergraphElement : private internal::GraphElement {
	friend class Hypergraph;
	friend class GraphListBase;
	friend class internal::GraphList<AdjHypergraphElement>;

private:
	//! The associated hyperedge or hypernode.
	GraphElement* m_element;

	//! The corresponding adjacency entry.
	/**
	 * Note: For instance, if this AdjHypergraphElement is associated with
	 * a hypernode v, its element is an hyperedge e (incident with v) then
	 * the twin AdjHypergraphElement is the one associated with e and
	 * containing element v. This is a different from m_twin in Graph_d.h!
	 */
	adjHypergraphEntry m_twin;

	//! The (unique) index of the adjacency entry.
	int m_index;

	//! Constructs an adjacency element for a given hyper{node,edge}.
	explicit AdjHypergraphElement(GraphElement* pElement)
		: m_element(pElement), m_twin(nullptr), m_index(0) { }

	//! Constructs an adjacency entry for a given hyper{node,edge} and index.
	AdjHypergraphElement(GraphElement* pElement, int pIndex)
		: m_element(pElement), m_twin(nullptr), m_index(pIndex) { }

public:
	//! Returns the index of this adjacency element.
	int index() const { return m_index; }

	//! Returns the element associated with this adjacency entry.
	/**
	 * It is always a pointer to an instance of HypernodeElement or
	 * HyperedgeElement, you are expected to explicitly overtype it
	 * to the correct type.
	 */
	GraphElement* element() const { return m_element; }

	//! Returns the pointer to a twin adjacency list.
	adjHypergraphEntry twin() const { return m_twin; }

	//! Returns the successor in the adjacency list.
	adjHypergraphEntry succ() const { return static_cast<adjHypergraphEntry>(m_next); }

	//! Returns the predecessor in the adjacency list.
	adjHypergraphEntry pred() const { return static_cast<adjHypergraphEntry>(m_prev); }

	//! Returns the cyclic successor in the adjacency list.
	adjHypergraphEntry cyclicSucc() const;

	//! Returns the cyclic predecessor in the adjacency list.
	adjHypergraphEntry cyclicPred() const;

	OGDF_NEW_DELETE;
};

//! Class for the representation of hyperedges.
class OGDF_EXPORT HyperedgeElement : private internal::GraphElement {
	friend class Hypergraph;
	friend class GraphListBase;
	friend class internal::GraphList<HyperedgeElement>;

private:
	//! The adjacency list of the hyperedge.
	internal::GraphList<AdjHypergraphElement> m_adjHypernodes;

	//! The (unique) index of the hyperedge.
	int m_index;

	//! The number of incidend hypernodes.
	int m_cardinality;

	//! The hypergraph containing the hyperedge (if any).
	Hypergraph* m_hypergraph;

	//! Constructs an hyperedge element between hypernodes.
	/**
	 * @param pIndex is the index of the hyperedge.
	 */
	explicit HyperedgeElement(int pIndex)
		: m_index(pIndex), m_cardinality(0), m_hypergraph(nullptr) { }

public:
	//! Returns the index of a hyperedge.
	int index() const { return m_index; }

	//! Returns the number of incident hypernodes.
	int cardinality() const { return m_cardinality; }

	//! Returns the hypergraph containing the hyperedge.
	Hypergraph* hypergraph() const { return m_hypergraph; }

	//! Returns the first entry in the adjaceny list.
	adjHypergraphEntry firstAdj() const { return m_adjHypernodes.head(); }

	//! Returns the last entry in the adjacency list.
	adjHypergraphEntry lastAdj() const { return m_adjHypernodes.tail(); }

	//! Returns the incident hypernodes of the hyperedge.
	internal::GraphList<AdjHypergraphElement> incidentHypernodes() const { return m_adjHypernodes; }

	//! Returns a list with all incident hypernodes of the hyperedge.
	template<class NODELIST>
	void allHypernodes(NODELIST& hypernodes) const {
		hypernodes.clear();
		for (adjHypergraphEntry adj = firstAdj(); adj; adj = adj->succ()) {
			hypernodes.pushBack(reinterpret_cast<hypernode>(adj->element()));
		}
	}

	//! Returns true iff \p v is incident to the hyperedge.
	bool incident(hypernode v) const {
		for (adjHypergraphEntry adj = firstAdj(); adj; adj = adj->succ()) {
			if (reinterpret_cast<hypernode>(adj->element()) == v) {
				return true;
			}
		}
		return false;
	}

	//! Returns the successor in the list of all hyperedges.
	hyperedge succ() const { return static_cast<hyperedge>(m_next); }

	//! Returns the predecessor in the list of all hyperedges.
	hyperedge pred() const { return static_cast<hyperedge>(m_prev); }

	//! Equality operator.
	bool operator==(const hyperedge e) const {
		return e->index() == m_index && e->hypergraph() == m_hypergraph;
	}

	friend std::ostream& operator<<(std::ostream& os, ogdf::hyperedge e);

	OGDF_NEW_DELETE;
};

//! Class for the representation of hypernodes.
class OGDF_EXPORT HypernodeElement : private internal::GraphElement {
	friend class Hypergraph;
	friend class GraphListBase;
	friend class internal::GraphList<HypernodeElement>;

public:
	//! The type of hypernodes.
	enum class Type {
		normal = 0x0000001, //! Default type.
		dummy = 0x0000002, //! Temporary hypernode.
		OR = 0x0000003, //! Electric circuit: OR gate.
		BUF = 0x0000004, //! Electric circuit: Buffer gate (iscas85).
		AND = 0x0000005, //! Electric circuit: AND gate.
		NOR = 0x0000006, //! Electric circuit: NOR gate.
		NOT = 0x0000007, //! Electric circuit: NOT gate.
		XOR = 0x0000008, //! Electric circuit: XOR gate.
		DFF = 0x0000009, //! Electric circuit: D-Flip-Flop gate (max500nodes).
		NAND = 0x0000010, //! Electric circuit: NAND gate.
		INPUT = 0x0000011, //! Electric circuit: Input.
		OUTPUT = 0x0000012 //! Electric circuit: Output.
	};

private:
	//! The adjacency list of the hypernode.
	internal::GraphList<AdjHypergraphElement> m_adjHyperedges;

	//! The (unique) index of the hypernode.
	int m_index;

	//! The number of incident hyperedges.
	int m_degree;

	//! The type of the hypernode.
	Type m_type;

	//! The hypergraph containing the hypernode (if any).
	Hypergraph* m_hypergraph;

	//! Constructor.
	explicit HypernodeElement(int pIndex)
		: m_index(pIndex), m_degree(0), m_type(Type::normal), m_hypergraph(nullptr) { }

	//! Constructor.
	HypernodeElement(int pIndex, Type pType)
		: m_index(pIndex), m_degree(0), m_type(pType), m_hypergraph(nullptr) { }

public:
	//! Returns the (unique) hypernode index.
	int index() const { return m_index; }

	//! Returns the hypernode degree.
	int degree() const { return m_degree; }

	//! Returns the hypergraph containing the hypernode.
	Hypergraph* hypergraph() const { return m_hypergraph; }

	//! Returns the type of hypernode.
	Type type() const { return m_type; }

	//! Sets the type of hypernode.
	void type(Type pType) { m_type = pType; }

	//! Returns the first entry in the adjaceny list.
	adjHypergraphEntry firstAdj() const { return m_adjHyperedges.head(); }

	//! Returns the last entry in the adjacency list.
	adjHypergraphEntry lastAdj() const { return m_adjHyperedges.tail(); }

	//! Returns a list with all incident hyperedges of the hypernode.
	template<class NODELIST>
	void allHyperedges(NODELIST& hyperedges) const {
		hyperedges.clear();
		for (adjHypergraphEntry adj = firstAdj(); adj; adj = adj->succ()) {
			hyperedges.pushBack(reinterpret_cast<hyperedge>(adj->element()));
		}
	}

	//! Returns true iff \p v is adjacent to the hypernode.
	bool adjacent(hypernode v) const {
		for (adjHypergraphEntry adj = firstAdj(); adj; adj = adj->succ()) {
			if (reinterpret_cast<hyperedge>(adj->element())->incident(v)) {
				return true;
			}
		}
		return false;
	}

	//! Returns the successor in the list of all hypernodes.
	hypernode succ() const { return static_cast<hypernode>(m_next); }

	//! Returns the predecessor in the list of all hypernodes.
	hypernode pred() const { return static_cast<hypernode>(m_prev); }

	//! Equality operator.
	bool operator==(const hypernode v) const {
		return v->index() == m_index && v->hypergraph() == m_hypergraph;
	}

	OGDF_NEW_DELETE;
};

//! Registry for nodes and edges of a hypergraph.
template<typename Key>
class HypergraphRegistry
	: public RegistryBase<Key*, HypergraphRegistry<Key>, internal::GraphIterator<Key*>> {
	Hypergraph* m_pGraph;
	int* m_nextKeyIndex;

public:
	using iterator = internal::GraphIterator<Key*>;

	//! Constructor.
	HypergraphRegistry(Hypergraph* graph, int* nextKeyIndex)
		: m_pGraph(graph), m_nextKeyIndex(nextKeyIndex) { }

	static inline int keyToIndex(Key* key) { return key->index(); }

	bool isKeyAssociated(Key* key) const {
		if (key == nullptr) {
			return false;
		}
#ifdef OGDF_DEBUG
		if (key->hypergraph() == m_pGraph) {
			OGDF_ASSERT(keyToIndex(key) < this->getArraySize());
			return true;
		} else {
			return false;
		}
#else
		return true;
#endif
	}

	int calculateArraySize(int add) const { return calculateTableSize(*m_nextKeyIndex + add); }

	int maxKeyIndex() const { return (*m_nextKeyIndex) - 1; }

	//! Returns a pointer to the associated hypergraph
	Hypergraph* graphOf() const { return m_pGraph; }
};

OGDF_EXPORT HypergraphRegistry<HypernodeElement>::iterator begin(
		const HypergraphRegistry<HypernodeElement>& self);

OGDF_EXPORT HypergraphRegistry<HypernodeElement>::iterator end(
		const HypergraphRegistry<HypernodeElement>& self);

OGDF_EXPORT HypergraphRegistry<HyperedgeElement>::iterator begin(
		const HypergraphRegistry<HyperedgeElement>& self);

OGDF_EXPORT HypergraphRegistry<HyperedgeElement>::iterator end(
		const HypergraphRegistry<HyperedgeElement>& self);

//! RegisteredArray for nodes and edges of a hypergraph.
template<typename Key, typename Value, bool WithDefault, typename Registry = HypergraphRegistry<Key>>
class HypergraphRegisteredArray : public RegisteredArray<Registry, Value, WithDefault, Hypergraph> {
	using RA = RegisteredArray<Registry, Value, WithDefault, Hypergraph>;

public:
	using RA::RA;

	//! Returns a pointer to the associated hypergraph.
	Hypergraph* hypergraphOf() const {
		if (RA::registeredAt() == nullptr) {
			return nullptr;
		} else {
			return RA::registeredAt()->graphOf();
		}
	}
};

#define OGDF_DECL_REG_ARRAY_TYPE(v, c) HypergraphRegisteredArray<HypernodeElement, v, c>
OGDF_DECL_REG_ARRAY(HypernodeArray)
#undef OGDF_DECL_REG_ARRAY_TYPE

#define OGDF_DECL_REG_ARRAY_TYPE(v, c) HypergraphRegisteredArray<HyperedgeElement, v, c>
OGDF_DECL_REG_ARRAY(HyperedgeArray)
#undef OGDF_DECL_REG_ARRAY_TYPE

class OGDF_EXPORT HypergraphObserver;

class OGDF_EXPORT Hypergraph : public Observable<HypergraphObserver, Hypergraph> {
	//! The registered hypernode arrays
	HypergraphRegistry<HypernodeElement> m_regHypernodeArrays;

	//! The registered hyperedge arrays
	HypergraphRegistry<HyperedgeElement> m_regHyperedgeArrays;

	//! The list of all hypernodes.
	internal::GraphList<HypernodeElement> m_hypernodes;

	//! The list of all hyperedges.
	internal::GraphList<HyperedgeElement> m_hyperedges;

	//! The number of hypernodes in the hypergraph.
	int m_nHypernodes;

	//! The number of hyperedges in the hypergraph.
	int m_nHyperedges;

	//! The Index that will be assigned to the next created hypernode.
	int m_hypernodeIdCount;

	//! The Index that will be assigned to the next created hyperedge.
	int m_hyperedgeIdCount;

public:
	//! Constructs an empty hypergraph.
	Hypergraph();

	//! Constructs a hypergraph that is a copy of \p H.
	Hypergraph(const Hypergraph& H);

	//! Destructor.
	~Hypergraph();

	//! Returns true iff the hypergraph is empty (ie. contains no hypernodes).
	bool empty() const { return m_nHypernodes == 0; }

	//! Returns the list of all hypernodes.
	const internal::GraphList<HypernodeElement>& hypernodes() const { return m_hypernodes; }

	//! Returns the list of all hyperedges.
	const internal::GraphList<HyperedgeElement>& hyperedges() const { return m_hyperedges; }

	//! Returns the number of hypernodes in the hypergraph.
	int numberOfHypernodes() const { return m_nHypernodes; }

	//! Returns the number of hyperedges in the hypergraph.
	int numberOfHyperedges() const { return m_nHyperedges; }

	//! Returns the largest used hypernode index.
	int maxHypernodeIndex() const { return m_hypernodeIdCount - 1; }

	//! Returns the largest used hyperedge index.
	int maxHyperedgeIndex() const { return m_hyperedgeIdCount - 1; }

	//! Returns the first hypernode in the list of all hypernodes.
	hypernode firstHypernode() const { return m_hypernodes.head(); }

	//! Returns the last hypernode in the list of all hypernodes.
	hypernode lastHypernode() const { return m_hypernodes.tail(); }

	//! Returns the first hyperedge in the list of all hyperedges.
	hyperedge firstHyperedge() const { return m_hyperedges.head(); }

	//! Returns the last hyperedge in the list of all hyperedges.
	hyperedge lastHyperEdge() const { return m_hyperedges.tail(); }

	//! Creates a new hypernode and returns it.
	hypernode newHypernode();

	//! Creates a new hypernode with given \p pIndex and returns it.
	hypernode newHypernode(int pIndex);

	//! Creates a new hypernode with given \p pType and returns it.
	hypernode newHypernode(HypernodeElement::Type pType);

	//! Creates a new hypernode with given \p pIndex and \p pType and returns it.
	hypernode newHypernode(int pIndex, HypernodeElement::Type pType);

	//! Creates a new hyperedge btween \p hypernodes and returns it.
	/**
	 * @param hypernodes are the hypernodes of the newly created hyperedge.
	 * @return the newly created hyperedge.
	 */
	hyperedge newHyperedge(List<hypernode>& hypernodes);

	//! Creates a new hyperedge between \p hypernodes and returns it.
	/**
	 * @param pIndex is the unique index of the newly created hyperedge.
	 * @param hypernodes are the hypernodes of the newly created hyperedge.
	 * @return the newly created hyperedge.
	 */
	hyperedge newHyperedge(int pIndex, List<hypernode>& hypernodes);

	//! Removes hypernode \p v and all incident hyperedges from the hypergraph.
	/**
	 * @param v is the hypernode that will be deleted.
	 */
	void delHypernode(hypernode v);

	//! Removes hyperedge \p e from the hypergraph.
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
	template<class LIST>
	void allHypernodes(LIST& hypernodeList) const {
		hypernodeList.clear();
		for (hypernode v = m_hypernodes.head(); v; v = v->succ()) {
			hypernodeList.pushBack(v);
		}
	}

	//! Returns a list with all hyperedges of the hypergraph.
	template<class LIST>
	void allHyperedges(LIST& hyperedgeList) const {
		hyperedgeList.clear();
		for (hyperedge e = m_hyperedges.head(); e; e = e->succ()) {
			hyperedgeList.pushBack(e);
		}
	}

	//! Reads hypergraph in bench format from the input stream.
	void readBenchHypergraph(std::istream& is);

	//! Reads hypergraph in bench format from the file.
	void readBenchHypergraph(const char* filename);

	//! Reads hypergraph in pla format from the input stream.
	void readPlaHypergraph(std::istream& is);

	//! Reads hypergraph in pla format from the file.
	void loadPlaHypergraph(const char* fileName);

	//! Checks the consistency of the data structure.
	bool consistency() const;

	//! Returns a reference to the registry of hypernode arrays associated with this hypergraph.
	HypergraphRegistry<HypernodeElement>& hypernodeRegistry() { return m_regHypernodeArrays; }

	//! Returns a const reference to the registry of hypernode arrays associated with this hypergraph.
	const HypergraphRegistry<HypernodeElement>& hypernodeRegistry() const {
		return m_regHypernodeArrays;
	}

	operator const HypergraphRegistry<HypernodeElement>&() const { return m_regHypernodeArrays; }

	//! Returns a reference to the registry of hyperedge arrays associated with this hypergraph.
	HypergraphRegistry<HyperedgeElement>& hyperedgeRegistry() { return m_regHyperedgeArrays; }

	//! Returns a const reference to the registry of hyperedge arrays associated with this hypergraph.
	const HypergraphRegistry<HyperedgeElement>& hyperedgeRegistry() const {
		return m_regHyperedgeArrays;
	}

	operator const HypergraphRegistry<HyperedgeElement>&() const { return m_regHyperedgeArrays; }

	Hypergraph& operator=(const Hypergraph& H);

	friend std::ostream& operator<<(std::ostream& os, ogdf::Hypergraph& H);

	friend std::istream& operator>>(std::istream& is, ogdf::Hypergraph& H);

	OGDF_MALLOC_NEW_DELETE;

private:
	void initArrays();

	int nextEntry(char* buffer, int from, string stop);

	HypernodeElement::Type gateType(string gate);
};

}
