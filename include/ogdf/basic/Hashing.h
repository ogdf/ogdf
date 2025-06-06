/** \file
 * \brief Declaration of classes used for hashing.
 *
 * Declares HashingBase and HashElementBase, and declares and implements
 * classes Hashing, HashElement, HashConstIterator.
 *
 * \author Carsten Gutwenger
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

#include <ogdf/basic/basic.h>
#include <ogdf/basic/memory.h>

#include <cmath>
#include <cstddef>
#include <limits>
#include <string>

namespace ogdf {

//! Base class for elements within a hash table.
/**
 * This class realizes only chaining of elements and maintianing hash values
 * for rehashing.
 */
class OGDF_EXPORT HashElementBase {
	friend class HashingBase;

	HashElementBase* m_next; //!< The successor in the list.
	size_t m_hashValue; //!< The hash value.

public:
	//! Creates a hash element with hash value \p hashValue.
	explicit HashElementBase(size_t hashValue) : m_next(nullptr), m_hashValue(hashValue) { }

	//! Returns the successor to this element in the list.
	HashElementBase* next() const { return m_next; }

	//! Returns the hash value of this element.
	size_t hashValue() const { return m_hashValue; }

	OGDF_NEW_DELETE
};

/**
 * \brief Base class for hashing with chaining and table doubling.
 *
 * The actual hashing is provided by the parameterized class Hashing<K,I>
 * which derives from HashingBase.
 */
class OGDF_EXPORT HashingBase {
protected:
	int m_tableSize; //!< The current table size.
	int m_hashMask; //!< The current table size minus one.
	int m_minTableSize; //!< The minimal table size.
	int m_tableSizeLow; //!< The minimal number of elements at this table size.
	int m_tableSizeHigh; //!< The maximal number of elements at this table size.
	int m_count; //!< The current number of elements.
	HashElementBase** m_table; //!< The hash table (an array of lists).

public:
	//! Creates a hash table with minimum table size \p minTableSize.
	explicit HashingBase(int minTableSize);

	//! Copy constructor.
	HashingBase(const HashingBase& H);

	//! Destruction
	virtual ~HashingBase();

	//! Resizes the hash table to \p newTableSize.
	void resize(int newTableSize);

	//! Inserts a new element \p pElement into the hash table.
	void insert(HashElementBase* pElement);

	//! Removes the element \p pElement from the hash table.
	void del(HashElementBase* pElement);

	//! Removes all elements from the hash table.
	void clear();

	//! Assignment operator.
	HashingBase& operator=(const HashingBase& H);

	//! Returns the number of elements in the hash table.
	int size() const { return m_count; }

	//! Returns if the hash table is empty
	int empty() const { return m_count == 0; }

	/**
	 * \brief Returns the first element in the list for elements with hash value \p hashValue.
	 *
	 * This is the list #m_table[\ hashValue & #m_hashMask].
	 */
	HashElementBase* firstListElement(size_t hashValue) const {
		return *(m_table + (hashValue & m_hashMask));
	}

	/**
	 * \brief Returns the first element in the list of all elements in the hash table.
	 *
	 * This function is used by hash iterators for iterating over all elements
	 * in the hash table.
	 * @param pList is assigned the list containing the first element.
	 * \return a pointer to the first element or \c nullptr if hash table is empty.
	 */
	HashElementBase* firstElement(HashElementBase*** pList) const;

	/**
	 * \brief Returns the successor of \p pElement in the list of all elements in the hash table.
	 *
	 * This function is used by hash iterators for iterating over all elements
	 * in the hash table.
	 * @param pList is assigned the list containing the first element.
	 * @param pElement points to an element in the has table.
	 * \return a pointer to the first element or \c nullptr if hash table is empty.
	 */
	HashElementBase* nextElement(HashElementBase*** pList, HashElementBase* pElement) const;

protected:
	//! Deletes all elements in hash table (but does not free m_table!).
	void destroyAll();

	/**
	 * \brief Called to delete hash element.
	 *
	 * This must be done in Hashing<K,I> since only this class knows the actual
	 * element type; alternatively, HashElementBase could have a virtual destructor.
	 */
	virtual void destroy(HashElementBase* pElement) = 0;

	//! Called to create a copy of the element \p pElement.
	virtual HashElementBase* copy(HashElementBase* pElement) const = 0;

private:
	//! Initializes the table for given table size.
	void init(int tableSize);

	//! Copies all elements from \p H to this hash table.
	void copyAll(const HashingBase& H);
};

/**
 * \brief Representation of elements in a hash table.
 *
 * This class adds key and information members to HashElementBase. The two
 * template parameters are \a K for the type of keys and \p I for the type
 * of information.
 */
template<class K, class I>
class HashElement : public HashElementBase {
	K m_key; //!< The key value.
	I m_info; //!< The information value.

public:
	//! Creates a hash element with given hash value, key, and information.
	HashElement(size_t hashValue, const K& key, const I& info)
		: HashElementBase(hashValue), m_key(key), m_info(info) { }

	//! Returns the successor element in the list.
	HashElement<K, I>* next() const { return (HashElement<K, I>*)HashElementBase::next(); }

	//! Returns the key value.
	const K& key() const { return m_key; }

	//! Returns the information value.
	const I& info() const { return m_info; }

	//! Returns a refeence to the information value.
	I& info() { return m_info; }

	OGDF_NEW_DELETE
};


template<class K, class I, class H>
class HashConstIterator;

/**
 * \brief Default hash functions.
 *
 * This class implements a default hash function for various
 * basic data types.
 *
 * \see Hashing, HashArray, HashArray2D
 */
template<class K>
class DefHashFunc {
public:
	//! Returns the hash value of \p key.
	size_t hash(const K& key) const { return size_t(key); }
};

//! Specialized default hash function for pointer types.
template<>
class DefHashFunc<void*> {
public:
	size_t hash(const void*& key) const { return size_t(key && 0xffffffff); }
};

//! Specialized default hash function for double.
template<>
class DefHashFunc<double> {
public:
	size_t hash(const double& key) const {
		int dummy;
		return (size_t)((double)std::numeric_limits<size_t>::max() * frexp(key, &dummy));
	}
};

//! Specialized default hash function for string.
template<>
class DefHashFunc<string> {
public:
	size_t hash(const string& key) const;
};

//! %Hashing with chaining and table doubling.
/**
 * @ingroup containers
 *
 * The class Hashing<K,I> implements a hashing table which realizes a
 * mapping from a key type \a K to an information type \a I.
 *
 * The class requires three template parameters:
 * @tparam K is the type of keys.
 * @tparam I is the type of information.
 * @tparam H is the hash function type.
 * The hash function type argument is optional; its default uses the class
 * DefHashFunc.
 *
 * Hash function classes have to define `int hash(const E &key)`
 * or `int hash(E key)`.
 */
template<class K, class I, class H = DefHashFunc<K>>
class Hashing : private HashingBase {
	friend class HashConstIterator<K, I, H>;
	H m_hashFunc; //!< The hash function.

public:
	//! The type of const-iterators for hash tables.
	using const_iterator = HashConstIterator<K, I, H>;

	//! Creates a hash table for given initial table size \p minTableSize.
	explicit Hashing(int minTableSize = 256, const H& hashFunc = H())
		: HashingBase(minTableSize), m_hashFunc(hashFunc) { }

	//! Copy constructor.
	Hashing(const Hashing<K, I, H>& h) = default;

	//! Destruction
	~Hashing() { HashingBase::destroyAll(); }

	//! Returns the number of elements in the hash table.
	int size() const { return HashingBase::size(); }

	//! Returns true iff the table is empty, i.e., contains no elements.
	bool empty() const { return HashingBase::size() == 0; }

	//! Returns true iff the hash table contains an element with key \p key.
	bool member(const K& key) const { return lookup(key) != nullptr; }

	//! Returns an hash iterator to the first element in the list of all elements.
	HashConstIterator<K, I, H> begin() const;

	//! Returns the hash element with key \p key in the hash table; returns \c nullptr if no such element exists.
	HashElement<K, I>* lookup(const K& key) const {
		HashElement<K, I>* pElement = (HashElement<K, I>*)firstListElement(m_hashFunc.hash(key));
		for (; pElement; pElement = pElement->next()) {
			if (pElement->key() == key) {
				return pElement;
			}
		}

		return nullptr;
	}

	//! Assignment operator.
	Hashing<K, I, H>& operator=(const Hashing<K, I, H>& hashing) = default;

	/**
	 * \brief Inserts a new element with key \p key and information \p info into the hash table.
	 *
	 * The new element will only be inserted if no element with key \p key is
	 * already contained; if such an element already exists the information of
	 * this element will be changed to \p info.
	 */
	HashElement<K, I>* insert(const K& key, const I& info) {
		HashElement<K, I>* pElement = lookup(key);

		if (pElement) {
			pElement->info() = info;
		} else {
			HashingBase::insert(pElement = new HashElement<K, I>(m_hashFunc.hash(key), key, info));
		}

		return pElement;
	}

	/**
	 * \brief Inserts a new element with key \p key and information \p info into the hash table.
	 *
	 * The new element will only be inserted if no element with key \p key is
	 * already contained; if such an element already exists the information of
	 * this element remains unchanged.
	 */
	HashElement<K, I>* insertByNeed(const K& key, const I& info) {
		HashElement<K, I>* pElement = lookup(key);

		if (!pElement) {
			HashingBase::insert(pElement = new HashElement<K, I>(m_hashFunc.hash(key), key, info));
		}

		return pElement;
	}

	/**
	 * \brief Inserts a new element with key \p key and information \p info into the hash table.
	 *
	 * This is a faster version of insert() that assumes that no element with key
	 * \p key is already contained in the hash table.
	 */
	HashElement<K, I>* fastInsert(const K& key, const I& info) {
		HashElement<K, I>* pElement = new HashElement<K, I>(m_hashFunc.hash(key), key, info);
		HashingBase::insert(pElement);
		return pElement;
	}

	//! Removes the element with key \p key from the hash table (does nothing if no such element).
	void del(const K& key) {
		HashElement<K, I>* pElement = lookup(key);
		if (pElement) {
			HashingBase::del(pElement);
			delete pElement;
		}
	}

	//! Removes all elements from the hash table.
	void clear() { HashingBase::clear(); }

protected:
	/**
	 * \brief Returns the first element in the list of all elements in the hash table.
	 *
	 * This function is used by hash iterators for iterating over all elements
	 * in the hash table.
	 * @param pList is assigned the list containing the first element.
	 * \return a pointer to the first element or \c nullptr if hash table is empty.
	 */
	HashElement<K, I>* firstElement(HashElement<K, I>*** pList) const {
		return (HashElement<K, I>*)(HashingBase::firstElement((HashElementBase***)pList));
	}

	/**
	 * \brief Returns the successor of \p pElement in the list of all elements in the hash table.
	 *
	 * This function is used by hash iterators for iterating over all elements
	 * in the hash table.
	 * @param pList is assigned the list containing the first element.
	 * @param pElement points to an element in the has table.
	 * \return a pointer to the first element or \c nullptr if hash table is empty.
	 */
	HashElement<K, I>* nextElement(HashElement<K, I>*** pList, HashElement<K, I>* pElement) const {
		return (HashElement<K, I>*)(HashingBase::nextElement((HashElementBase***)pList, pElement));
	}

private:
	//! Deletes hash element \p pElement.
	virtual void destroy(HashElementBase* pElement) override {
		delete (HashElement<K, I>*)(pElement);
	}

	//! Returns a copy of hash element \p pElement.
	virtual HashElementBase* copy(HashElementBase* pElement) const override {
		HashElement<K, I>* pX = (HashElement<K, I>*)(pElement);
		return new HashElement<K, I>(pX->hashValue(), pX->key(), pX->info());
	}
};

/**
 * \brief Iterators for hash tables.
 *
 * This class implements an iterator for iterating over all elements in
 * a hash table. Hash iterators are provided by Hashing<K,I>::begin().
 *
 * <H3>Example</H3>
 * The following code snippet demonstrates how to iterate over all elements
 * of a hash table. First, the example fills a hash table with a tiny
 * German&ndash;English dictionary, and then it iterates over the elements
 * and prints the entries.
 * \code
 *   Hashing<string,string> H;
 *
 *   H.fastInsert("Hund","dog");
 *   H.fastInsert("Katze","cat");
 *   H.fastInsert("Maus","mouse");
 *
 *   HashConstIterator<string,string> it;
 *   for(it = H.begin(); it.valid(); ++it)
 *     std::cout << it.key() << " -> " << it.info() << std::endl;
 * \endcode
 */
template<class K, class I, class H = DefHashFunc<K>>
class HashConstIterator {
	HashElement<K, I>* m_pElement; //!< The hash element to which the iterator points.
	HashElement<K, I>** m_pList; //!< The list containg the hash element.
	const Hashing<K, I, H>* m_pHashing; //!< The associated hash table.

public:
	//! Creates a hash iterator pointing to no element.
	HashConstIterator() : m_pElement(nullptr), m_pList(nullptr), m_pHashing(nullptr) { }

	//! Creates a hash iterator pointing to element \p pElement in list \p pList of hash table \p pHashing.
	HashConstIterator(HashElement<K, I>* pElement, HashElement<K, I>** pList,
			const Hashing<K, I, H>* pHashing)
		: m_pElement(pElement), m_pList(pList), m_pHashing(pHashing) { }

	//! Copy constructor.
	HashConstIterator(const HashConstIterator<K, I, H>& it)
		: m_pElement(it.m_pElement), m_pList(it.m_pList), m_pHashing(it.m_pHashing) { }

	//! Assignment operator.
	HashConstIterator& operator=(const HashConstIterator<K, I, H>& it) {
		m_pElement = it.m_pElement;
		m_pList = it.m_pList;
		m_pHashing = it.m_pHashing;
		return *this;
	}

	//! Returns true if the hash iterator points to an element.
	bool valid() const { return m_pElement != nullptr; }

	//! Returns the key of the hash element pointed to.
	const K& key() const { return m_pElement->key(); }

	//! Returns the information of the hash element pointed to.
	const I& info() const { return m_pElement->info(); }

	//! Equality operator.
	friend bool operator==(const HashConstIterator<K, I, H>& it1,
			const HashConstIterator<K, I, H>& it2) {
		return it1.m_pElement == it2.m_pElement;
	}

	//! Inequality operator.
	friend bool operator!=(const HashConstIterator<K, I, H>& it1,
			const HashConstIterator<K, I, H>& it2) {
		return it1.m_pElement != it2.m_pElement;
	}

	//! Moves this hash iterator to the next element (iterator gets invalid if no more elements).
	HashConstIterator<K, I, H>& operator++() {
		m_pElement = m_pHashing->nextElement(&m_pList, m_pElement);
		return *this;
	}
};

template<class K, class I, class H>
inline HashConstIterator<K, I, H> Hashing<K, I, H>::begin() const {
	HashElement<K, I>*pElement, **pList;
	pElement = firstElement(&pList);
	return HashConstIterator<K, I, H>(pElement, pList, this);
}

}
