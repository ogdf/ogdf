/*
 * $Revision: 3851 $
 *
 * last checkin:
 *   $Author: klein $
 *   $Date: 2013-11-20 07:26:56 +0100 (Mi, 20. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Declarations for Comparer objects.
 *
 * \author Markus Chimani, Carsten Gutwenger
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

#ifndef OGDF_COMPARER_H
#define OGDF_COMPARER_H

#include <ogdf/basic/basic.h>

namespace ogdf {

//--------------------------------------------------------------------
// A comparer interface has to define
// bool less (const E &x, const E &y);
// bool leq  (const E &x, const E &y);
// bool equal(const E &x, const E &y);
// bool geq (const E &x, const E &y);
// bool greater  (const E &x, const E &y);
//
// "const E &" can be replaced by "E"
//--------------------------------------------------------------------

//! Standard comparer (valid as a static comparer).
/**
 * Standard comparers are used by some sorting and searching methods.
 * The implementation of the generic class only provides dummies that
 * always throw a NoStdComparerException.
 *
 * The compare operations are static, hence the StdComparer cannot
 * only be used as a comparer object, but also as a static comparer
 * when required.
 *
 * You need to specialize this class for types you want to use with
 * sorting and searching methods like quicksort and binary search. There
 * already exist specializations for several standard types. If your type
 * already provides compare operators, you can use the macro #OGDF_STD_COMPARER
 * to automatically generate the specialization based on these operators.
 */
template<typename E> class StdComparer
{
public:
	static bool less(const E &x, const E &y) { OGDF_THROW(NoStdComparerException); }
	static bool leq(const E &x, const E &y) { OGDF_THROW(NoStdComparerException); }
	static bool greater(const E &x, const E &y) { OGDF_THROW(NoStdComparerException); }
	static bool geq(const E &x, const E &y) { OGDF_THROW(NoStdComparerException); }
	static bool equal(const E &x, const E &y) { OGDF_THROW(NoStdComparerException); }
};

//! Generates a specialization of the standard static comparer for \a type based on compare operators.
#define OGDF_STD_COMPARER(type) \
	template<> class StdComparer<type> \
	{ \
	public: \
		static bool less   (const type &x, const type &y) { return x <  y; } \
		static bool leq    (const type &x, const type &y) { return x <= y; } \
		static bool greater(const type &x, const type &y) { return x >  y; } \
		static bool geq    (const type &x, const type &y) { return x >= y; } \
		static bool equal  (const type &x, const type &y) { return x == y; } \
	};

OGDF_STD_COMPARER(short)
OGDF_STD_COMPARER(int)
OGDF_STD_COMPARER(float)
OGDF_STD_COMPARER(double)

//! Generates a specialization of the standard static comparer for booleans.
template<> class StdComparer<bool> {
public:
	static bool less   (const bool &x, const bool &y) { return !x &&  y; }
	static bool leq    (const bool &x, const bool &y) { return !x ||  y; }
	static bool greater(const bool &x, const bool &y) { return  x && !y; }
	static bool geq    (const bool &x, const bool &y) { return  x || !y; }
	static bool equal  (const bool &x, const bool &y) { return  x ==  y; }
};


//! A static comparer which compares the target of pointers ("content"), instead of the pointer's adresses.
/**
 * For the comparison of the contents, you may give your own static comparer
 */
template<class CONTENTTYPE, class STATICCONTENTCOMPARER = StdComparer<CONTENTTYPE> >
class TargetComparer {
	typedef CONTENTTYPE* CONTENTPOINTER;
public:
	static bool less   (const CONTENTPOINTER &x, const CONTENTPOINTER &y) { return STATICCONTENTCOMPARER::less   (*x,*y); }
	static bool leq    (const CONTENTPOINTER &x, const CONTENTPOINTER &y) { return STATICCONTENTCOMPARER::leq    (*x,*y); }
	static bool greater(const CONTENTPOINTER &x, const CONTENTPOINTER &y) { return STATICCONTENTCOMPARER::greater(*x,*y); }
	static bool geq    (const CONTENTPOINTER &x, const CONTENTPOINTER &y) { return STATICCONTENTCOMPARER::geq    (*x,*y); }
	static bool equal  (const CONTENTPOINTER &x, const CONTENTPOINTER &y) { return STATICCONTENTCOMPARER::equal  (*x,*y); }
};


//! Add this macro to your class to turn it into a full comparer.
/**
 * It is assumed that your class has a method "compare(const type &x, const type &y)", which
 * returns 0 if the two elements are equal, a negative value if \a x is smaller, and a positive
 * value if \a x is greater.
 *
 * Note: If the compare function of your class requires no additional data other than the
 * two elements to compare, your should usually use the more general #OGDF_AUGMENT_STATICCOMPARER:
 * A static comparer is also always valid as a normal comparer.
 *
 * Usage in Definition:
 * \code
 * class MyComparer {
 * private:
 *   Oracle oracle;
 * public:
 *   MyComparer(Oracle o) : oracle(o) {}
 *   int compare(const MyStuff& x1, const MyStuff& x2) const {
 *      return ... //compare x1 with x2, using oracle
 *   }
 *   OGDF_AUGMENT_COMPARER(MyStuff)
 * }
 * \endcode
 *
 * Use the Comparer:
 * \code
 * MyStuff a=...;
 * MyStuff b=...;
 * Oracle or;
 * MyComparer comp(or);
 * if( comp.less(a,b) )
 *   ... // do something
 * ...
 * Array<MyStuff> ay(10);
 * ... // fill array
 * ay.quicksort(comp); // sort the array using the MyComparer comp
 * \endcode
 */
#define OGDF_AUGMENT_COMPARER(type) \
	public: \
	bool less(const type &x, const type &y) const { return compare(x,y) < 0; } \
	bool leq(const type &x, const type &y) const { return compare(x,y) <= 0; } \
	bool greater(const type &x, const type &y) const { return compare(x,y) > 0; } \
	bool geq(const type &x, const type &y) const { return compare(x,y) >= 0; } \
	bool equal(const type &x, const type &y) const { return compare(x,y) == 0; }

//! Add this macro to your class to turn it into a full static comparer.
/**
 * It is assumed that your class has a *static* method "compare(const type &x, const type &y)", which
 * returns 0 if the two elements are equal, a negative value if \a x is smaller, and a positive
 * value if \a x is greater.
 *
 * Note: You should use this macro instead of #OGDF_AUGMENT_COMPARER whenever your compare function
 * requires no additional data stored in the object, other than the two elements to compare.
 * A static comparer is also always valid as a normal comparer.
 *
 * Usage in Definition:
 * \code
 * class MyComparer {
 * public:
 *    static int comparer(const MyStuff& x1, const MyStuff& x2) {
 *       return ... //compare x1 with x2
 *    }
 *    OGDF_AUGMENT_STATICCOMPARER(MyStuff)
 * }
 * \endcode
 *
 * Use the Comparer:
 * \code
 * MyStuff a=...;
 * MyStuff b=...;
 * MyComparer comp;
 * if( MyComparer.less(a,b) ) // use it statically on the class
 *    ... // do something
 * if( comp.less(a,b) ) // use it on the object
 *    ... // do something
 * ...
 * Array<MyStuff> ay(10);
 * ... // fill array
 * ay.quicksort(comp); // sort the array using the MyComparer comp
 * \endcode
 */
#define OGDF_AUGMENT_STATICCOMPARER(type) \
	public: \
	static bool less(const type &x, const type &y) { return compare(x,y) < 0; } \
	static bool leq(const type &x, const type &y) { return compare(x,y) <= 0; } \
	static bool greater(const type &x, const type &y) { return compare(x,y) > 0; } \
	static bool geq(const type &x, const type &y) { return compare(x,y) >= 0; } \
	static bool equal(const type &x, const type &y) { return compare(x,y) == 0; }


//! Abstract base class for comparer classes.
/**
 * The parameterized class \a VComparer<E> is an abstract base class for
 * encapsulating compare functions for type \a E. Implementations derive
 * from this class and implement at least the compare() method.
 *
 * The methods of this class are all \a virtual, which comes with a
 * certain performance penalty. Its advantage is that if you require
 * multiple Comparers for the same class \a E, functions using
 * compareres on \a E are not generated multiple times, which means
 * smaller code.
 *
 * If size is not an issue, but speed is, use a Comparer with
 * non-virtual functions. You may want to use the convenience classes
 * StdComparer and TargetComparer, or the convenience macros
 * #OGDF_AUGMENT_COMPARER, #OGDF_AUGMENT_STATICCOMPARER, #OGDF_STD_COMPARER to
 * obtain non-virtual classes with few effort.
 */
template<class E> class VComparer {
public:
	//! Initializes a comparer.
	VComparer() { }

	virtual ~VComparer() { }

	//! Compares \a x and \a y and returns the result as an integer.
	/** The returns value is
	 *  - < 0 iff x < y,
	 *  - = 0 iff x = y,
	 *  - > 0 iff x > y
	 */
	virtual int compare(const E &x, const E &y) const = 0;

	//! Returns true iff \a x < \a y
	virtual bool less(const E &x, const E &y) const {
		return compare(x,y) < 0;
	}

	//! Returns true iff \a x <= \a y
	virtual bool leq(const E &x, const E &y) const {
		return compare(x,y) <= 0;
	}

	//! Returns true iff \a x > \a y
	virtual bool greater(const E &x, const E &y) const {
		return compare(x,y) > 0;
	}

	//! Returns true iff \a x >= \a y
	virtual bool geq(const E &x, const E &y) const {
		return compare(x,y) >= 0;
	}

	//! Returns true iff \a x = \a y
	virtual bool equal(const E &x, const E &y) const {
		return compare(x,y) == 0;
	}
}; // class VComparer


//! Augments any data elements of type \a X with keys of type \a Priority. This class is also its own Comparer
/**
 * Also defines comparator function using the keys.
 * This class is intended as a helpful convenience class for using with BinaryHeapSimple, Top10Heap,..
 */
 template<class X, class Priority=double>
 class Prioritized {

	 X x;
	 Priority p;

 public:
	 //! Constructor of empty element. Be careful!
	 Prioritized() : x(0), p(0) { }

	 //! Constructor using a key/value pair
	 Prioritized(X xt, Priority pt) : x(xt),p(pt) { }

	 //! Copy-constructor
	 Prioritized(const Prioritized& P) : x(P.x),p(P.p) { }

	 //! Returns the key of the element
	 Priority priority() const { return p; }

	 //! Returns the data of the element
	 X item() const { return x;}

	 //! Sets priority
	 void setPriority(Priority pp) { p = pp; }

	 //! Sets value x
	 void setItem(X item) { x=item; }

	 //! Comparison oprator based on the compare-operator for the key type (\a Priority)
	 bool operator<(const Prioritized<X,Priority>& P) const { return p<P.p; }

	 //! Comparison oprator based on the compare-operator for the key type (\a Priority)
	 bool operator<=(const Prioritized<X,Priority>& P) const { return p<=P.p; }

	 //! Comparison oprator based on the compare-operator for the key type (\a Priority)
	 bool operator>(const Prioritized<X,Priority>& P) const { return p>P.p; }

	 //! Comparison oprator based on the compare-operator for the key type (\a Priority)
	 bool operator>=(const Prioritized<X,Priority>& P) const { return p>=P.p; }

	 //! Comparison oprator based on the compare-operator for the key type (\a Priority)
	 bool operator==(const Prioritized<X,Priority>& P) const { return p==P.p; }

	 //! Comparison oprator based on the compare-operator for the key type (\a Priority)
	 bool operator!=(const Prioritized<X,Priority>& P) const { return p!=P.p; }
 };

template<class X, class Priority> class StdComparer< Prioritized<X,Priority> >
{
public:
	static bool less   (const Prioritized<X,Priority> &x, const Prioritized<X,Priority> &y) { return x <  y; }
	static bool leq    (const Prioritized<X,Priority> &x, const Prioritized<X,Priority> &y) { return x <= y; }
	static bool greater(const Prioritized<X,Priority> &x, const Prioritized<X,Priority> &y) { return x >  y; }
	static bool geq    (const Prioritized<X,Priority> &x, const Prioritized<X,Priority> &y) { return x >= y; }
	static bool equal  (const Prioritized<X,Priority> &x, const Prioritized<X,Priority> &y) { return x == y; }
};


} //namespace

#endif /*OGF_COMPARER_H*/
