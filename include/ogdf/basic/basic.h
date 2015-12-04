/** \file
 * \brief Basic declarations, included by all source files.
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

#pragma once

#include <ogdf/internal/basic/config.h>


//---------------------------------------------------------
// assertions
//---------------------------------------------------------

#ifdef OGDF_DEBUG
#include <cassert>
#define OGDF_ASSERT(expr) assert(expr);
#define OGDF_ASSERT_IF(minLevel,expr) \
	if (int(ogdf::debugLevel) >= int(minLevel)) { assert(expr); } else { }
#define OGDF_SET_DEBUG_LEVEL(level) ogdf::debugLevel = level;

#define OGDF_DEBUG_OUTPUT(kind, str) \
	std::cerr << std::endl << "OGDF " << kind << " from " << __FILE__ << ":" << __LINE__ << "(" << __func__ << "): " << str << std::endl;

#else
#define OGDF_ASSERT(expr)
#define OGDF_ASSERT_IF(minLevel,expr)
#define OGDF_SET_DEBUG_LEVEL(level)
#define OGDF_DEBUG_OUTPUT(kind, str)
#endif

//---------------------------------------------------------
// deprecation
//---------------------------------------------------------


#if __cplusplus >= 201402L
#define OGDF_DEPRECATED [[deprecated]]
#else

#ifdef _MSC_VER
#define OGDF_DEPRECATED __declspec(deprecated)
#elif defined(__GNUC__)
#define OGDF_DEPRECATED __attribute__ ((deprecated))
#else
#define OGDF_DEPRECATED
#endif

#endif



//---------------------------------------------------------
// macros for optimization
//---------------------------------------------------------

// Visual C++ compiler
#ifdef _MSC_VER

#define OGDF_LIKELY(x)    (x)
#define OGDF_UNLIKELY(x)  (x)

#define OGDF_DECL_ALIGN(b) __declspec(align(b))
#define OGDF_DECL_THREAD __declspec(thread)


// GNU gcc compiler (also Intel compiler)
#elif defined(__GNUC__)
//// make sure that SIZE_MAX gets defined
//#define __STDC_LIMIT_MACROS

#define OGDF_LIKELY(x)    __builtin_expect((x),1)
#define OGDF_UNLIKELY(x)  __builtin_expect((x),0)

#define OGDF_DECL_ALIGN(b) __attribute__ ((aligned(b)))
#define OGDF_DECL_THREAD __thread


// other compiler
#else
#define OGDF_LIKELY(x)    (x)
#define OGDF_UNLIKELY(x)  (x)

#define OGDF_DECL_ALIGN(b)
#endif

#ifndef __SIZEOF_POINTER__
#ifdef _M_X64
#define __SIZEOF_POINTER__ 8
#else
#define __SIZEOF_POINTER__ 4
#endif
#endif


//---------------------------------------------------------
// common includes
//---------------------------------------------------------

// stdlib
#include <cstdint>
#include <cmath>
#include <ctime>
#include <fstream>
#include <algorithm>
#include <limits>

//! The namespace for all OGDF objects.
namespace ogdf {

using std::ifstream;		// from <fstream>
using std::ofstream;		// from <fstream>
using std::min;				// from <algorithm>
using std::max;				// from <algorithm>
using std::numeric_limits;	// from <limits>

#ifndef OGDF_DLL

/**
 *  The class Initialization is used for initializing global variables.
 *  You should never create instances of it!
*/
class Initialization {
	static int s_count;

public:
	Initialization();
	~Initialization();
};

static Initialization s_ogdfInitializer;

#endif


	// forward declarations
	template<class E> class List;



	enum Direction { before, after };

	/**
	 * @addtogroup random
	 */
	//@{

	//! Returns a random value suitable as initial seed for a random number engine.
	/**
	 * <H3>Thread Safety</H3>
	 * This functions is thread-safe.
	 */
	OGDF_EXPORT long unsigned int randomSeed();

	//! Sets the seed for functions like randomSeed(), randomNumber(), randomDouble().
	OGDF_EXPORT void setSeed(int val);

	//! Returns random integer between low and high (including).
	/**
	 * <H3>Thread Safety</H3>
	 * This functions is thread-safe.
	 */
	OGDF_EXPORT int randomNumber(int low, int high);

	//! Returns a random double value from the interval [\a low, \a high).
	/**
	 * <H3>Thread Safety</H3>
	 * This functions is thread-safe.
	 */
	OGDF_EXPORT double randomDouble(double low, double high);

	//! Returns a random double value from the normal distribution
	//! with mean m and standard deviation sd
	inline double randomDoubleNormal(double m, double sd)
	{
		double x1, y1, w;

		do {
			double rndVal = randomDouble(0,1);
			x1 = 2.0 * rndVal - 1.0;
			rndVal = randomDouble(0,1);
			double x2 = 2.0 * rndVal - 1.0;
			w = x1*x1 + x2*x2;
		} while (w >= 1.0);

		w = sqrt((-2.0 * log(w))/w) ;
		y1 = x1*w;

		return(m + y1*sd);
	}

	//@}


	//! Returns used CPU time from T to current time and assigns current time to T.
	/**
	 * @ingroup date-time
	 */
	OGDF_EXPORT double usedTime(double& T);

	//! \a doDestruction() returns false if a data type does not require to
	//! call its destructor (e.g. build-in data types).
	template<class E>inline bool doDestruction(const E *) { return true; }

	// specializations
	template<>inline bool doDestruction(const char *) { return false; }
	template<>inline bool doDestruction<int>(const int *) { return false; }
	template<>inline bool doDestruction<double>(const double *) { return false; }


	//! Removes trailing space, horizontal and vertical tab, feed, newline, and carriage return  from \a str.
	OGDF_EXPORT void removeTrailingWhitespace(string &str);

	//! Compares the two strings \a str1 and \a str2, ignoring the case of characters.
	OGDF_EXPORT bool equalIgnoreCase(const string &str1, const string &str2);

	//! Tests if \a prefix is a prefix of \a str, ignoring the case of characters.
	OGDF_EXPORT bool prefixIgnoreCase(const string &prefix, const string &str);

	//@}
	/**
	* @addtogroup container-functions
	*/
	//@{

	//! Searches for the position of \a x in container \a C; returns -1 if not found.
	/**
	 * Positions are number 0, 1, 2, ... The function uses the equality operator for comparing elements.
	 *
	 * \param C is a container containing elements of type \a T.
	 * \param x is the element to search for.
	 * \return the position of the first occurrence of \a x in \a C (positions start with 0), or -1 if
	 *         \a x is not in \a C.
	 */
	template< typename CONTAINER, typename T>
	int searchPos(const CONTAINER &C, const T &x)
	{
		int pos = 0;
		for (const T &y : C) {
			if (x == y)
				return pos;
			++pos;
		}

		return -1;
	}

	//@}


#ifdef OGDF_DEBUG
	/** We maintain a debug level in debug versions indicating how many
	 *  internal checks (usually assertions) are done.
	 *  Usage: Set the variable ogdf::debugLevel using the macro
	 *   OGDF_SET_DEBUG_LEVEL(level) to the desired level
	 *   in the calling code (e.g. main()). The debugLevel can be set
	 *   to a higher level for critical parts (e.g., where you assume a bug)
	 *   ensuring that other parts are not too slow.
	 */
	enum DebugLevel {
		dlMinimal, dlExtendedChecking, dlConsistencyChecks, dlHeavyChecks
	};
	extern DebugLevel debugLevel;
#endif


//! Abstract base class for bucket functions.
/**
 * The parameterized class \a BucketFunc<E> is an abstract base class
 * for bucket functions. Derived classes have to implement \a getBucket().
 * Bucket functions are used by bucket sort functions for container types.
 */
template<class E> class BucketFunc
{
public:
	virtual ~BucketFunc() { }

	//! Returns the bucket of \a x.
	virtual int getBucket(const E &x) = 0;
};

using std::stoi;
using std::stoll;
using std::stoul;
using std::stoull;

using std::stof;
using std::stod;
using std::stold;

using std::to_string;

} // end namespace ogdf


/**
 * @name C++11 support
 * OGDF uses some new functions defined in the C++11 standard; if these functions are not
 * supported by the compiler, they are defined by OGDF.
 */
//@{


//@}
