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

#include <ogdf/internal/basic/config.h>

//! @name Assertions (only active in debug builds)
//! @{

//! Assert condition \a expr. See doc/build.md for more information.
//! @ingroup macros
#define OGDF_ASSERT(expr)
//! Assert condition \a expr if debug level is at least \a minLevel.
//! @ingroup macros
#define OGDF_ASSERT_IF(minLevel,expr)
//! Set debug level to \a level.
//! @ingroup macros
#define OGDF_SET_DEBUG_LEVEL(level)

#ifdef OGDF_DEBUG
# undef OGDF_ASSERT
# ifndef OGDF_USE_ASSERT_EXCEPTIONS
#  include <cassert>
#  define OGDF_ASSERT(expr) assert(expr)
# else
#  include <stdexcept>
#  include <sstream>

namespace ogdf {
/**
 * A trivial exception for failed assertions.
 * Only available if the macro OGDF_USE_ASSERT_EXCEPTIONS is defined.
 */
class AssertionFailed : public std::runtime_error {
	using std::runtime_error::runtime_error;
};
}

#  define OGDF_ASSERT(expr) do { \
	if (!(expr)) { \
		std::stringstream ogdf_assert_ss; \
		ogdf_assert_ss \
		 << "OGDF assertion `" #expr "' failed at " __FILE__ ":" \
		 << __LINE__ \
		 << "(" << OGDF_FUNCTION_NAME << ")"; \
		ogdf::get_stacktrace(ogdf_assert_ss); \
		throw ogdf::AssertionFailed(ogdf_assert_ss.str()); \
	} } while (false)
# endif
# undef OGDF_ASSERT_IF
# define OGDF_ASSERT_IF(minLevel,expr) do { \
	if (int(ogdf::debugLevel) >= int(minLevel)) { \
		OGDF_ASSERT(expr); \
	} } while (false)
# undef OGDF_SET_DEBUG_LEVEL
# define OGDF_SET_DEBUG_LEVEL(level) ogdf::debugLevel = level
#endif

//! @}

//---------------------------------------------------------
// common includes
//---------------------------------------------------------

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
public:
	Initialization();
	~Initialization();
};

// This has to be in the header file. Being in the cpp file does not
// guarantee that it is constructed (with all linkers).
static Initialization s_ogdfInitializer;

#endif

#ifdef OGDF_USE_ASSERT_EXCEPTIONS
	//! Output a mangled stack backtrace of the caller function to stream
	extern void get_stacktrace(std::ostream &);
#endif

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

		return m + y1 * sd;
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
