/*
 * $Revision: 3804 $
 *
 * last checkin:
 *   $Author: chimani $
 *   $Date: 2013-10-29 12:02:20 +0100 (Di, 29. Okt 2013) $
 ***************************************************************/

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


#ifdef _MSC_VER
#pragma once
#endif

#ifndef OGDF_BASIC_H
#define OGDF_BASIC_H


/**
 * \mainpage The Open Graph Drawing Framework
 *
 * \section sec_intro Introduction
 * The Open Graph Drawing Framework (OGDF) is a C++ library containing
 * implementations of various graph drawing algorithms. The library is self
 * contained; optionally, additional packages like LP-solvers are required
 * for some implementations.
 *
 * Here, you find the library's code documentation. For more general information
 * on OGDF see http://www.ogdf.net. There, you can also find further explanations,
 * how-tos, and example code.
 *
 * The OGDF project is a cooperation between
 * - [Chair of Algorithm Engineering](http://ls11-www.cs.uni-dortmund.de/), Faculty of Computer Science, TU Dortmund, Germany
 * - [Theoretical Computer Science](http://www.cs.uos.de/theoinf), Department of Mathematics/Computer Science, Uni Osnabrück, Germany
 * - [Chair of Prof. J&uuml;nger](http://www.informatik.uni-koeln.de/ls_juenger/), Department of Computer Science, University of Cologne, Germany
 * - [University of Sydney](http://sydney.edu.au/engineering/it/), Australia
 * - [oreas GmbH](http://www.oreas.com/), Cologne, Germany
 */


#include <ogdf/internal/basic/config.h>


// include windows.h on Windows systems
#if defined(OGDF_SYSTEM_WINDOWS) || defined(__CYGWIN__)
#define WIN32_EXTRA_LEAN
#define WIN32_LEAN_AND_MEAN
#undef NOMINMAX
#define NOMINMAX
#include <windows.h>
#endif


//---------------------------------------------------------
// assertions
//---------------------------------------------------------

#ifdef OGDF_DEBUG
#include <assert.h>
#define OGDF_ASSERT(expr) assert(expr);
#define OGDF_ASSERT_IF(minLevel,expr) \
	if (int(ogdf::debugLevel) >= int(minLevel)) { assert(expr); } else { }
#define OGDF_SET_DEBUG_LEVEL(level) ogdf::debugLevel = level;

#else
#define OGDF_ASSERT(expr)
#define OGDF_ASSERT_IF(minLevel,expr)
#define OGDF_SET_DEBUG_LEVEL(level)
#endif


//---------------------------------------------------------
// macros for optimization
//---------------------------------------------------------

// Visual C++ compiler
#ifdef _MSC_VER

#define OGDF_LIKELY(x)    (x)
#define OGDF_UNLIKELY(x)  (x)

#ifdef OGDF_DEBUG
#define OGDF_NODEFAULT    default: assert(0);
#else
#define OGDF_NODEFAULT    default: __assume(0);
#endif

#define OGDF_DECL_ALIGN(b) __declspec(align(b))
#define OGDF_DECL_THREAD __declspec(thread)


// GNU gcc compiler (also Intel compiler)
#elif defined(__GNUC__)
//// make sure that SIZE_MAX gets defined
//#define __STDC_LIMIT_MACROS
//#include <stdint.h>

#define OGDF_LIKELY(x)    __builtin_expect((x),1)
#define OGDF_UNLIKELY(x)  __builtin_expect((x),0)
#define OGDF_NODEFAULT    default: ;

#define OGDF_DECL_ALIGN(b) __attribute__ ((aligned(b)))
#define OGDF_DECL_THREAD __thread


// other compiler
#else
#define OGDF_LIKELY(x)    (x)
#define OGDF_UNLIKELY(x)  (x)
#define OGDF_NODEFAULT

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
// define data types with known size
//---------------------------------------------------------

#if defined(_MSC_VER)

typedef unsigned __int8  __uint8;
typedef unsigned __int16 __uint16;
typedef unsigned __int32 __uint32;
typedef unsigned __int64 __uint64;

#else

#undef __int8
#undef __int16
#undef __int32
#undef __int64

typedef signed char        __int8;
typedef short              __int16;
typedef int                __int32;
typedef long long          __int64;
typedef unsigned char      __uint8;
typedef unsigned short     __uint16;
typedef unsigned int       __uint32;
typedef unsigned long long __uint64;
#endif


//---------------------------------------------------------
// common includes
//---------------------------------------------------------

// stdlib
#include <cmath>
#include <ctime>
#include <fstream>
#include <algorithm>
#include <limits>

using std::ifstream;		// from <fstream>
using std::ofstream;		// from <fstream>
using std::min;				// from <algorithm>
using std::max;				// from <algorithm>
using std::numeric_limits;	// from <limits>

#ifdef OGDF_SYSTEM_UNIX
#include <stdint.h>
#endif
// make sure that SIZE_MAX gets defined
#ifndef SIZE_MAX
#define SIZE_MAX ((size_t)-1)
#endif


// ogdf
#include <ogdf/basic/Logger.h>
#include <ogdf/basic/exceptions.h>
#include <ogdf/basic/System.h>
#include <ogdf/basic/memory.h>
#include <ogdf/basic/comparer.h>



//! The namespace for all OGDF objects.
namespace ogdf {

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


	/**
	 * @name Global basic functions
	 */
	//@{

	// forward declarations
	template<class E> class List;


	enum Direction { before, after };

	//! Returns random integer between low and high (including).
	inline int randomNumber(int low, int high) {
#if RAND_MAX == 32767
		// We get only 15 random bits on some systems (Windows, Solaris)!
		int r1 = (rand() & ((1 << 16) - 1));
		int r2 = (rand() & ((1 << 16) - 1));
		int r = (r1 << 15) | r2;
#else
		int r = rand();
#endif
		return low + (r % (high-low+1));
	}

	//! Returns random double value between low and high.
	inline double randomDouble(double low, double high) {
		double val = low +(rand()*(high-low))/RAND_MAX;
		OGDF_ASSERT(val >= low && val <= high);
		return val;
	}

	//! Returns a random double value from the normal distribution
	//! with mean m and standard deviation sd
	inline double randomDoubleNormal(double m, double sd)
	{
		double x1, x2, y1, w, rndVal;

		do {
			rndVal = randomDouble(0,1);
			x1 = 2.0 * rndVal - 1.0;
			rndVal = randomDouble(0,1);
			x2 = 2.0 * rndVal - 1.0;
			w = x1*x1 + x2*x2;
		} while (w >= 1.0);

		w = sqrt((-2.0 * log(w))/w) ;
		y1 = x1*w;

		return(m + y1*sd);
	}



	//! Returns used CPU time from T to current time and assigns
	//! current time to T.
	OGDF_EXPORT double usedTime(double& T);

	//! \a doDestruction() returns false if a data type does not require to
	//! call its destructor (e.g. build-in data types).
	template<class E>inline bool doDestruction(const E *) { return true; }

	// specializations
	template<>inline bool doDestruction(const char *) { return false; }
	template<>inline bool doDestruction<int>(const int *) { return false; }
	template<>inline bool doDestruction<double>(const double *) { return false; }


	//! Compares the two strings \a str1 and \a str2, ignoring the case of characters.
	OGDF_EXPORT bool equalIgnoreCase(const string &str1, const string &str2);

	//! Tests if \a prefix is a prefix of \a str, ignoring the case of characters.
	OGDF_EXPORT bool prefixIgnoreCase(const string &prefix, const string &str);

	//@}


	/**
	 * @name Files and directories
	 */
	//@{

	//! The type of an entry in a directory.
	enum FileType {
		ftEntry,     /**< file or directory */
		ftFile,      /**< file */
		ftDirectory  /**< directory */
	};

	//! Returns true iff \a fileName is a regular file (not a directory).
	OGDF_EXPORT bool isFile(const char *fileName);

	//! Returns true iff \a fileName is a directory.
	OGDF_EXPORT bool isDirectory(const char *fileName);

	//! Changes current directory to \a dirName; returns true if successful.
	OGDF_EXPORT bool changeDir(const char *dirName);

	//! Returns in \a files the list of files in directory \a dirName.
	/** The optional argument \a pattern can be used to filter files.
	 *
	 *  \pre \a dirName is a directory
	 */
	OGDF_EXPORT void getFiles(const char *dirName,
		List<string> &files,
		const char *pattern = "*");

	//! Appends to \a files the list of files in directory \a dirName.
	/** The optional argument \a pattern can be used to filter files.
	 *
	 *  \pre \a dirName is a directory
	 */
	OGDF_EXPORT void getFilesAppend(const char *dirName,
		List<string> &files,
		const char *pattern = "*");


	//! Returns in \a subdirs the list of directories contained in directory \a dirName.
	/** The optional argument \a pattern can be used to filter files.
	 *
	 *  \pre \a dirName is a directory
	 */
	OGDF_EXPORT void getSubdirs(const char *dirName,
		List<string> &subdirs,
		const char *pattern = "*");

	//! Appends to \a subdirs the list of directories contained in directory \a dirName.
	/** The optional argument \a pattern can be used to filter files.
	 *
	 *  \pre \a dirName is a directory
	 */
	OGDF_EXPORT void getSubdirsAppend(const char *dirName,
		List<string> &subdirs,
		const char *pattern = "*");


	//! Returns in \a entries the list of all entries contained in directory \a dirName.
	/** Entries may be files or directories. The optional argument \a pattern
	 *  can be used to filter files.
	 *
	 *  \pre \a dirName is a directory
	 */
	OGDF_EXPORT void getEntries(const char *dirName,
		List<string> &entries,
		const char *pattern = "*");

	//! Appends to \a entries the list of all entries contained in directory \a dirName.
	/** Entries may be files or directories. The optional argument \a pattern
	 *  can be used to filter files.
	 *
	 *  \pre \a dirName is a directory
	 */
	OGDF_EXPORT void getEntriesAppend(const char *dirName,
		List<string> &entries,
		const char *pattern = "*");


	//! Returns in \a entries the list of all entries of type \a t contained in directory \a dirName.
	/** The optional argument \a pattern can be used to filter files.
	 *
	 *  \pre \a dirName is a directory
	 */
	OGDF_EXPORT void getEntries(const char *dirName,
		FileType t,
		List<string> &entries,
		const char *pattern = "*");

	//! Appends to \a entries the list of all entries of type \a t contained in directory \a dirName.
	/** The optional argument \a pattern can be used to filter files.
	 *
	 *  \pre \a dirName is a directory
	 */
	OGDF_EXPORT void getEntriesAppend(const char *dirName,
		FileType t,
		List<string> &entries,
		const char *pattern = "*");

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



/**
 * @name Atomic operations
 */
//@{

#ifdef OGDF_SYSTEM_WINDOWS

#define OGDF_MEMORY_BARRIER MemoryBarrier()

//! Atomically decrements (decreases by one) the variable to which \a pX points.
/**
 * @param pX points to the variable to be decremented.
 * @return The resulting decremented value.
 */
inline __int32 atomicDec(__int32 volatile *pX) {
	return (__int32)InterlockedDecrement((LONG volatile *)pX);
}

//! Atomically decrements (decreases by one) the variable to which \a pX points.
/**
 * @param pX points to the variable to be decremented.
 * @return The resulting decremented value.
 */
inline __int64 atomicDec(__int64 volatile *pX) {
	return InterlockedDecrement64(pX);
}

//! Atomically increments (increases by one) the variable to which \a pX points.
/**
 * @param pX points to the variable to be incremented.
 * @return The resulting incremented value.
 */
inline __int32 atomicInc(__int32 volatile *pX) {
	return (__int32)InterlockedIncrement((LONG volatile *)pX);
}

//! Atomically increments (increases by one) the variable to which \a pX points.
/**
 * @param pX points to the variable to be incremented.
 * @return The resulting incremented value.
 */
inline __int64 atomicInc(__int64 volatile *pX) {
	return InterlockedIncrement64(pX);
}


//! Atomically sets the variable pointed to by \a pX to \a value and returns its previous value.
/**
 * @param pX    points to the variable to be modified.
 * @param value is the value to which the variable is set.
 * @return The previous value of the modified variable.
 */
inline __int32 atomicExchange(__int32 volatile *pX, __int32 value) {
	return InterlockedExchange((LONG volatile *)pX, (LONG)value);
}

//! Atomically sets the variable pointed to by \a pX to \a value and returns its previous value.
/**
 * @param pX    points to the variable to be modified.
 * @param value is the value to which the variable is set.
 * @return The previous value of the modified variable.
 */
inline __int64 atomicExchange(__int64 volatile *pX, __int64 value) {
	return InterlockedExchange64(pX, value);
}

//! Atomically sets the variable pointed to by \a pX to \a value and returns its previous value.
/**
 * @param pX    points to the variable to be modified.
 * @param value is the value to which the variable is set.
 * @return The previous value of the modified variable.
 */
template<typename T>
inline T *atomicExchange(T * volatile *pX, T *value) {
	return (T *)InterlockedExchangePointer((PVOID volatile *)pX, value);
}


#if defined(_M_AMD64)
//! Atomically subtracts \a value from the variable to which \a pX points.
/**
 * @param pX    points to the variable to be modified.
 * @param value is the value to be subtracted form the variable.
 * @return The resulting value of the modified variable.
 */
inline __int32 atomicSub(__int32 volatile *pX, __int32 value) {
	return (__int32)InterlockedAdd((LONG volatile *)pX, -value);
}

//! Atomically subtracts \a value from the variable to which \a pX points.
/**
 * @param pX    points to the variable to be modified.
 * @param value is the value to be subtracted form the variable.
 * @return The resulting value of the modified variable.
 */
inline __int64 atomicSub(__int64 volatile *pX, __int64 value) {
	return InterlockedAdd64(pX, -value);
}

//! Atomically adds \a value to the variable to which \a pX points.
/**
 * @param pX    points to the variable to be modified.
 * @param value is the value to be added to the variable.
 * @return The resulting value of the modified variable.
 */
inline __int32 atomicAdd(__int32 volatile *pX, __int32 value) {
	return (__int32)InterlockedAdd((LONG volatile *)pX, value);
}

//! Atomically adds \a value to the variable to which \a pX points.
/**
 * @param pX    points to the variable to be modified.
 * @param value is the value to be added to the variable.
 * @return The resulting value of the modified variable.
 */
inline __int64 atomicAdd(__int64 volatile *pX, __int64 value) {
	return InterlockedAdd64(pX, value);
}

#endif


#else

#define OGDF_MEMORY_BARRIER __sync_synchronize()

inline __int32 atomicDec(__int32 volatile *pX) {
	return  __sync_sub_and_fetch(pX, 1);
}

inline __int64 atomicDec(__int64 volatile *pX) {
	return  __sync_sub_and_fetch(pX, 1);
}

inline __int32 atomicInc(__int32 volatile *pX) {
	return  __sync_add_and_fetch(pX, 1);
}

inline __int64 atomicInc(__int64 volatile *pX) {
	return  __sync_add_and_fetch(pX, 1);
}

inline __int32 atomicSub(__int32 volatile *pX, __int32 value) {
	return __sync_sub_and_fetch(pX, value);
}

inline __int64 atomicSub(__int64 volatile *pX, __int64 value) {
	return __sync_sub_and_fetch(pX, value);
}

inline __int32 atomicAdd(__int32 volatile *pX, __int32 value) {
	return __sync_add_and_fetch(pX, value);
}

inline __int64 atomicAdd(__int64 volatile *pX, __int64 value) {
	return __sync_add_and_fetch(pX, value);
}

inline __int32 atomicExchange(__int32 volatile *pX, __int32 value) {
	return __sync_lock_test_and_set(pX, value);
}

inline __int64 atomicExchange(__int64 volatile *pX, __int64 value) {
	return __sync_lock_test_and_set(pX, value);
}

template<typename T>
inline T *atomicExchange(T * volatile *pX, T *value) {
	return __sync_lock_test_and_set(pX, value);
}

#endif
//@}

} // end namespace ogdf


/**
 * @name C++11 support
 * OGDF uses some new functions defined in the C++11 standard; if these functions are not
 * supported by the compiler, they are be defined by OGDF.
 */
//@{

#ifndef OGDF_HAVE_CPP11

int                stoi  (const string& _Str, size_t *_Idx = 0, int _Base = 10);
long long          stoll (const string& _Str, size_t *_Idx = 0, int _Base = 10);
unsigned long      stoul (const string& _Str, size_t *_Idx = 0, int _Base = 10);
unsigned long long stoull(const string& _Str, size_t *_Idx = 0, int _Base = 10);

float       stof (const string& _Str, size_t *_Idx = 0);
double      stod (const string& _Str, size_t *_Idx = 0);
long double stold(const string& _Str, size_t *_Idx = 0);

string to_string(long long _Val);
string to_string(unsigned long long _Val);
string to_string(long double _Val);

#else

using std::stoi;
using std::stoll;
using std::stoul;
using std::stoull;

using std::stof;
using std::stod;
using std::stold;

using std::to_string;

#endif


#if !defined(OGDF_HAVE_CPP11) || (defined(_MSC_VER) && (_MSC_VER < 1700))

inline string to_string(int           value) { return to_string( (long long)          value); }
inline string to_string(long          value) { return to_string( (long long)          value); }
inline string to_string(unsigned int  value) { return to_string( (unsigned long long) value); }
inline string to_string(unsigned long value) { return to_string( (unsigned long long) value); }
inline string to_string(float         value) { return to_string( (long double)        value); }
inline string to_string(double        value) { return to_string( (long double)        value); }

#endif

// in C++11 we can directly pass a string as filename
#ifdef OGDF_HAVE_CPP11
#define OGDF_STRING_OPEN(filename) (filename)
#else
#define OGDF_STRING_OPEN(filename) (filename).c_str()
#endif
//@}


#endif
