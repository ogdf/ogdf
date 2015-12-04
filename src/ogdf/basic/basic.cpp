/** \file
 * \brief Implementation of basic functionality (incl. file and
 * directory handling)
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

#include <ogdf/basic/Thread.h>
#include <ogdf/basic/List.h>
#include <random>

// Windows includes
#ifdef OGDF_SYSTEM_WINDOWS

#define WIN32_EXTRA_LEAN
#define WIN32_LEAN_AND_MEAN
#undef NOMINMAX
#define NOMINMAX
#include <windows.h>

#include <direct.h>
#if defined(_MSC_VER) && defined(UNICODE)
#undef GetFileAttributes
#undef FindFirstFile
#undef FindNextFile
#define GetFileAttributes  GetFileAttributesA
#define FindFirstFile  FindFirstFileA
#define WIN32_FIND_DATA WIN32_FIND_DATAA
#define FindNextFile  FindNextFileA
#endif
#endif

#ifdef __BORLANDC__
#define _chdir chdir
#endif

// Unix includes
#ifdef OGDF_SYSTEM_UNIX
#include <cstring>
#include <unistd.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/stat.h>
#include <fnmatch.h>
#endif

#ifdef OGDF_DLL

#ifdef OGDF_SYSTEM_WINDOWS

#ifdef __MINGW32__
extern "C"
#endif
BOOL APIENTRY DllMain(HANDLE hModule, DWORD  ul_reason_for_call, LPVOID lpReserved)
{
	switch (ul_reason_for_call)
	{
	case DLL_PROCESS_ATTACH:
		ogdf::PoolMemoryAllocator::init();
		ogdf::System::init();
		break;

	case DLL_THREAD_ATTACH:
	case DLL_THREAD_DETACH:
		break;

	case DLL_PROCESS_DETACH:
		ogdf::PoolMemoryAllocator::cleanup();
		break;
	}
	return TRUE;
}

#else

void __attribute__ ((constructor)) my_load(void)
{
	ogdf::PoolMemoryAllocator::init();
	ogdf::System::init();
}

void __attribute__ ((destructor)) my_unload(void)
{
	ogdf::PoolMemoryAllocator::cleanup();
}

#endif

#else

namespace ogdf {

//static int variables are automatically initialized with 0
int Initialization::s_count;

Initialization::Initialization()
{
	if (s_count++ == 0) {
		ogdf::PoolMemoryAllocator::init();
		ogdf::System::init();
	}
}

Initialization::~Initialization()
{
	if (--s_count == 0) {
		ogdf::PoolMemoryAllocator::cleanup();
	}
}

} // namespace ogdf

#endif

namespace ogdf {

inline bool charCompareIgnoreCase(char a, char b)
{
	return (toupper(a) == toupper(b));
}

void removeTrailingWhitespace(std::string &str)
{
	std::size_t found = str.find_last_not_of(" \t\v\f\n\r");
	if (found != std::string::npos) {
		str.erase(found+1);
	} else { // string consists only of whitespacae
		str.clear();
	}
}

bool equalIgnoreCase(const string &str1, const string &str2)
{
	return (str1.size() == str2.size() &&
		std::equal(str1.begin(), str1.end(), str2.begin(), charCompareIgnoreCase));
}

bool prefixIgnoreCase(const string &prefix, const string &str)
{
	string::size_type len = prefix.length();
	return (str.size() >= len &&
		std::equal(prefix.begin(), prefix.end(), str.begin(), charCompareIgnoreCase));
}

// debug level (in debug build only)
#ifdef OGDF_DEBUG
DebugLevel debugLevel;
#endif

static std::mt19937 s_random;

#ifndef OGDF_MEMORY_POOL_NTS
static std::mutex s_randomMutex;
#endif

long unsigned int randomSeed()
{
#ifndef OGDF_MEMORY_POOL_NTS
	std::lock_guard<std::mutex> guard(s_randomMutex);
#endif
	return 7*s_random()+3;  // do not directly return seed, add a bit of variation
}

void setSeed(int val)
{
	s_random.seed(val);
}

int randomNumber(int low, int high)
{
	OGDF_ASSERT(low <= high);

	std::uniform_int_distribution<> dist(low,high);

#ifndef OGDF_MEMORY_POOL_NTS
	std::lock_guard<std::mutex> guard(s_randomMutex);
#endif
	return dist(s_random);
}

double randomDouble(double low, double high)
{
	OGDF_ASSERT(low <= high);

	std::uniform_real_distribution<> dist(low,high);

#ifndef OGDF_MEMORY_POOL_NTS
	std::lock_guard<std::mutex> guard(s_randomMutex);
#endif
	return dist(s_random);
}

double usedTime(double& T)
{
	double t = T;

#ifdef OGDF_SYSTEM_WINDOWS
	FILETIME creationTime;
	FILETIME exitTime;
	FILETIME kernelTime;
	FILETIME userTime;

	GetProcessTimes(GetCurrentProcess(), &creationTime, &exitTime, &kernelTime, &userTime);
	ULARGE_INTEGER user;
	user.LowPart = userTime.dwLowDateTime;
	user.HighPart = userTime.dwHighDateTime;
	T = double(user.QuadPart) * 0.0000001;

#else
	struct tms now;
	times (&now);
	T = (double)now.tms_utime / (double)sysconf(_SC_CLK_TCK);
#endif

	return T - t;
}

} // end namespace ogdf
