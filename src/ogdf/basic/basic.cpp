/*
 * $Revision: 3830 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 09:55:21 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

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
#include <time.h>

// Windows includes
#ifdef OGDF_SYSTEM_WINDOWS
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
#ifdef OGDF_USE_THREAD_POOL
		ogdf::Thread::initPool();
#endif
	}
}

Initialization::~Initialization()
{
	if (--s_count == 0) {
#ifdef OGDF_USE_THREAD_POOL
		ogdf::Thread::cleanupPool();
#endif
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


double usedTime(double& T)
{
	double t = T;

#ifdef OGDF_SYSTEM_WINDOWS
	FILETIME creationTime;
	FILETIME exitTime;
	FILETIME kernelTime;
	FILETIME userTime;

	BOOL res = GetProcessTimes(GetCurrentProcess(), &creationTime, &exitTime, &kernelTime, &userTime);
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


#ifdef OGDF_SYSTEM_WINDOWS

bool isFile(const char *fileName)
{
	DWORD att = GetFileAttributes(fileName);

	if (att == 0xffffffff) return false;
	return (att & FILE_ATTRIBUTE_DIRECTORY) == 0;
}


bool isDirectory(const char *fileName)
{
	DWORD att = GetFileAttributes(fileName);

	if (att == 0xffffffff) return false;
	return (att & FILE_ATTRIBUTE_DIRECTORY) != 0;
}


bool changeDir(const char *dirName)
{
	return (_chdir(dirName) == 0);
}


void getEntriesAppend(const char *dirName,
		FileType t,
		List<string> &entries,
		const char *pattern)
{
	OGDF_ASSERT(isDirectory(dirName));

	string filePattern = string(dirName) + "\\" + pattern;

	WIN32_FIND_DATA findData;
	HANDLE handle = FindFirstFile(filePattern.c_str(), &findData);

	if (handle != INVALID_HANDLE_VALUE)
	{
		do {
			DWORD isDir = (findData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY);
			if(isDir && (
				strcmp(findData.cFileName,".") == 0 ||
				strcmp(findData.cFileName,"..") == 0)
			)
				continue;

			if (t == ftEntry || (t == ftFile && !isDir) ||
				(t == ftDirectory && isDir))
			{
				entries.pushBack(findData.cFileName);
			}
		} while(FindNextFile(handle, &findData));

		FindClose(handle);
	}
}
#endif

#ifdef OGDF_SYSTEM_UNIX

bool isDirectory(const char *fname)
{
	struct stat stat_buf;

	if (stat(fname,&stat_buf) != 0)
		return false;
	return (stat_buf.st_mode & S_IFMT) == S_IFDIR;
}

bool isFile(const char *fname)
{
	struct stat stat_buf;

	if (stat(fname,&stat_buf) != 0)
		return false;
	return (stat_buf.st_mode & S_IFMT) == S_IFREG;
}

bool changeDir(const char *dirName)
{
	return (chdir(dirName) == 0);
}

void getEntriesAppend(const char *dirName,
	FileType t,
	List<string> &entries,
	const char *pattern)
{
	OGDF_ASSERT(isDirectory(dirName));

	DIR* dir_p = opendir(dirName);

	dirent* dir_e;
	while ( (dir_e = readdir(dir_p)) != NULL )
	{
		const char *fname = dir_e->d_name;
		if (pattern != 0 && fnmatch(pattern,fname,0)) continue;

		string fullName = string(dirName) + "/" + fname;

		bool isDir = isDirectory(fullName.c_str());
		if(isDir && (
			strcmp(fname,".") == 0 ||
			strcmp(fname,"..") == 0)
			)
			continue;

		if (t == ftEntry || (t == ftFile && !isDir) ||
			(t == ftDirectory && isDir))
		{
			entries.pushBack(fname);
		}
	}

	closedir(dir_p);
}
#endif


void getEntries(const char *dirName,
		FileType t,
		List<string> &entries,
		const char *pattern)
{
	entries.clear();
	getEntriesAppend(dirName, t, entries, pattern);
}


void getFiles(const char *dirName,
	List<string> &files,
	const char *pattern)
{
	getEntries(dirName, ftFile, files, pattern);
}


void getSubdirs(const char *dirName,
	List<string> &subdirs,
	const char *pattern)
{
	getEntries(dirName, ftDirectory, subdirs, pattern);
}


void getEntries(const char *dirName,
	List<string> &entries,
	const char *pattern)
{
	getEntries(dirName, ftEntry, entries, pattern);
}


void getFilesAppend(const char *dirName,
	List<string> &files,
	const char *pattern)
{
	getEntriesAppend(dirName, ftFile, files, pattern);
}


void getSubdirsAppend(const char *dirName,
	List<string> &subdirs,
	const char *pattern)
{
	getEntriesAppend(dirName, ftDirectory, subdirs, pattern);
}


void getEntriesAppend(const char *dirName,
	List<string> &entries,
	const char *pattern)
{
	getEntriesAppend(dirName, ftEntry, entries, pattern);
}


} // end namespace ogdf


//---------------------------------------------------------
// additional C++11 functions
//---------------------------------------------------------

#ifndef OGDF_HAVE_CPP11

#include <iomanip>
#include <sstream>


int stoi(
	const string& _Str,
	size_t *_Idx,
	int _Base)
{
	std::istringstream sstr(_Str);
	int val;
	sstr >> std::setbase(_Base) >> val;
	if(_Idx) *_Idx = sstr.tellg();
	return val;
}


long long stoll(
	const string& _Str,
	size_t *_Idx,
	int _Base)
{
	std::istringstream sstr(_Str);
	long long val;
	sstr >> std::setbase(_Base) >> val;
	if(_Idx) *_Idx = sstr.tellg();
	return val;
}


unsigned long stoul(
	const string& _Str,
	size_t *_Idx,
	int _Base)
{
	std::istringstream sstr(_Str);
	unsigned long val;
	sstr >> std::setbase(_Base) >> val;
	if(_Idx) *_Idx = sstr.tellg();
	return val;
}


unsigned long long stoull(
	const string& _Str,
	size_t *_Idx,
	int _Base)
{
	std::istringstream sstr(_Str);
	unsigned long long val;
	sstr >> std::setbase(_Base) >> val;
	if(_Idx) *_Idx = sstr.tellg();
	return val;
}


float stof(
	const string& _Str,
	size_t *_Idx)
{
	std::istringstream sstr(_Str);
	float val;
	sstr >> val;
	if(_Idx) *_Idx = sstr.tellg();
	return val;
}


double stod(
	const string& _Str,
	size_t *_Idx)
{
	std::istringstream sstr(_Str);
	double val;
	sstr >> val;
	if(_Idx) *_Idx = sstr.tellg();
	return val;
}


long double stold(
	const string& _Str,
	size_t *_Idx)
{
	std::istringstream sstr(_Str);
	long double val;
	sstr >> val;
	if(_Idx) *_Idx = sstr.tellg();
	return val;
}


string to_string(long long _Val)
{
	std::ostringstream sstr;
	sstr << _Val;
	return sstr.str();
}


string to_string(unsigned long long _Val)
{
	std::ostringstream sstr;
	sstr << _Val;
	return sstr.str();
}


string to_string(long double _Val)
{
	std::ostringstream sstr;
	sstr << _Val;
	return sstr.str();
}

#endif
