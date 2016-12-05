/** \file
 * \brief Declaration of deprecated file system functions
 *
 * \author Stephan Beyer (responsible only for deprecation)
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

#include <ogdf/basic/List.h>

namespace ogdf {

/**
 * @addtogroup file-system
 */
//@{

//! The type of an entry in a directory.
enum FileType {
	ftEntry,     /**< file or directory */
	ftFile,      /**< file */
	ftDirectory  /**< directory */
};

//! Returns true iff \a fileName is a regular file (not a directory).
OGDF_DEPRECATED("Please use another library (or C++17) for filesystem functions.")
OGDF_EXPORT bool isFile(const char *fileName);

//! Returns true iff \a fileName is a directory.
OGDF_DEPRECATED("Please use another library (or C++17) for filesystem functions.")
OGDF_EXPORT bool isDirectory(const char *fileName);

//! Changes current directory to \a dirName; returns true if successful.
OGDF_DEPRECATED("Please use another library (or C++17) for filesystem functions.")
OGDF_EXPORT bool changeDir(const char *dirName);

//! Returns in \a files the list of files in directory \a dirName.
/** The optional argument \a pattern can be used to filter files.
 *
 *  \pre \a dirName is a directory
 */
OGDF_DEPRECATED("Please use another library (or C++17) for filesystem functions.")
OGDF_EXPORT void getFiles(const char *dirName,
                          List<string> &files,
                          const char *pattern = "*");

//! Appends to \a files the list of files in directory \a dirName.
/** The optional argument \a pattern can be used to filter files.
 *
 *  \pre \a dirName is a directory
 */
OGDF_DEPRECATED("Please use another library (or C++17) for filesystem functions.")
OGDF_EXPORT void getFilesAppend(const char *dirName,
                                List<string> &files,
                                const char *pattern = "*");

//! Returns in \a subdirs the list of directories contained in directory \a dirName.
/** The optional argument \a pattern can be used to filter files.
 *
 *  \pre \a dirName is a directory
 */
OGDF_DEPRECATED("Please use another library (or C++17) for filesystem functions.")
OGDF_EXPORT void getSubdirs(const char *dirName,
                            List<string> &subdirs,
                            const char *pattern = "*");

//! Appends to \a subdirs the list of directories contained in directory \a dirName.
/** The optional argument \a pattern can be used to filter files.
 *
 *  \pre \a dirName is a directory
 */
OGDF_DEPRECATED("Please use another library (or C++17) for filesystem functions.")
OGDF_EXPORT void getSubdirsAppend(const char *dirName,
                                  List<string> &subdirs,
                                  const char *pattern = "*");

//! Returns in \a entries the list of all entries contained in directory \a dirName.
/** Entries may be files or directories. The optional argument \a pattern
 *  can be used to filter files.
 *
 *  \pre \a dirName is a directory
 */
OGDF_DEPRECATED("Please use another library (or C++17) for filesystem functions.")
OGDF_EXPORT void getEntries(const char *dirName,
                            List<string> &entries,
                            const char *pattern = "*");

//! Appends to \a entries the list of all entries contained in directory \a dirName.
/** Entries may be files or directories. The optional argument \a pattern
 *  can be used to filter files.
 *
 *  \pre \a dirName is a directory
 */
OGDF_DEPRECATED("Please use another library (or C++17) for filesystem functions.")
OGDF_EXPORT void getEntriesAppend(const char *dirName,
                                  List<string> &entries,
                                  const char *pattern = "*");

//! Returns in \a entries the list of all entries of type \a t contained in directory \a dirName.
/** The optional argument \a pattern can be used to filter files.
 *
 *  \pre \a dirName is a directory
 */
OGDF_DEPRECATED("Please use another library (or C++17) for filesystem functions.")
OGDF_EXPORT void getEntries(const char *dirName,
                            FileType t,
                            List<string> &entries,
                            const char *pattern = "*");

//! Appends to \a entries the list of all entries of type \a t contained in directory \a dirName.
/** The optional argument \a pattern can be used to filter files.
 *
 *  \pre \a dirName is a directory
 */
OGDF_DEPRECATED("Please use another library (or C++17) for filesystem functions.")
OGDF_EXPORT void getEntriesAppend(const char *dirName,
                                  FileType t,
                                  List<string> &entries,
                                  const char *pattern = "*");

//@}

}
