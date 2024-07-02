/** \file
 * \brief Utility macros for declaring copy and move constructors and assignment operations.
 *
 * \author Simon D. Fink
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

/* This needs to be synchronized with doc/ogdf-doxygen.cfg:EXPAND_AS_DEFINED */

//! Explicitly disables (deletes) copy construction and assignment for class \p cls.
#define OGDF_NO_COPY(cls)          \
	cls(const cls& copy) = delete; \
	cls& operator=(const cls& copy) = delete;

//! Explicitly disables (deletes) move construction and assignment for class \p cls.
#define OGDF_NO_MOVE(cls)     \
	cls(cls&& move) = delete; \
	cls& operator=(cls&& move) = delete;

//! Explicitly provides default copy construction and assignment for class \p cls.
#define OGDF_DEFAULT_COPY(cls)      \
	cls(const cls& copy) = default; \
	cls& operator=(const cls& copy) = default;

//! Explicitly provides default move construction and assignment for class \p cls.
#define OGDF_DEFAULT_MOVE(cls)          \
	cls(cls&& move) noexcept = default; \
	cls& operator=(cls&& move) noexcept = default;

//! Declares the copy constructor for class \p cls.
#define OGDF_COPY_CONSTR(cls) cls(const cls& copy)
//! Declares the copy assignment operation for class \p cls. Don't forget to return \p *this;
#define OGDF_COPY_OP(cls) cls& operator=(const cls& copy)
//! Declares the move constructor for class \p cls.
#define OGDF_MOVE_CONSTR(cls) cls(cls&& move) noexcept : cls()
//! Declares the move assignment operation for class \p cls.
#define OGDF_MOVE_OP(cls) cls& operator=(cls&& move) noexcept
//! Declares the swap function for class \p cls.
#define OGDF_SWAP_OP(cls) friend void swap(cls& first, cls& second) noexcept

//! Provide move construct/assign and copy assign given there is a copy constructor and swap.
/**
 * Requires custom implementations of OGDF_COPY_CONSTR and OGDF_SWAP_OP (may not be the std::swap via move to temporary)
 * and automatically provides OGDF_COPY_OP, OGDF_MOVE_CONSTR, and OGDF_MOVE_OP.
 *
 * See https://stackoverflow.com/a/3279550 for more details on this idiom.
 * Note that your swap implementation is expected to be noexcept (as declared e.g. by OGDF_SWAP_OP)
 * as all methods declared by this macro are also noexcept (see https://stackoverflow.com/a/7628345).
 */
#define OGDF_COPY_MOVE_BY_SWAP(cls)                       \
	/*! Assign this cls to be a copy of \p copy_by_value.
	    Internally, this will use the \ref cls ## (const cls&) "copy constructor" to create a temporary copy
	    of the argument \p copy_by_value (as it is passed by value) and then \ref swap(cls&, cls&) "swap" this object
	    instance with the temporary copy.
	    If the assigned-from object can be moved, the move constructor will be automatically used instead of copying.
	    Note that this method thereby covers both copy-assignment and move-assignment.
	    See https://stackoverflow.com/a/3279550 for more details on this idiom.
	    This method is noexcept as a potentially-throwing copy constructor call is made within the context
	    of the caller (see https://stackoverflow.com/a/7628345) and swap is expected to be noexcept. */ \
	cls& operator=(cls copy_by_value) noexcept {          \
		swap(*this, copy_by_value);                       \
		return *this;                                     \
	}                                                     \
	OGDF_MOVE_CONSTR(cls) { swap(*this, move); }
