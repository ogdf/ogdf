#pragma once

//! Disables (deletes) copy construction and assignment for class \p cls.
#define NO_COPY(cls)               \
	cls(const cls& copy) = delete; \
	cls& operator=(const cls& copy) = delete;

//! Disables (deletes) move construction and assignment for class \p cls.
#define NO_MOVE(cls)          \
	cls(cls&& move) = delete; \
	cls& operator=(cls&& move) = delete;

//! Provides default copy construction and assignment for class \p cls.
#define DEFAULT_COPY(cls)           \
	cls(const cls& copy) = default; \
	cls& operator=(const cls& copy) = default;

//! Provides default move construction and assignment for class \p cls.
#define DEFAULT_MOVE(cls)      \
	cls(cls&& move) = default; \
	cls& operator=(cls&& move) = default;

//! Disables copy construction and assignment for class \p cls by throwing a runtime_error.
#define NO_RT_COPY(cls)                                                              \
	cls(const cls& copy) { throw std::runtime_error("class cls is not copyable"); }; \
	cls& operator=(const cls& copy) { throw std::runtime_error("class cls is not copyable"); };

//! Disables move construction and assignment for class \p cls by throwing a runtime_error.
#define NO_RT_MOVE(cls)                                                        \
	cls(cls&& move) { throw std::runtime_error("class cls is not movable"); }; \
	cls& operator=(cls&& move) { throw std::runtime_error("class cls is not movable"); };

//! Declares the copy constructor for class \p cls.
#define COPY_CONSTR(cls) cls(const cls& copy)
//! Declares the copy assignment operation for class \p cls.
#define COPY_OP(cls) cls& operator=(const cls& copy)
//! Declares the move constructor for class \p cls.
#define MOVE_CONSTR(cls) cls(cls&& move) noexcept : cls()
//! Declares the move assignment operation for class \p cls.
#define MOVE_OP(cls) cls& operator=(cls&& move) noexcept
//! Declares the swap function for class \p cls outside of the definition of \p cls.
#define SWAP_OP(cls) friend void swap(cls& first, cls& second) noexcept

//! Provide move construct/assign and copy assign given there is a copy constructor and swap.
/**
 * Only requires custom implementations of COPY_CONSTR and SWAP_OP and automatically provides
 * COPY_OP, MOVE_CONSTR, and MOVE_OP.
 */
#define COPY_MOVE_BY_SWAP(cls)          \
	cls& operator=(cls copy_by_value) { \
		using std::swap;                \
		swap(*this, copy_by_value);     \
	}                                   \
	MOVE_CONSTR(cls) {                  \
		using std::swap;                \
		swap(*this, move);              \
	}                                   \
	MOVE_OP(cls) {                      \
		using std::swap;                \
		swap(*this, move);              \
	}
