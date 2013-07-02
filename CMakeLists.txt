#
# CMake listfile to specify the build process, see:
# http://www.cmake.org/cmake/help/documentation.html
#

# Require CMake 2.8.9.
cmake_minimum_required(VERSION 2.8.9)

# Check for polluted source tree.
if(EXISTS ${CMAKE_SOURCE_DIR}/CMakeCache.txt OR
	EXISTS ${CMAKE_SOURCE_DIR}/CMakeFiles)
	message(FATAL_ERROR
		"Source directory is polluted:"
		"\n  * remove CMakeCache.txt"
		"\n  * remove CMakeFiles directory")
endif()

# Suppress in-source builds.
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
	message(FATAL_ERROR
		"CMake generation is not allowed within the source directory:"
		"\n  * mkdir build"
		"\n  * cd build"
		"\n  * Unix-like: cmake -G \"Unix Makefiles\" .."
		"\n  * Windows: cmake -G \"Visual Studio 10\" ..")
endif()

# Include macros.
include(${CMAKE_CURRENT_SOURCE_DIR}/config/source_dirs.cmake)

# Enable C++11
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR EMSCRIPTEN)
	list(APPEND CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	list(APPEND CMAKE_CXX_FLAGS "-std=gnu++0x -fno-operator-names")
endif()

# Enable Emscripten, if requested (-DEMSCRIPTEN=1).
if(EMSCRIPTEN)
	if(NOT IS_DIRECTORY "$ENV{EMSCRIPTEN}")
		message(FATAL_ERROR "Invalid environment variable EMSCRIPTEN.")
	endif()
	# Setup tools.
	set(CMAKE_C_COMPILER "emcc")
	set(CMAKE_CXX_COMPILER "em++")
	set(CMAKE_AR "emar")
	set(CMAKE_RANLIB "emranlib")
	# Setup commands for static libraries.
	set(CMAKE_C_ARCHIVE_CREATE "${CMAKE_C_COMPILER} -o <TARGET> -emit-llvm <LINK_FLAGS> <OBJECTS>")
	set(CMAKE_CXX_ARCHIVE_CREATE "${CMAKE_CXX_COMPILER} -o <TARGET> -emit-llvm <LINK_FLAGS> <OBJECTS>")
	# Setup file suffixes for static libraries.
	set(CMAKE_STATIC_LIBRARY_SUFFIX_C ".bc")
	set(CMAKE_STATIC_LIBRARY_SUFFIX_CXX ".bc")
	# Warn on missing symbols.
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -s WARN_ON_UNDEFINED_SYMBOLS=1")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s WARN_ON_UNDEFINED_SYMBOLS=1")
	# Setup find root path.
	set(CMAKE_FIND_ROOT_PATH "$ENV{EMSCRIPTEN}")
	set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
	set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
	set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE BOTH)
	set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)
	# Overwrite some things to make it compile.
	add_definitions(-Dsrandom=srand -Drandom=rand)
	set(OPENMP_ENABLED FALSE)
	set(DOGDF_MEMORY_MANAGER "OGDF_MEMORY_MALLOC_TS")
endif()

# Setup root project.
project(OGDF)

# Change some defaults.
if(MSVC)
	# /MP Enable build with multiple processes.
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
	# _VARIADIC_MAX=10 Set max. amount of variadic arguments.
	# _CRT_SECURE_NO_WARNINGS Disable CRT warnings.
	# _SCL_SECURE_NO_WARNINGS=10 Disable SCL warnings.
	add_definitions(
		"/D_VARIADIC_MAX=10"
		"/D_CRT_SECURE_NO_WARNINGS"
		"/D_SCL_SECURE_NO_WARNINGS")
endif()

# Setup OpenMP, if requested.
set(OPENMP_ENABLED TRUE CACHE BOOL "Enable OpenMP")
if(OPENMP_ENABLED)
	find_package(OpenMP REQUIRED)
	if(OPENMP_FOUND)
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	endif()
endif()

# We can't use add_subdirectory() here, because of OGDF's special directory
# structure. Thus we'll have to include projects (order is important!).
include(${CMAKE_CURRENT_SOURCE_DIR}/config/coin.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/config/ogdf.cmake)
if(NOT EMSCRIPTEN)
	include(${CMAKE_CURRENT_SOURCE_DIR}/config/ogdf-test.cmake)
endif()
