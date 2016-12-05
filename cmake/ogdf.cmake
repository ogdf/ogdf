# OGDF (library only) CMake configuration

# cache configuration
option(BUILD_SHARED_LIBS "Whether to build shared libraries instead of static ones." OFF)
set(OGDF_MEMORY_MANAGER "POOL_TS" CACHE STRING "Memory manager to be used.")
set_property(CACHE OGDF_MEMORY_MANAGER PROPERTY STRINGS POOL_TS POOL_NTS MALLOC_TS)
option(OGDF_DEBUG "Whether to include OGDF assertions in Debug mode (increased runtime)." ON)
mark_as_advanced(OGDF_DEBUG)
option(OGDF_USE_ASSERT_EXCEPTIONS "Whether to throw an exception on failed assertions." OFF)
set(OGDF_USE_ASSERT_EXCEPTIONS_WITH_STACK_TRACE "OFF" CACHE STRING "Which library (libdw, libbdf, libunwind) to use in case a stack trace should be written to a failed assertion exceptions's what(). Library must be found by CMake to be able to use it.")
set_property(CACHE OGDF_USE_ASSERT_EXCEPTIONS_WITH_STACK_TRACE PROPERTY STRINGS "OFF")

# find available packages for stack traces
find_package(Libdw)
if(LIBDW_FOUND)
  set_property(CACHE OGDF_USE_ASSERT_EXCEPTIONS_WITH_STACK_TRACE APPEND PROPERTY STRINGS "ON_LIBDW")
endif()
find_package(Libbfd)
if(LIBBFD_FOUND)
  set_property(CACHE OGDF_USE_ASSERT_EXCEPTIONS_WITH_STACK_TRACE APPEND PROPERTY STRINGS "ON_LIBBFD")
endif()
find_package(Libunwind)
if(LIBUNWIND_FOUND)
  set_property(CACHE OGDF_USE_ASSERT_EXCEPTIONS_WITH_STACK_TRACE APPEND PROPERTY STRINGS "ON_LIBUNWIND")
endif()
set(OGDF_EXTRA_CXX_FLAGS "${available_default_warning_flags}" CACHE STRING "Extra compiler flags for compiling OGDF and tests.")
mark_as_advanced(OGDF_EXTRA_CXX_FLAGS)

# compilation
file(GLOB_RECURSE OGDF_SOURCES src/ogdf/*.cpp include/ogdf/*.h)
if(NOT OGDF_COMPILE_LEGACY)
  file(GLOB_RECURSE OGDF_LEGACY_SOURCES src/ogdf/legacy/*.cpp)
  foreach(legacyFile ${OGDF_LEGACY_SOURCES})
    list(REMOVE_ITEM OGDF_SOURCES ${legacyFile})
  endforeach()
endif()
add_library(OGDF ${OGDF_SOURCES})
target_link_libraries(OGDF COIN)
group_files(OGDF_SOURCES "ogdf")
target_compile_features(OGDF PUBLIC cxx_range_for)
if(COIN_EXTERNAL_SOLVER_INCLUDE_DIRECTORIES)
  target_include_directories(OGDF SYSTEM PUBLIC ${COIN_EXTERNAL_SOLVER_INCLUDE_DIRECTORIES})
endif()
set_property(TARGET OGDF APPEND_STRING PROPERTY COMPILE_FLAGS ${OGDF_EXTRA_CXX_FLAGS})

# set OGDF_INSTALL for shared libraries
if(BUILD_SHARED_LIBS)
  target_compile_definitions(OGDF PRIVATE OGDF_INSTALL)
endif()

# special variable handling for stack trace feature
if(OGDF_DEBUG AND OGDF_USE_ASSERT_EXCEPTIONS)
endif()

# autogen header variables for debug mode
set(SHOW_STACKTRACE 0)
if(OGDF_DEBUG AND OGDF_USE_ASSERT_EXCEPTIONS)
  include(check-pretty-function)
  if(compiler_has_pretty_function)
    set(OGDF_FUNCTION_NAME "__PRETTY_FUNCTION__")
  else() # fallback to C++11 standard
    set(OGDF_FUNCTION_NAME "__func__")
  endif()
  if(OGDF_USE_ASSERT_EXCEPTIONS_WITH_STACK_TRACE MATCHES LIBDW)
    set(SHOW_STACKTRACE 1)
    set(BACKWARD_HAS_DW 1)
  elseif(OGDF_USE_ASSERT_EXCEPTIONS_WITH_STACK_TRACE MATCHES LIBBFD)
    set(SHOW_STACKTRACE 1)
    set(BACKWARD_HAS_BFD 1)
  elseif(OGDF_USE_ASSERT_EXCEPTIONS_WITH_STACK_TRACE MATCHES LIBUNWIND)
    set(SHOW_STACKTRACE 1)
    set(BACKWARD_HAS_UNWIND 1)
  endif()
endif()

# autogen header variables if libs are shared
if(BUILD_SHARED_LIBS)
  set(OGDF_DLL 1)
endif()

# autogen header variables for SSE3
include(check-sse3)
if(has_sse3_intrin)
  set(OGDF_SSE3_EXTENSIONS <intrin.h>)
elseif(has_sse3_pmmintrin)
  set(OGDF_SSE3_EXTENSIONS <pmmintrin.h>)
else()
  message(STATUS "SSE3 could not be activated")
endif()

# add stack trace include paths
if(BACKWARD_HAS_DW)
  target_include_directories(OGDF SYSTEM PUBLIC ${LIBDW_INCLUDE_DIR})
elseif(BACKWARD_HAS_BFD)
  target_include_directories(OGDF SYSTEM PUBLIC ${LIBBFD_INCLUDE_DIR} ${LIBDL_INCLUDE_DIR})
elseif(BACKWARD_HAS_UNWIND)
  target_include_directories(OGDF SYSTEM PUBLIC ${LIBUNWIND_INCLUDE_DIR})
endif()
