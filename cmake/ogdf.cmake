# OGDF (library only) CMake configuration

# cache configuration
option(BUILD_SHARED_LIBS "Whether to build shared libraries instead of static ones." OFF)
set(OGDF_MEMORY_MANAGER "POOL_TS" CACHE STRING "Memory manager to be used.")
set_property(CACHE OGDF_MEMORY_MANAGER PROPERTY STRINGS POOL_TS POOL_NTS MALLOC_TS)
set(OGDF_DEBUG_MODE "REGULAR" CACHE STRING "Whether to use (heavy) OGDF assertions in debug mode.")
set_property(CACHE OGDF_DEBUG_MODE PROPERTY STRINGS NONE REGULAR HEAVY)
mark_as_advanced(OGDF_DEBUG_MODE)
option(OGDF_USE_ASSERT_EXCEPTIONS "Whether to throw an exception on failed assertions." OFF)
set(OGDF_USE_ASSERT_EXCEPTIONS_WITH_STACK_TRACE "OFF" CACHE
    STRING "Which library (libdw, libbdf, libunwind) to use in case a stack trace should be written \
    to a failed assertion exceptions's what(). Library must be found by CMake to be able to use it.")
if(OGDF_USE_ASSERT_EXCEPTIONS)
  set_property(CACHE OGDF_USE_ASSERT_EXCEPTIONS_WITH_STACK_TRACE PROPERTY STRINGS "OFF")
else()
  unset(OGDF_USE_ASSERT_EXCEPTIONS_WITH_STACK_TRACE CACHE)
endif()
option(OGDF_WARNING_ERRORS "Whether to treat compiler warnings as errors; may break compilation!" OFF)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang" AND OGDF_MEMORY_MANAGER STREQUAL MALLOC_TS)
  option(OGDF_LEAK_CHECK "Whether to use the address sanitizer for the MALLOC_TS memory manager \
  (and for COIN)." OFF)
else()
  unset(OGDF_LEAK_CHECK CACHE)
endif()

# sets debug mode if it was not explicitly set to NONE and we are in a Debug or Multiconfig build (or the passed config argument is Debug)
function(set_debug_mode)
  set(OGDF_DEBUG OFF PARENT_SCOPE)
  if(OGDF_DEBUG_MODE STREQUAL NONE)
    return()
  elseif(ARGC GREATER 0)
    if(ARGV0 MATCHES Debug)
      set(OGDF_DEBUG ON PARENT_SCOPE)
    endif()
  else()
    if(MULTICONFIG_BUILD OR CMAKE_BUILD_TYPE MATCHES Debug)
      set(OGDF_DEBUG ON PARENT_SCOPE)
    endif()
  endif()
endfunction()

if(OGDF_DEBUG_MODE STREQUAL HEAVY)
  set(OGDF_HEAVY_DEBUG ON)
endif()
set_debug_mode()

# find available packages for stack traces
if(OGDF_USE_ASSERT_EXCEPTIONS)
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
endif()

# find CGAL if enabled
if (OGDF_INCLUDE_CGAL)
  find_package(CGAL REQUIRED COMPONENTS Core)
  find_package(OpenMP REQUIRED)
  set(extra_flags "${extra_flags} ${OpenMP_CXX_FLAGS}")
endif()

set(extra_flags_desc "Extra compiler flags for compiling OGDF, tests, and examples")
set(OGDF_EXTRA_CXX_FLAGS "${available_default_warning_flags}" CACHE
    STRING "${extra_flags_desc}.")
set(OGDF_EXTRA_CXX_FLAGS_DEBUG "${available_default_warning_flags_debug}" CACHE
    STRING "${extra_flags_desc} applied only when compiling in debug mode.")
set(OGDF_EXTRA_CXX_FLAGS_RELEASE "${available_default_warning_flags_release}" CACHE
    STRING "${extra_flags_desc} applied only when not compiling in debug mode.")
mark_as_advanced(OGDF_EXTRA_CXX_FLAGS)
mark_as_advanced(OGDF_EXTRA_CXX_FLAGS_DEBUG)
mark_as_advanced(OGDF_EXTRA_CXX_FLAGS_RELEASE)

# static analysis using clang-tidy
option(OGDF_ENABLE_CLANG_TIDY "Enable static analysis using clang-tidy" OFF)

if(OGDF_ENABLE_CLANG_TIDY)
  find_program(CLANG_TIDY clang-tidy)
  if(CLANG_TIDY)
    set(CMAKE_CXX_CLANG_TIDY ${CLANG_TIDY};-extra-arg=-Wno-unknown-warning-option)
  else()
    message(WARNING "clang-tidy not found!")
  endif()
endif()

# compilation
file(GLOB_RECURSE ogdf_headers include/ogdf/*.h)
file(GLOB_RECURSE ogdf_sources src/ogdf/*.cpp)
add_library(OGDF "${ogdf_sources}")
set_property(TARGET OGDF PROPERTY VERSION ${PROJECT_VERSION})
target_link_libraries(OGDF PUBLIC COIN)
group_files(ogdf_sources "ogdf")
group_files(ogdf_headers "ogdf")
target_compile_features(OGDF PUBLIC cxx_std_${CMAKE_CXX_STANDARD})

target_include_directories(OGDF PUBLIC # for the autogen header
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include/ogdf-$<IF:$<CONFIG:Debug>,debug,release>>
  $<INSTALL_INTERFACE:include/ogdf-$<IF:$<CONFIG:Debug>,debug,release>>)
target_include_directories(OGDF PUBLIC # for the general include files
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
if(COIN_EXTERNAL_SOLVER_INCLUDE_DIRECTORIES)
  target_include_directories(OGDF SYSTEM PUBLIC ${COIN_EXTERNAL_SOLVER_INCLUDE_DIRECTORIES})
endif()

function (add_ogdf_extra_flags TARGET_NAME)
  set(extra_flags ${OGDF_EXTRA_CXX_FLAGS})
  if(CMAKE_BUILD_TYPE STREQUAL Debug)
    set(extra_flags  "${extra_flags} ${OGDF_EXTRA_CXX_FLAGS_DEBUG}")
  else()
    set(extra_flags  "${extra_flags} ${OGDF_EXTRA_CXX_FLAGS_RELEASE}")
  endif()
  if(OGDF_WARNING_ERRORS)
    set(extra_flags "${warnings_as_errors_flag} ${extra_flags}")
  endif()
  if(OGDF_LEAK_CHECK)
    set(leak_flag_list "-fsanitize=address" "-fno-omit-frame-pointer")

    # apple silicon needs explicit linking with the right asan library
    if(APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")
      if(BUILD_SHARED_LIBS)
        list(APPEND leak_flag_list "-shared-libsan")
      else()
        # list(APPEND leak_flag_list "-static-libsan") # not available
        message(WARNING "ASAN requested via OGDF_LEAK_CHECK=ON, but static ASAN is not supported on macOS arm64. ASAN flags will be skipped.")
        set(leak_flag_list "")
      endif()
    endif()

    string(REPLACE ";" " " leak_flags "${leak_flag_list}")
    set(extra_flags "${extra_flags} ${leak_flags}")
    set_property(TARGET ${TARGET_NAME} APPEND_STRING PROPERTY LINK_FLAGS " ${leak_flags} ")
    # If OGDF is compiled with ASAN, compile COIN with ASAN as well.
    # This avoids container-overflow false positives.
    target_compile_options(COIN PRIVATE ${leak_flag_list})
    target_link_options(COIN PRIVATE ${leak_flag_list})
  endif()
  set_property(TARGET ${TARGET_NAME} APPEND_STRING PROPERTY COMPILE_FLAGS " ${extra_flags} ")
endfunction()

add_ogdf_extra_flags(OGDF)

# do not count warnings in external libraries as errors
file(GLOB_RECURSE lib_headers include/ogdf/lib/*.h)
file(GLOB_RECURSE lib_sources src/ogdf/lib/*.cpp)
file(GLOB_RECURSE bandit_headers test/include/bandit/*.h)
set(tinydir_headers "test/include/tinydir.h")
set(lib_sources "${lib_sources};${lib_headers};${bandit_headers};${tinydir_headers}")
set_source_files_properties(${lib_sources} PROPERTIES COMPILE_FLAGS " ${warnings_not_as_errors_flag} ")

# set OGDF_INSTALL and default visibility to hidden for shared libraries
if(BUILD_SHARED_LIBS)
  target_compile_definitions(OGDF PRIVATE OGDF_INSTALL)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(OGDF PRIVATE "-fvisibility=hidden" "-fvisibility-inlines-hidden")
  endif()
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

# autogen header variables for mallinfo2
include(CheckCXXSymbolExists)
CHECK_CXX_SYMBOL_EXISTS(mallinfo2 "malloc.h" OGDF_HAS_MALLINFO2)

# autogen header variables for SSE3
include(check-sse3)
if(has_sse3_intrin)
  set(OGDF_SSE3_EXTENSIONS <intrin.h>)
elseif(has_sse3_pmmintrin)
  set(OGDF_SSE3_EXTENSIONS <pmmintrin.h>)
else()
  message(STATUS "SSE3 could not be activated")
endif()

# autogen header variables for Linux-specific CPU_SET, etc.
include(check-cpu-macros)
if(has_linux_cpu_macros)
  set(OGDF_HAS_LINUX_CPU_MACROS 1)
endif()

if(CMAKE_COMPILER_IS_GNUCXX OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" AND NOT "${CMAKE_CXX_SIMULATE_ID}" STREQUAL "MSVC"))
  target_link_libraries(OGDF PUBLIC pthread)
endif()

# add stack trace settings
if(BACKWARD_HAS_DW)
  target_include_directories(OGDF SYSTEM PUBLIC ${LIBDW_INCLUDE_DIR})
elseif(BACKWARD_HAS_BFD)
  target_include_directories(OGDF SYSTEM PUBLIC ${LIBBFD_INCLUDE_DIR} ${LIBDL_INCLUDE_DIR})
elseif(BACKWARD_HAS_UNWIND)
  target_include_directories(OGDF SYSTEM PUBLIC ${LIBUNWIND_INCLUDE_DIR})
endif()
if(SHOW_STACKTRACE)
  if(BACKWARD_HAS_DW)
    target_link_libraries(OGDF PUBLIC ${LIBDW_LIBRARY})
  elseif(BACKWARD_HAS_BFD)
    target_link_libraries(OGDF PUBLIC ${LIBBFD_LIBRARY} ${LIBDL_LIBRARY})
  elseif(BACKWARD_HAS_UNWIND)
    target_link_libraries(OGDF PUBLIC ${LIBUNWIND_LIBRARY})
  endif()
endif()

# create autogen header(s)
function(create_autogen_header CURRENT_CONFIG)
  set_debug_mode(${CURRENT_CONFIG})
  if(CURRENT_CONFIG MATCHES Debug)
    configure_file("${module_dir}/config_autogen.h.in" "${PROJECT_BINARY_DIR}/include/ogdf-debug/ogdf/basic/internal/config_autogen.h")
  else()
    configure_file("${module_dir}/config_autogen.h.in" "${PROJECT_BINARY_DIR}/include/ogdf-release/ogdf/basic/internal/config_autogen.h")
  endif()
endfunction()

if(CMAKE_CONFIGURATION_TYPES)
  foreach(entry IN LISTS CMAKE_CONFIGURATION_TYPES)
    create_autogen_header(${entry})
  endforeach()
  set_debug_mode()
else()
  create_autogen_header(${CMAKE_BUILD_TYPE})
endif()

# installation
configure_file(cmake/ogdf-config.cmake "${PROJECT_BINARY_DIR}/ogdf-config.cmake" @ONLY)
install(TARGETS OGDF EXPORT OgdfTargets COMPONENT OGDF)
install(DIRECTORY "${PROJECT_BINARY_DIR}/include/" include/ogdf # copy everything *inside* the former dir and the latter dir itself
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  COMPONENT OGDFheaders
  FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.hpp"
    PATTERN "*.inc")
install(EXPORT OgdfTargets DESTINATION ${CMAKE_INSTALL_DATADIR}/ogdf)
install(FILES
        "${PROJECT_BINARY_DIR}/ogdf-config.cmake"
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindLibbfd.cmake"
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindLibdw.cmake"
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindLibunwind.cmake"
        DESTINATION ${CMAKE_INSTALL_DATADIR}/ogdf)
