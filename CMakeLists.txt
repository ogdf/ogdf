cmake_minimum_required(VERSION 3.10...4.0)
project(OGDF-PROJECT LANGUAGES CXX VERSION "2023.09")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_DEBUG_POSTFIX "-debug")

set(module_dir "${PROJECT_SOURCE_DIR}/cmake")
list(INSERT CMAKE_MODULE_PATH 0 "${module_dir}" )

# build type configuration
set(BUILD_TYPES Debug Release)
if(CMAKE_CONFIGURATION_TYPES)
  set(MULTICONFIG_BUILD 1)
  set(CMAKE_CONFIGURATION_TYPES ${BUILD_TYPES} CACHE STRING "" FORCE)
else()
  set(MULTICONFIG_BUILD 0)
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build to be created." FORCE)
  endif()
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${BUILD_TYPES})
endif()

set(OGDF_INCLUDE_CGAL OFF CACHE BOOL "Indicates whether components that require CGAL ({src,include}/ogdf/geometric) should be built. Requires OpenMP")

include(compiler-specifics)
include(group-files)
include(coin)
include(ogdf)

include(make-user-target)
include(tests)
include(examples)

include(doc)

add_custom_target(build-all DEPENDS OGDF tests examples)

# target documentation
message(STATUS "The default target builds OGDF (and dependencies like COIN).")
if(DOXYGEN_FOUND AND DOC_INSTALL)
  message(STATUS "Because DOC_INSTALL is set, it also generates the OGDF documentation (in-source).")
endif()
message(STATUS "")
message(STATUS "Other important targets:")
if(DOXYGEN_FOUND)
  if (DOC_INSTALL)
    message(STATUS "        OGDF: build OGDF library only")
  endif()
  message(STATUS "         doc: build doxygen documentation (in-source)")
endif()
message(STATUS "       tests: build tests")
message(STATUS "    examples: build examples")
message(STATUS "   build-all: build OGDF, tests, examples")
