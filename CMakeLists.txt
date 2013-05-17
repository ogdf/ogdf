#
# CMake listfile to specify the build process, see:
# http://www.cmake.org/cmake/help/documentation.html
#
project(OGDF)

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

# We can't use add_subdirectory() here, because of OGDF's special directory
# structure. Thus we'll have to include projects (order is important!).
include(${CMAKE_CURRENT_SOURCE_DIR}/config/coin.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/config/ogdf.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/config/ogdf-test.cmake)
