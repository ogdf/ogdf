#
# CMake file to specify the build process, see:
# http://www.cmake.org/cmake/help/documentation.html
#

# Add OGDF test target.
file(GLOB_RECURSE OGDF_TEST_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    "test/*.cpp" "test/*.h")
add_executable(ogdf-test ${OGDF_TEST_SOURCES})
set(OGDF_TEST_DEFINES "${OGDF_DEFINES}")
target_link_libraries(ogdf-test ogdf)
set_target_properties(ogdf-test PROPERTIES
    COMPILE_DEFINITIONS "${OGDF_TEST_DEFINES}")
