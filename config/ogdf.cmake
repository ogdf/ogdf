#
# CMake file to specify the build process, see:
# http://www.cmake.org/cmake/help/documentation.html
#

# Add OGDF options.
option(OGDF_DLL "Build OGDF as shared library" FALSE)

# Add OGDF target.
source_dirs(OGDF_SOURCES
    "src/ogdf"
    "include/ogdf")
include_directories("include")
if(OGDF_DLL)
    add_library(ogdf SHARED ${OGDF_SOURCES})
    set(OGDF_DEFINES
        "${COIN_DEFINES}"
        "OGDF_DLL")
else()
    add_library(ogdf STATIC ${OGDF_SOURCES})
    set(OGDF_DEFINES
        "${COIN_DEFINES}")
endif()
if(WIN32)
    target_link_libraries(ogdf "Psapi.lib")
endif()
target_link_libraries(ogdf coin)
set_target_properties(ogdf PROPERTIES
    COMPILE_DEFINITIONS "${OGDF_DEFINES}")
