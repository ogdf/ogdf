#
# CMake file to specify the build process, see:
# http://www.cmake.org/cmake/help/documentation.html
#

# Add COIN target.
file(GLOB_RECURSE COIN_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    "src/coin/*.cpp" "src/coin/*.h" "include/coin/*.h" "include/coin/*.hpp")
list(REMOVE_ITEM COIN_SOURCES "src/coin/Osi/OsiGrbSolverInterface.cpp")
list(REMOVE_ITEM COIN_SOURCES "src/coin/Osi/OsiCpxSolverInterface.cpp")
include_directories("include/coin")
add_library(coin STATIC ${COIN_SOURCES})
set(COIN_DEFINES
    "CLP_BUILD"
    "COINUTILS_BUILD"
    "OSI_BUILD"
    "SYMPHONY_BUILD"
    "__OSI_CLP__"
    "COMPILE_IN_CG"
    "COMPILE_IN_CP"
    "COMPILE_IN_LP"
    "COMPILE_IN_TM"
    "USE_CGL_CUTS")
if(WIN32)
    list(APPEND COIN_DEFINES
        "_CRT_SECURE_NO_WARNINGS"
        "_SCL_SECURE_NO_WARNINGS")
else()
    list(APPEND OGDF_DEFINES
        "HAVE_CONFIG_H")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/coinstuff/config.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/coin/config.h
        COPYONLY)
endif()
set_target_properties(coin PROPERTIES
    COMPILE_DEFINITIONS "${COIN_DEFINES}")
    