#
# CMake file to specify the build process, see:
# http://www.cmake.org/cmake/help/documentation.html
#

# Add COIN properties.
set(COIN_ENABLED TRUE CACHE BOOL "Enable COIN and ABACUS")
set(COIN_DEFAULT_SOLVER "CLP" CACHE STRING "Specifies the default solver used by COIN")
set_property(CACHE COIN_DEFAULT_SOLVER PROPERTY STRINGS 
	"CLP" 
	"SYM" 
	"CPX" 
	"GRB")
#TODO: NYI.
#option(COIN_EXTERNAL_SOLVERS "External solvers" "")
#option(COIN_SOLVER_INCLUDES "Additional include directories" "")

# Create config.h
if(COIN_ENABLED)
	set(COIN_ENABLED 1)
endif()
set(COIN_DEFAULT_SOLVER "COIN_OSI_${COIN_DEFAULT_SOLVER}")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/coinstuff/config.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/coin/config.h
	COPYONLY)

# Add COIN target.
source_dirs(COIN_SOURCES
	"src/coin"
	"include/coin")
list(REMOVE_ITEM COIN_SOURCES "src/coin/Osi/OsiGrbSolverInterface.cpp")
list(REMOVE_ITEM COIN_SOURCES "src/coin/Osi/OsiCpxSolverInterface.cpp")
add_definitions(
	"-DCLP_BUILD"
	"-DCOINUTILS_BUILD"
	"-DOSI_BUILD"
	"-DSYMPHONY_BUILD"
	"-D__OSI_CLP__"
	"-DCOMPILE_IN_CG"
	"-DCOMPILE_IN_CP"
	"-DCOMPILE_IN_LP"
	"-DCOMPILE_IN_TM"
	"-DUSE_CGL_CUTS")
if(NOT WIN32)
	#TODO: this isn't used by makeVC*Proj.py and breaks the build. intended?
	add_definitions("-DHAVE_CONFIG_H")
endif()
include_directories("include/coin")
if(COIN_ENABLED)
	add_library(coin STATIC ${COIN_SOURCES})
else()
	add_library(coin STATIC ${COIN_SOURCES} EXCLUDE_FROM_ALL)
endif()
