#
# CMake file to specify the build process, see:
# http://www.cmake.org/cmake/help/documentation.html
#

# Add OGDF properties.
set(OGDF_DLL FALSE CACHE BOOL "Build OGDF as shared library")
set(OGDF_MEMORY_MANAGER "" CACHE STRING "Specifies the memory manager used by OGDF")
set_property(CACHE OGDF_MEMORY_MANAGER PROPERTY STRINGS 
	"" 
	"OGDF_MEMORY_POOL_TS" 
	"OGDF_MEMORY_POOL_NTS" 
	"OGDF_MEMORY_MALLOC_TS")
	
# Create config_autogen.h
if(NOT OGDF_MEMORY_MANAGER STREQUAL "")
	set(OGDF_MEMORY_MANAGER "#define ${OGDF_MEMORY_MANAGER}")
endif()
if(OGDF_DLL)
	set(OGDF_DLL 1)
else()
	set(OGDF_DLL 0)
endif()
if(NOT WIN32)
	set(WIN32 0)
endif()
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/ogdfstuff/config_autogen.h.in
	${CMAKE_CURRENT_SOURCE_DIR}/include/ogdf/internal/config_autogen.h
	@ONLY)

# Add OGDF target.
source_dirs(OGDF_SOURCES
	"src/ogdf"
	"include/ogdf")
foreach(SOURCE ${OGDF_SOURCES})
	get_filename_component(NAME_OF_SOURCE ${SOURCE} NAME)
	# Filter ignored files
	if(${NAME_OF_SOURCE} MATCHES "^[_\\.]" OR ${SOURCE} MATCHES "ogdf/legacy")
		list(REMOVE_ITEM OGDF_SOURCES ${SOURCE})
	endif()
endforeach(SOURCE)
include_directories("include")
if(OGDF_DLL)
	add_library(ogdf SHARED ${OGDF_SOURCES})
else()
	add_library(ogdf STATIC ${OGDF_SOURCES})
endif()
if(WIN32)
	target_link_libraries(ogdf "Psapi.lib")
endif()
if(COIN_ENABLED)
	target_link_libraries(ogdf coin)
endif()
