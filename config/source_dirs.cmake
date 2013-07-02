#
# CMake macros file, see:
# http://www.cmake.org/cmake/help/documentation.html
#

# This macro globs C++ sources and groups them into source groups to reflect the
# directory structure.
macro(source_dirs SOURCES)
    # Loop over all optional arguments.
    foreach(SOURCE_DIR ${ARGN})
        # Glob source files from directory.
        #message(STATUS "Globbing sources from ${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE_DIR}/")
        file(GLOB_RECURSE SOURCES_OF_DIR RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
            "${SOURCE_DIR}/*.cpp" "${SOURCE_DIR}/*.hh" "${SOURCE_DIR}/*.c"
            "${SOURCE_DIR}/*.hpp" "${SOURCE_DIR}/*.cc" "${SOURCE_DIR}/*.h" "${SOURCE_DIR}/*.inc")

        # Add to list of sources.
        set(${SOURCES} ${${SOURCES}} ${SOURCES_OF_DIR})

        # Add to list of group directories.
        foreach(SOURCE ${SOURCES_OF_DIR})
            get_filename_component(DIR_OF_SOURCE ${SOURCE} PATH)
            list(APPEND GROUP_DIRS ${DIR_OF_SOURCE})
        endforeach(SOURCE)
    endforeach(SOURCE_DIR)

    # Make group directories unique and loop them.
    list(REMOVE_DUPLICATES GROUP_DIRS)
    foreach(GROUP_DIR ${GROUP_DIRS})
        # Collect sources of group.
        set(SOURCES_OF_GROUP)
        foreach(SOURCE ${${SOURCES}})
            get_filename_component(SOURCE_DIR ${SOURCE} PATH)
            if(${SOURCE_DIR} STREQUAL ${GROUP_DIR})
                #message(STATUS "${SOURCE} -> ${GROUP_DIR}")
                set(SOURCES_OF_GROUP ${SOURCES_OF_GROUP} ${SOURCE})
            endif()
        endforeach(SOURCE)

        # Translate directory name to group name.
        string(REGEX REPLACE "^src(/[^/]*)?" "Source Files" GROUP_NAME "${GROUP_DIR}")
        string(REGEX REPLACE "^include(/[^/]*)?" "Header Files" GROUP_NAME "${GROUP_NAME}")
        string(REGEX REPLACE "^test" "Test Files" GROUP_NAME "${GROUP_NAME}")
        string(REGEX REPLACE "/" "\\\\" GROUP_NAME "${GROUP_NAME}")

        # Set source group.
        source_group(${GROUP_NAME} FILES ${SOURCES_OF_GROUP})
    endforeach(GROUP_DIR)
endmacro(source_dirs)
