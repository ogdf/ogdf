# Compilation of examples

option(BUILD_EXAMPLES "Whether to build the examples (in-source) for the documentation" OFF)

if(BUILD_EXAMPLES)
  file(GLOB_RECURSE example_sources doc/examples/*.cpp)
  foreach(source ${example_sources})
    get_filename_component(target ${source} NAME_WE)
    get_filename_component(targetdir ${source} DIRECTORY)
    add_executable(ex-${target} ${source})
    set_property(TARGET ex-${target} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${targetdir}")
    make_user_executable(ex-${target})
  endforeach()
endif()
