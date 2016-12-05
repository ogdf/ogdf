# Tests-related CMake configuration

# cache configuration
option(OGDF_SEPARATE_TESTS "Whether to build separate test executables (used for continuous integration)" OFF)

function(make_test_executable TARGET)
  target_include_directories(${TARGET} BEFORE PUBLIC test/include)
  make_user_executable(${TARGET})
endfunction()

file(GLOB_RECURSE TEST_SOURCES test/src/*.cpp)
group_files(TEST_SOURCES "test")
if(OGDF_SEPARATE_TESTS)
  SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/test/bin")
  foreach(SOURCE_FILE ${TEST_SOURCES})
    get_filename_component(TARGET ${SOURCE_FILE} NAME_WE)
    if(NOT ${TARGET} STREQUAL "main")
      add_executable(${TARGET} ${SOURCE_FILE} "test/src/main.cpp")
      make_test_executable(${TARGET})
    endif()
  endforeach()
else()
  add_executable(test-ogdf ${TEST_SOURCES})
  make_test_executable(test-ogdf)
endif()
