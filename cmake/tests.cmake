# Tests-related CMake configuration

# cache configuration
option(OGDF_SEPARATE_TESTS "Whether to build separate test executables (used for continuous integration)" OFF)

function(make_test_executable TARGET)
  target_include_directories(${TARGET} BEFORE PUBLIC test/include)
  make_user_executable(${TARGET})
  addOgdfExtraFlags(${TARGET})
endfunction()

file(GLOB_RECURSE TEST_SOURCES test/src/*.cpp)
group_files(TEST_SOURCES "test")
if(OGDF_SEPARATE_TESTS)
  add_custom_target(tests)
  SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/test/bin")
  foreach(SOURCE_FILE ${TEST_SOURCES})
    get_filename_component(TARGET ${SOURCE_FILE} NAME_WE)
    if(NOT ${TARGET} STREQUAL "main")
      add_executable(test-${TARGET} EXCLUDE_FROM_ALL ${SOURCE_FILE} "test/src/main.cpp")
      add_dependencies(tests test-${TARGET})
      make_test_executable(test-${TARGET})
    endif()
  endforeach()
else()
  add_executable(tests EXCLUDE_FROM_ALL ${TEST_SOURCES})
  make_test_executable(tests)
endif()
