# CMake configuration that is related to special compilers

if(MSVC)
  # a hack to make recursive DFS work on big graphs with MSVC
  add_definitions(/bigobj)
  set(CMAKE_CXX_STACK_SIZE "10000000")

  # speed up builds by using parallel compilation & faster debug linking
  add_definitions(/MP)
  set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} /Debug:fastlink")
  string(REGEX REPLACE "/Z[iI7]" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Z7")

  # COIN has no DLL exports, must hence always be compiled as a static library
  set(COIN_LIBRARY_TYPE STATIC)
endif()

# use native arch (ie, activate things like SSE)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  # cannot use add_definitions() here because it does not work with check-sse3.cmake
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

# set default warning flags for OGDF and tests
set(available_default_warning_flags "")
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  set(available_default_warning_flags "-Wall -Wextra -Wno-unused-parameter -Werror -Wno-error=deprecated-declarations -Wno-unknown-pragmas -Wno-error=sign-compare -Wno-error=conversion")
  if(CMAKE_CXX_COMPILER_ID MATCHES Clang)
    set(available_default_warning_flags "${available_default_warning_flags} -Wno-error=zero-length-array -Wno-error=uninitialized")
  else()
    set(available_default_warning_flags "${available_default_warning_flags} -Wno-error=maybe-uninitialized")
    if(NOT CMAKE_BUILD_TYPE STREQUAL Debug)
      set(available_default_warning_flags "${available_default_warning_flags} -Wno-error=unused-but-set-variable")
    endif()
  endif()
  if(NOT CMAKE_BUILD_TYPE STREQUAL Debug)
    set(available_default_warning_flags "${available_default_warning_flags} -Wno-error=unused-variable")
  endif()
elseif(MSVC)
  set(available_default_warning_flags "/W3 /WX /wd4018 /wd4068 /wd4101 /wd4244 /wd4250 /wd4267 /wd4373 /wd4800 /wd4996")
  # this has to be explained because MSVC is so cryptic:
  # /W3 sets the warning level of MSVC to 3 (all warnings except informational warnings),
  # /WX tells MSVC to treat all warnings as errors,
  # /wd<code> disables the warning with the specific code,
  #     4018 = signed/unsigned mismatch
  #     4068 = unknown pragmas
  #     4101 = unused variable
  #     4244, 4267 = implicit conversion
  #     4250 = class inherits a member from another class via dominance
  #     4373 = behavior in old MSVC versions is different to C++ standard
  #     4800 = bool conversion from int or pointers
  #     4996 = deprecated declaration
endif()
