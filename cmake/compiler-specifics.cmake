# CMake configuration that is related to special compilers

if(MSVC)
  # speed up builds by using parallel compilation & faster debug linking;
  # do not define min/max macros in windows.h
  add_definitions(/MP /DNOMINMAX)
  set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} /Debug:fastlink")
  string(REGEX REPLACE "/Z[iI7]" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Z7 /bigobj")

  # COIN has no DLL exports, must hence always be compiled as a static library
  set(COIN_LIBRARY_TYPE STATIC)

  if(CMAKE_CXX_COMPILER_LAUNCHER STREQUAL "ccache.exe")
    message(STATUS "Configuring MSVC build to use ccache.exe")

    # https://github.com/ccache/ccache/wiki/MS-Visual-Studio#usage-with-cmake
    find_program(ccache_exe ccache)
    if(ccache_exe)
      file(COPY_FILE
        ${ccache_exe} ${CMAKE_BINARY_DIR}/cl.exe
        ONLY_IF_DIFFERENT)

      # By default visual studio generators will use /Zi which is not compatible
      # with ccache for whatever reason.
      message(STATUS "Setting MSVC debug information format to 'Embedded'")
      set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<$<CONFIG:Debug,RelWithDebInfo>:Embedded>")

      set(CMAKE_VS_GLOBALS
        "CLToolExe=cl.exe"
        "CLToolPath=${CMAKE_BINARY_DIR}"
        "TrackFileAccess=false"
        "UseMultiToolTask=true"
        "DebugInformationFormat=OldStyle"
      )
    else()
      message(FATAL_ERROR "Could not find ccache.exe")
    endif()
  endif()
endif()

# use specified or native arch (ie, activate things like SSE)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang" AND NOT ${CMAKE_SYSTEM_PROCESSOR} MATCHES "^arm")
  set(OGDF_ARCH "native" CACHE STRING "Target CPU (micro)architecture passed to the compiler via `-march`.")
  mark_as_advanced(OGDF_ARCH)

  # cannot use add_definitions() here because it does not work with check-sse3.cmake
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=${OGDF_ARCH}")
else()
  unset(OGDF_ARCH CACHE)
endif()

# set default warning flags for OGDF and tests
set(available_default_warning_flags "")
set(available_default_warning_flags_debug "")
set(available_default_warning_flags_release "")
set(warnings_as_errors_flag "")
set(warnings_not_as_errors_flag "")

# Test for Windows first to allow ClangCl builds on Windows.
if(MSVC)
  set(available_default_warning_flags "/W3 /wd4018 /wd4068 /wd4101 /wd4244 /wd4250 /wd4267 /wd4373 /wd4800")
  # this has to be explained because MSVC is so cryptic:
  # /W3 sets the warning level of MSVC to 3 (all warnings except informational warnings),
  # /wd<code> disables the warning with the specific code,
  #     4018 = signed/unsigned mismatch
  #     4068 = unknown pragmas
  #     4101 = unused variable
  #     4244, 4267 = implicit conversion
  #     4250 = class inherits a member from another class via dominance
  #     4373 = behavior in old MSVC versions is different to C++ standard
  #     4800 = bool conversion from int or pointers
  set(warnings_as_errors_flag "/WX")
  set(warnings_not_as_errors_flag "/WX-")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  set(available_default_warning_flags "-Wall -Wextra -Wno-unused-parameter -Wno-unknown-pragmas \
      -Wno-error=sign-compare -Wno-error=conversion -Wno-error=strict-aliasing")
  if(CMAKE_CXX_COMPILER_ID MATCHES Clang)
    if(OGDF_INCLUDE_CGAL AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10.0)
      message(WARNING "Cannot compile with CGAL enabled with clang<10.0")
    endif()
    set(available_default_warning_flags "${available_default_warning_flags} -Wno-error=zero-length-array -Wno-error=uninitialized")
  else()
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.0)
      set(available_default_warning_flags "${available_default_warning_flags} -Wshadow")
    else()
      set(available_default_warning_flags "${available_default_warning_flags} -Wno-error=array-bounds")
    endif()

    # GCC seems to detect many false positives for the "stringop-overread" warning.
    # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=99578
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 11.0 AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 12.2)
      set(available_default_warning_flags "${available_default_warning_flags} -Wno-error=stringop-overread")
    endif()
    set(available_default_warning_flags "${available_default_warning_flags} -Wno-error=maybe-uninitialized -Wno-error=unused-but-set-parameter")
    set(available_default_warning_flags_release "-Wno-error=unused-but-set-variable -Wno-error=strict-overflow")
  endif()
  set(available_default_warning_flags_release "${available_default_warning_flags_release} -Wno-error=unused-variable")
  set(warnings_as_errors_flag "-Werror")
  set(warnings_not_as_errors_flag "-Wno-error")
endif()
