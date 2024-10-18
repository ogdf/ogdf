# This file is to be included by user code using find_package()

if("@OGDF_USE_ASSERT_EXCEPTIONS_WITH_STACK_TRACE@" STREQUAL "ON_LIBDW")
    find_package(Libdw)
elseif("@OGDF_USE_ASSERT_EXCEPTIONS_WITH_STACK_TRACE@" STREQUAL "ON_LIBBFD")
    find_package(Libbfd)
elseif("@OGDF_USE_ASSERT_EXCEPTIONS_WITH_STACK_TRACE@" STREQUAL "ON_LIBUNWIND")
    find_package(Libunwind)
endif()

include("${CMAKE_CURRENT_LIST_DIR}/CoinTargets.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/OgdfTargets.cmake")

set(OGDF_INCLUDE_DIRS $<TARGET_PROPERTY:OGDF,INTERFACE_INCLUDE_DIRECTORIES>)
