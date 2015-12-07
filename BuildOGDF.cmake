option(BUILD_OGDF_SHARED "Build OGDF Library as dynamic linked libraries" ON)
include(ExternalProject)
if (BUILD_OGDF_SHARED)
    ExternalProject_Add(ogdf
            GIT_REPOSITORY https://github.com/ogdf/ogdf.git
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/deps -DBUILD_SHARED_LIBS=ON -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            BUILD_COMMAND ${MAKE}
            INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/deps
            PREFIX ${CMAKE_CURRENT_BINARY_DIR}/deps
    )
else (BUILD_OGDF_SHARED)
    ExternalProject_Add(ogdf
            GIT_REPOSITORY https://github.com/ogdf/ogdf.git
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/deps -DBUILD_SHARED_LIBS=OFF -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            BUILD_COMMAND ${MAKE}
            INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/deps
            PREFIX ${CMAKE_CURRENT_BINARY_DIR}/deps
    )
endif (BUILD_OGDF_SHARED)
ExternalProject_Get_Property(ogdf install_dir)
set(OGDF_INCLUDE_DIR ${install_dir}/include)
if (BUILD_OGDF_SHARED)
    if (WIN32)
        set(COIN_LIBRARY_PATH ${install_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}COIN.dll)
    elseif(APPLE)
        set(COIN_LIBRARY_PATH ${install_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}COIN.dylib)
    else()
        set(COIN_LIBRARY_PATH ${install_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}COIN.so)
    endif()
    if (WIN32)
        set(OGDF_LIBRARY_PATH ${install_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}OGDF.dll)
    elseif(APPLE)
        set(OGDF_LIBRARY_PATH ${install_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}OGDF.dylib)
    else()
        set(OGDF_LIBRARY_PATH ${install_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}OGDF.so)
    endif()
else (BUILD_OGDF_SHARED)
    if (WIN32)
        set(COIN_LIBRARY_PATH ${install_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}COIN.lib)
    elseif(APPLE)
        set(COIN_LIBRARY_PATH ${install_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}COIN.a)
    else()
        set(COIN_LIBRARY_PATH ${install_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}COIN.a)
    endif()
    if (WIN32)
        set(OGDF_LIBRARY_PATH ${install_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}OGDF.lib)
    elseif(APPLE)
        set(OGDF_LIBRARY_PATH ${install_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}OGDF.a)
    else()
        set(OGDF_LIBRARY_PATH ${install_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}OGDF.a)
    endif()
endif (BUILD_OGDF_SHARED)
set(COIN_LIBRARY coin_lib)
set(OGDF_LIBRARY ogdf_lib)
add_library(${COIN_LIBRARY} UNKNOWN IMPORTED)
add_library(${OGDF_LIBRARY} UNKNOWN IMPORTED)
set_property(TARGET ${COIN_LIBRARY} PROPERTY IMPORTED_LOCATION ${COIN_LIBRARY_PATH})
set_property(TARGET ${OGDF_LIBRARY} PROPERTY IMPORTED_LOCATION ${OGDF_LIBRARY_PATH})
add_dependencies(${COIN_LIBRARY} ogdf)
add_dependencies(${OGDF_LIBRARY} ogdf)

MESSAGE( STATUS "OGDF_INCLUDE_DIR:               " ${OGDF_INCLUDE_DIR} )
MESSAGE( STATUS "COIN_LIBRARY_PATH:              " ${COIN_LIBRARY_PATH} )
MESSAGE( STATUS "OGDF_LIBRARY_PATH:              " ${OGDF_LIBRARY_PATH} )
MESSAGE( STATUS "COIN_LIBRARY:                   " ${COIN_LIBRARY} )
MESSAGE( STATUS "OGDF_LIBRARY:                   " ${OGDF_LIBRARY} )