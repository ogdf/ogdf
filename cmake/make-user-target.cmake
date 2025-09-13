# Compilation of executables that do or do not use OGDF
function(make_some_target TARGET include_path)
  add_ogdf_extra_flags(${TARGET})
  target_compile_features(${TARGET} PUBLIC cxx_std_${CMAKE_CXX_STANDARD})
  target_include_directories(${TARGET} BEFORE PUBLIC ${include_path})
endfunction()

# Compilation of executables that use OGDF
function(make_user_target TARGET)
  make_some_target(${TARGET} ${PROJECT_BINARY_DIR}/include)
  target_link_libraries(${TARGET} OGDF)
  set_target_properties(${TARGET} PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})

  # link CGAL if enabled
  if(${OGDF_INCLUDE_CGAL})
    target_link_libraries(${TARGET} CGAL::CGAL CGAL::CGAL_Core)
  endif()
endfunction()
