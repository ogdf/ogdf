# Documentation-related CMake configuration

find_package(Doxygen)
set(DOC_DIR "${CMAKE_SOURCE_DIR}/doc")
add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${DOC_DIR}/ogdf-doxygen.cfg WORKING_DIRECTORY ${DOC_DIR})
