
Classes in this directory require building with CGAL. The following requirements must be met:

 * CMake flag `OGDF_INCLUDE_CGAL=ON`
 * CGAL >= 5.0
 * Boost >= 1.48
 * C++17 compliant compiler

Note that, for clang, version 10+ is required, because CGAL forces a `-frounding-math` flag which clang-9 does not understand.
