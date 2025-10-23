find_path(LIBUNWIND_INCLUDE_DIR libunwind.h PATH_SUFFIXES libunwind)
find_library(LIBUNWIND_LIBRARY unwind)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Libunwind DEFAULT_MSG LIBUNWIND_LIBRARY LIBUNWIND_INCLUDE_DIR)
mark_as_advanced(LIBUNWIND_INCLUDE_DIR LIBUNWIND_LIBRARY)

add_library(Libunwind UNKNOWN IMPORTED)
set_target_properties(Libunwind PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${LIBUNWIND_INCLUDE_DIR}"
    IMPORTED_LOCATION "${LIBUNWIND_LIBRARY}"
)
