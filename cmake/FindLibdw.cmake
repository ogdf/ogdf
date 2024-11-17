find_path(LIBDW_INCLUDE_DIR "elfutils/libdw.h" "elfutils/libdwfl.h")
find_library(LIBDW_LIBRARY dw)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Libdw DEFAULT_MSG LIBDW_LIBRARY LIBDW_INCLUDE_DIR)
mark_as_advanced(LIBDW_INCLUDE_DIR LIBDW_LIBRARY)

add_library(Libdw UNKNOWN IMPORTED)
set_target_properties(Libdw PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${LIBDW_INCLUDE_DIR}"
        IMPORTED_LOCATION "${LIBDW_LIBRARY}"
)
