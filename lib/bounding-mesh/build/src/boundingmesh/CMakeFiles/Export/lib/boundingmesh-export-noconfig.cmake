#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "boundingmesh" for configuration ""
set_property(TARGET boundingmesh APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(boundingmesh PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libboundingmesh.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS boundingmesh )
list(APPEND _IMPORT_CHECK_FILES_FOR_boundingmesh "${_IMPORT_PREFIX}/lib/libboundingmesh.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
