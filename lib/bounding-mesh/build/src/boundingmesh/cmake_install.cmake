# Install script for directory: /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/boundingmesh.h"
    "/home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/ContractionUtils.h"
    "/home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/Decimator.h"
    "/home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/FileUtils.h"
    "/home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/Mesh.h"
    "/home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/MetricGenerator.h"
    "/home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/Primitives.h"
    "/home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/Segmenter.h"
    "/home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/SegmenterDownsampling.h"
    "/home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/SegmenterSimple.h"
    "/home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/Split.h"
    "/home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/VoxelSubset.h"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/ty/Physim/bounding-mesh-0.2/build/libboundingmesh.a")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE FILE FILES "/home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh/boundingmesh-config.cmake")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/boundingmesh-export.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/boundingmesh-export.cmake"
         "/home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh/CMakeFiles/Export/lib/boundingmesh-export.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/boundingmesh-export-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/boundingmesh-export.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE FILE FILES "/home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh/CMakeFiles/Export/lib/boundingmesh-export.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^()$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE FILE FILES "/home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh/CMakeFiles/Export/lib/boundingmesh-export-noconfig.cmake")
  endif()
endif()

