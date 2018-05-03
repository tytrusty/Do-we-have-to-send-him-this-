cmake_minimum_required(VERSION 3.1)
project(heat)
message(STATUS "SOURCE DIR IS: " ${PROJECT_SOURCE_DIR})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(CGAL)
if(CGAL_FOUND)
    add_definitions(-DCGAL_AVAILABLE -DCGAL_DISABLE_ROUNDING_MATH_CHECK)
    include_directories(${CGAL_INCLUDE_DIRS})
    link_directories(${CGAL_LIBRARIES_DIR})
    message(STATUS ${CGAL_INCLUDE_DIRS})
    message(STATUS ${CGAL_LIBRARIES_DIR})
    set(CGAL_LIBRARIES "CGAL;CGAL_Core;gmp;mpfr") # CGAL
    set(HAVE_CGAL True)
else()
    set(HAVE_CGAL False)
    message(STATUS "Configure package without CGAL!")
endif()


# libigl
option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR       "Use AntTweakBar"    OFF)
option(LIBIGL_WITH_CGAL              "Use CGAL"           OFF)
option(LIBIGL_WITH_COMISO            "Use CoMiso"         OFF)
option(LIBIGL_WITH_CORK              "Use Cork"           OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"         OFF)
option(LIBIGL_WITH_LIM               "Use LIM"            OFF)
option(LIBIGL_WITH_MATLAB            "Use Matlab"         OFF)
option(LIBIGL_WITH_MOSEK             "Use MOSEK"          OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          ON)
option(LIBIGL_WITH_PNG               "Use PNG"            OFF)
option(LIBIGL_WITH_PYTHON            "Use Python"         OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"       OFF)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_XML               "Use XML"            OFF)

find_package(LIBIGL REQUIRED QUIET)

include_directories(${PROJECT_SOURCE_DIR}/lib/bounding-mesh/src/boundingmesh )

# Add your project files
file(GLOB SRCFILES *.cpp)
file(GLOB SRCBOUND ${PROJECT_SOURCE_DIR}/lib/bounding-mesh/src/boundingmesh/*.cpp )
add_executable(${PROJECT_NAME}_bin ${SRCFILES} ${SRCBOUND})
MESSAGE(${CMAKE_SOURCE_DIR})
target_include_directories(${PROJECT_NAME}_bin PUBLIC "" ${CGAL_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::opengl_glfw igl::opengl_glfw_imgui ${CGAL_LIBRARIES})


