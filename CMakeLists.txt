cmake_minimum_required( VERSION 3.10 )

project( simulation_setup )

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#Optimization parameters
set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_CXX_FLAGS "-O2")

# Armadillo linking
find_package(Armadillo REQUIRED)

# For gnuplot
# find_package(Boost COMPONENTS iostreams filesystem system REQUIRED)
find_package(Boost REQUIRED COMPONENTS system iostreams filesystem)

# GMSH linking
find_library(GMSH_LIB gmsh)
if(NOT GMSH_LIB)
  message(FATAL_ERROR "Could not find libgmsh")
endif()

find_path(GMSH_INC gmsh.h)
if(NOT GMSH_INC)
  message(FATAL_ERROR "Could not find gmsh.h")
endif()

# OpenMP linking
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# FLTK linking
find_package(FLTK REQUIRED)

# OPENGL linking
find_package(OpenGL REQUIRED)

# CGAL 
find_package(CGAL)

# Include directories
include_directories(
    ./include )

# Set libraries
set(ALL_LIBS ${ARMADILLO_LIBRARIES} ${Boost_LIBRARIES} ${FLTK_LIBRARIES}
  ${OPENGL_LIBRARIES} ${GMSH_LIB} CGAL::CGAL Boost::iostreams)

set(SOURCES
    ./src/mesh2D.cpp
    ./src/structured2D_mesh.cpp
    ./src/rpim2D.cpp
    ./src/boundary_conditions.cpp
    ./src/pointcloud_rpim.cpp
    ./src/support_domain.cpp
    ./src/geometry_model.cpp
    ./src/shape_function.cpp
    ./src/strain.cpp
    ./src/material.cpp
    ./src/gq_triangle_rpim.cpp
    ./src/gq_line_rpim.cpp
    ./src/loading_conditions.cpp
    ./src/point_loads.cpp
    ./src/kd_trees.cpp
    )
    

# target
add_executable(main main.cpp ${SOURCES})

target_include_directories(main PRIVATE ${ARMADILLO_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS} ${GMSH_INC} ${FLTK_INCLUDE_DIRS} )
    
target_link_libraries(main ${ALL_LIBS})