cmake_minimum_required(VERSION 3.18)
project(mini_ami VERSION 0.1 LANGUAGES CXX)

# --------------------
# C++ Standard
# --------------------
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --------------------
# Build Type & Flags
# --------------------
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()
add_compile_options(-O3)
add_compile_definitions(USE_MPI)

# --------------------
# External Dependencies
# --------------------
find_package(Boost 1.71 REQUIRED COMPONENTS filesystem system)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# --------------------
# Executable
# --------------------
add_executable(${PROJECT_NAME}
  src/mini_main.cpp
  src/dlr_ami.hpp
  src/dlr_helper.cpp
  #src/sparse_sampler_sigma.cpp
  #src/helper.cpp
  #src/molecule_solver.cpp
  #src/energy.cpp
)
target_compile_options(${PROJECT_NAME} PRIVATE
  -O3
  -march=native
  -funroll-loops
  -ffast-math
)
# --------------------
# Include Directories
# --------------------
target_include_directories(${PROJECT_NAME} PRIVATE
  ${Boost_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../../../libami/install/include
  ${CMAKE_CURRENT_SOURCE_DIR}/../libamicalc/install/include
  ${CMAKE_CURRENT_SOURCE_DIR}/../libamigraph/install/include
  /mnt/c/Users/Rayan/Research/leblanc_codes/AMI/cppdlr/include
)

# --------------------
# Link Libraries
# --------------------
target_link_libraries(${PROJECT_NAME} PRIVATE
  # Ami stack (in dependency order)
  ${CMAKE_CURRENT_SOURCE_DIR}/../libamigraph/install/lib/libamigraph.so
  ${CMAKE_CURRENT_SOURCE_DIR}/../libamicalc/install/lib/libamicalc.so
  ${CMAKE_CURRENT_SOURCE_DIR}/../../../libami/install/lib/libami.so

  # cppdlr C API static archives
  /mnt/c/Users/Rayan/Research/leblanc_codes/AMI/cppdlr/lib/libcppdlr_c.a
  /mnt/c/Users/Rayan/Research/leblanc_codes/AMI/cppdlr/lib/libnda_c.a
  /mnt/c/Users/Rayan/Research/leblanc_codes/AMI/cppdlr/lib/libh5_c.a
  /mnt/c/Users/Rayan/Research/leblanc_codes/AMI/cppdlr/lib/libfmt.a

  # BLAS and LAPACK
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}

  # Standard and Boost
  m
  stdc++fs
  Boost::filesystem
  Boost::system
)

# --------------------
# Installation
# --------------------
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)
