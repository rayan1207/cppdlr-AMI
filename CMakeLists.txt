cmake_minimum_required(VERSION 3.20)
project(mini_ami VERSION 0.1 LANGUAGES CXX)

# ---- C++ standard ----
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---- Build type ----
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()







# ===========================
#   Dependencies
# ===========================
find_package(MPI REQUIRED)
find_package(Boost 1.89 CONFIG REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(cppdlr REQUIRED CONFIG)


# ===========================
#   Local AMI stack
# ===========================
set(AMIGRAPH_ROOT   "${CMAKE_CURRENT_SOURCE_DIR}/../libamigraph/install" CACHE PATH "amigraph install root")
set(AMICALC_ROOT    "${CMAKE_CURRENT_SOURCE_DIR}/../libamicalc/install"  CACHE PATH "amicalc install root")
set(AMI_ROOT        "${CMAKE_CURRENT_SOURCE_DIR}/../../../libami/install" CACHE PATH "ami install root")


set(AMIGRAPH_INCLUDE "${AMIGRAPH_ROOT}/include")
set(AMICALC_INCLUDE  "${AMICALC_ROOT}/include")
set(AMI_INCLUDE      "${AMI_ROOT}/include")


find_library(AMIGRAPH_LIB amigraph PATHS "${AMIGRAPH_ROOT}/lib" REQUIRED)
find_library(AMICALC_LIB  amicalc  PATHS "${AMICALC_ROOT}/lib"  REQUIRED)
find_library(AMI_LIB      ami      PATHS "${AMI_ROOT}/lib"      REQUIRED)


# ===========================
#   Executable
# ===========================
add_executable(${PROJECT_NAME}
  src/mini_main.cpp
  src/dlr_ami.hpp
  src/helper.cpp
  src/AmiDlr_func.cpp
  src/momenta_sampler.cpp
  src/ggm_mDLR.cpp
)

target_compile_definitions(${PROJECT_NAME} PRIVATE USE_MPI)

# Conda Torch headers must come first
target_include_directories(${PROJECT_NAME} BEFORE PRIVATE
  ${TORCH_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  ${MPI_INCLUDE_PATH}
  ${AMIGRAPH_INCLUDE}
 ${AMICALC_INCLUDE}
  ${AMI_INCLUDE}

)

# ===========================
#   Link libraries
# ===========================
target_link_libraries(${PROJECT_NAME} PRIVATE
  cppdlr::cppdlr_c
  MPI::MPI_CXX
  Boost::boost
  Eigen3::Eigen
  ${AMIGRAPH_LIB}
  ${AMICALC_LIB}
  ${AMI_LIB}
)

# ===========================
#   RPATH
# ===========================
if(DEFINED ENV{CPPDLR_ROOT})
  list(APPEND CMAKE_BUILD_RPATH "$ENV{CPPDLR_ROOT}/lib")
  list(APPEND CMAKE_INSTALL_RPATH "$ENV{CPPDLR_ROOT}/lib")
endif()

list(APPEND CMAKE_BUILD_RPATH
  "${AMIGRAPH_ROOT}/lib"
  "${AMICALC_ROOT}/lib"
  "${AMI_ROOT}/lib"

)
list(APPEND CMAKE_INSTALL_RPATH
  "${AMIGRAPH_ROOT}/lib"
  "${AMICALC_ROOT}/lib"
  "${AMI_ROOT}/lib"

)

set(CMAKE_MACOSX_RPATH ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# ===========================
#   Install
# ===========================
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)
