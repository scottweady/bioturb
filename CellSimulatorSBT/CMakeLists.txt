cmake_minimum_required(VERSION 3.17)
project(CellSimulator)

# global settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
option(EIGEN_USE_BLAS "Use BLAS as the backend for various Eigen calls" ON)

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
# this does not change the compiler
# to use mpi compiler, do cmake -D CMAKE_CXX_COMPILER=mpicxx
# linking executable with mpi
# target_link_libraries(exe PRIVATE OpenMP::OpenMP_CXX)
# target_link_libraries(exe PRIVATE MPI::MPI_CXX)

# compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-sign-compare")
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using intel, disable offloading
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qno-offload")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
endif()

if(EIGEN_USE_BLAS)
  add_compile_definitions("EIGEN_USE_BLAS")
  message("Using BLAS backend for Eigen")
endif()

add_compile_options("$<$<CONFIG:RELEASE>:-g>")
add_compile_options("$<$<CONFIG:RELEASE>:-O3>")
add_compile_options(-DSCTL_QUAD_T=__float128)

message(STATUS "variable is: " ${CMAKE_CXX_COMPILER_ID})

# Get VTK
find_package(VTK REQUIRED COMPONENTS IOXML HINTS ${SFTPATH}/lib/cmake
                                     $ENV{VTK_BASE} $ENV{VTK_ROOT})
message("Found VTK at ${VTK_DIR}")

# external libraries, header only
find_package(Eigen3 REQUIRED HINTS ${SFTPATH}/share/eigen3/cmake
             $ENV{EIGEN_BASE} $ENV{EIGEN_ROOT})
message("Found Eigen3 at ${Eigen3_DIR}")
find_package(Boost REQUIRED HINTS ${SFTPATH}/lib/cmake $ENV{BOOST_BASE}
             $ENV{BOOST_ROOT})
message("Found Boost at ${Boost_DIR}")

# Get Trilinos as one entity
find_package(Trilinos REQUIRED HINTS ${SFTPATH}/lib/cmake $ENV{TRILINOS_BASE}
             $ENV{TRILINOS_ROOT})

# clean up the too-long list
list(REMOVE_DUPLICATES Trilinos_INCLUDE_DIRS)
list(REMOVE_DUPLICATES Trilinos_TPL_INCLUDE_DIRS)

list(REVERSE Trilinos_LIBRARIES)
list(REMOVE_DUPLICATES Trilinos_LIBRARIES)
list(REVERSE Trilinos_LIBRARIES)

list(REVERSE Trilinos_TPL_LIBRARIES)
list(REMOVE_DUPLICATES Trilinos_TPL_LIBRARIES)
list(REVERSE Trilinos_TPL_LIBRARIES)

# add full path to Trilinos_LIBRARIES
foreach(X IN LISTS ${Trilinos_LIBRARIES})
    set(X ${Trilinos_LIBRARY_DIRS}/${Trilinos_LIBRARIES})
endforeach()

# Echo trilinos build info just for fun
MESSAGE("\nFound Trilinos!  Here are the details: ")
MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
MESSAGE("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
MESSAGE("End of Trilinos details\n")

include(FetchContent)

FetchContent_Declare(
  trng4
  GIT_REPOSITORY https://github.com/rabauke/trng4.git
  #GIT_SHALLOW ON
  GIT_TAG 1ef58a7)
FetchContent_Populate(trng4)
message(${trng4_SOURCE_DIR} ${trng4_BINARY_DIR})
add_subdirectory(${trng4_SOURCE_DIR} ${trng4_BINARY_DIR} EXCLUDE_FROM_ALL)

FetchContent_Declare(
  yamlcpp
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp
  GIT_SHALLOW ON
  GIT_TAG yaml-cpp-0.7.0)
FetchContent_Populate(yamlcpp)
message(${yamlcpp_SOURCE_DIR} ${yamlcpp_BINARY_DIR})
add_subdirectory(${yamlcpp_SOURCE_DIR} ${yamlcpp_BINARY_DIR} EXCLUDE_FROM_ALL)

message(
  STATUS
    "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS"
)
message(STATUS "SFTPATH: ${SFTPATH}")

enable_testing()

# internal submodules
set(SIMTOOLBOX ${CMAKE_SOURCE_DIR}/SimToolbox)
set(SRC ${CMAKE_SOURCE_DIR}/SRC)

add_subdirectory(STKFMM)

add_executable(
  CellSimulator.X
  ${SRC}/Config.cpp
  ${SRC}/main.cpp
  ${SRC}/GrowthTable.cpp
  ${SRC}/HydroTable.cpp
  ${SIMTOOLBOX}/Boundary/Boundary.cpp
  ${SIMTOOLBOX}/Constraint/BCQPSolver.cpp
  ${SIMTOOLBOX}/Constraint/ConstraintCollector.cpp
  ${SIMTOOLBOX}/Constraint/ConstraintOperator.cpp
  ${SIMTOOLBOX}/Constraint/ConstraintSolver.cpp
  ${SIMTOOLBOX}/Sylinder/SylinderConfig.cpp
  ${SIMTOOLBOX}/Trilinos/TpetraUtil.cpp
  ${SIMTOOLBOX}/Util/Base64.cpp)

target_compile_options(CellSimulator.X PRIVATE ${OpenMP_CXX_FLAGS})

target_compile_definitions(
  CellSimulator.X
  PRIVATE # FDPS Config
          PARTICLE_SIMULATOR_MPI_PARALLEL PARTICLE_SIMULATOR_THREAD_PARALLEL)

target_include_directories(
  CellSimulator.X
  PRIVATE ${CMAKE_SOURCE_DIR}
          ${SRC}
          ${SIMTOOLBOX}
          ${Trilinos_INCLUDE_DIRS})

target_link_libraries(
  CellSimulator.X
  PRIVATE ${Trilinos_LIBRARIES}
          ${Trilinos_TPL_LIBRARIES}
          ${STKFMM_LIBRARIES}
          ${PVFMM_LIB_DIR}/${PVFMM_STATIC_LIB}
          ${PVFMM_DEP_LIB}
          STKFMM_STATIC
          yaml-cpp
          VTK::IOXML
          OpenMP::OpenMP_CXX
          MPI::MPI_CXX
          trng4
          Eigen3::Eigen)

# for generated gitversion.h
add_custom_command(
  TARGET CellSimulator.X
  PRE_BUILD
  COMMAND
    ${CMAKE_COMMAND} -Dlocal_dir="${CMAKE_SOURCE_DIR}"
    -Doutput_dir="${CMAKE_BINARY_DIR}" -P
    "${CMAKE_SOURCE_DIR}/cmake/gitversion.cmake")

install(DIRECTORY ${CMAKE_SOURCE_DIR}/Run/postprocess
        DESTINATION ${CMAKE_INSTALL_PREFIX})
install(DIRECTORY ${CMAKE_SOURCE_DIR}/Run/result
        DESTINATION ${CMAKE_INSTALL_PREFIX})
install(DIRECTORY ${CMAKE_SOURCE_DIR}/Run/scripts
        DESTINATION ${CMAKE_INSTALL_PREFIX})
install(DIRECTORY ${CMAKE_SOURCE_DIR}/Example
        DESTINATION ${CMAKE_INSTALL_PREFIX})
install(FILES ${CMAKE_SOURCE_DIR}/Run/jobsub.slurm
        DESTINATION ${CMAKE_INSTALL_PREFIX})
install(FILES ${CMAKE_SOURCE_DIR}/Run/Submit.sh
        DESTINATION ${CMAKE_INSTALL_PREFIX})
install(FILES ${CMAKE_BINARY_DIR}/gitversion.txt
        DESTINATION ${CMAKE_INSTALL_PREFIX})
install(TARGETS CellSimulator.X DESTINATION ${CMAKE_INSTALL_PREFIX})
