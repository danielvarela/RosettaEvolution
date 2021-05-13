cmake_minimum_required(VERSION 3.9)
project(MiniFolding)

set (CMAKE_CXX_STANDARD 11)

# Rosetta flags for compatibility
############ Build Type, Mode, Link Rules, Compiler, Extras
SET(EXTRAS cxx11)
SET(MODE debug)
SET(LINK_RULES shared)
SET(CMAKE_BUILD_TYPE Debug)

SET(CC gcc)
SET(CXX g++)

ADD_DEFINITIONS(-DCXX11)
ADD_DEFINITIONS(-DPTR_STD)
REMOVE_DEFINITIONS(-DPTR_MODERN)
REMOVE_DEFINITIONS(-DPTR_BOOST)
SET(COMPILE_FLAGS "${COMPILE_FLAGS} -Wno-error -Wno-deprecated-declarations -O3 -lz")

# output binary directory 
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin/)
SET(CMAKE_INSTALL_PREFIX ../bin/)

# extra packages
find_package(OpenMP)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


find_package(Boost 1.61.0 COMPONENTS serialization mpi REQUIRED )
if(Boost_FOUND)
   message(STATUS "Boost target found: " ${BOOST_ROOT})
   message("Boost_FOUND" ${BOOST_FOUND})
   message("Boost_INCLUDE_DIRS" ${BOOST_INCLUDE_DIRS})  
   message("Boost_LIBRARY_DIRS" ${BOOST_LIBRARY_DIRS})
   message("Boost_LIBRARIES" ${BOOST_LIBRARIES})
   message("Boost found")
endif()

#SET(CMAKE_C_COMPILER mpicc)
#SET(CMAKE_CXX_COMPILER mpicxx)

REMOVE_DEFINITIONS(-DPTR_MODERN)
REMOVE_DEFINITIONS(-DPTR_BOOST)

ADD_DEFINITIONS(-DCXX11)
ADD_DEFINITIONS(-DPTR_STD)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin/)
SET(CMAKE_INSTALL_PREFIX ../bin/)


find_package(MPI)


# set(INCLUDE_ROSETTA_LIST 
#   /home/ulc/co/dvm/intel/ism/rm/rosetta2018/main/source/src
#   /home/ulc/co/dvm/intel/ism/rm/rosetta2018/main/source/external/include
#   #/opt/cesga/easybuild-cesga/software/Core/icc/2018.5.274/include/
#   /home/ulc/co/dvm/intel/ism/rm/rosetta2018/main/source/src/platform/linux/
#   /home/ulc/co/dvm/intel/ism/rm/rosetta2018/main/source/external/dbio
#   /home/ulc/co/dvm/intel/ism/rm/rosetta2018/main/source/external/libxml2/include
#   /usr/include
#   /usr/local/include)
	
# set(LINK_ROSETTA_LIST 
#   "/home/ulc/co/dvm/intel/ism/rm/rosetta2018/main/source/external/lib"
#   "/home/ulc/co/dvm/intel/ism/rm/rosetta2018/main/source/cmake/build_cxx11/"
#   "/home/ulc/co/dvm/intel/ism/rm/rosetta2018/main/source/src"
#   "/home/ulc/co/dvm/intel/ism/rm/rosetta2018/main/source/external"
#   "/usr/local/lib"
#   "/home/ulc/co/dvm/intel/ism/rm/rosetta2018/main/source/external/lib" 
#   "/home/ulc/co/dvm/intel/ism/rm/rosetta2018/main/source/src"
#   "/home/ulc/co/dvm/intel/ism/rm/rosetta2018/main/source/external"
#   "/opt/cesga/easybuild-cesga/software/Core/icc/2018.5.274/lib/"
#   )

set(INCLUDE_ROSETTA_LIST /home/dvarela/Code/rosetta/main/source/src
  /home/dvarela/Code/rosetta/main/source/external/include
  /home/dvarela/Code/rosetta/main/source/src/platform/linux/64/gcc/7
  /home/dvarela/Code/rosetta/main/source/src/platform/linux/64/gcc
  /home/dvarela/Code/rosetta/main/source/src/platform/linux/64
  /home/dvarela/Code/rosetta/main/source/src/platform/linux
#  /home/dvarela/Code/rosetta/main/source/external/boost_1_55_0
  /home/dvarela/Code/rosetta/main/source/external/libxml2/include
  /home/dvarela/Code/rosetta/main/source/external
  /home/dvarela/Code/rosetta/main/source/external/dbio
  /usr/include
  /usr/local/include)
	
set(LINK_ROSETTA_LIST "/home/dvarela/Code/rosetta/main/source/external/lib"
  "/home/dvarela/Code/rosetta/main/source/cmake/build_cxx11/"
  "/home/dvarela/Code/rosetta/main/source/src"
  "/home/dvarela/Code/rosetta/main/source/build/external/debug/linux/4.4/64/x86/gcc/7/static"
  "/home/dvarela/Code/rosetta/main/source/external"
  "/usr/lib/"
  "/usr/local/lib")



# My app source files
file(GLOB_RECURSE sources src/main.cpp src/Controller/*.cpp src/Algorithm/*.cpp src/Extra/*.cpp src/Movers/*.cpp src/MpiFiles/*.cpp ) 

set(_definitions MPI_ENABLED=true USE_CRYO_EM=false)
include_directories(${INCLUDE_ROSETTA_LIST})
include_directories(${MPI_INCLUDE_PATH})
include_directories(${BOOST_INCLUDE_DIRS})
link_directories(${LINK_ROSETTA_LIST})
add_executable(mpi_app ${sources})
target_compile_definitions(mpi_app PRIVATE USING_MPI=1)
set_target_properties(mpi_app PROPERTIES COMPILE_DEFINITIONS "${_definitions}")
target_compile_options(mpi_app PRIVATE ${Open_MP_CXX_FLAGS})
target_link_libraries(mpi_app devel protocols.8 protocols.7 protocols_e.6 protocols_d.6 protocols_c.6 protocols_b.6 protocols_a.6 protocols.1 protocols_h.5 protocols_g.5 protocols_f.5 protocols_e.5 protocols_d.5 protocols_c.5 protocols_b.5 protocols_a.5 protocols.4  protocols_a.2 protocols_b.2 protocols.3 core.5 core.4 core.3 core.2 core.1 basic numeric utility ObjexxFCL z cppdb xml2 sqlite3 cifparse mpi boost_date_time boost_serialization boost_mpi ${Open_MP_CXX_FLAGS} ${MPI_LIBRARIES})

set(_definitions MPI_ENABLED=true USE_CRYO_EM=false)
include_directories(${INCLUDE_ROSETTA_LIST})
include_directories(${MPI_INCLUDE_PATH})
include_directories(${BOOST_INCLUDE_DIRS})
link_directories(${LINK_ROSETTA_LIST})
add_executable(temp ${sources})
target_compile_definitions(temp PRIVATE USING_MPI=1)
set_target_properties(temp PROPERTIES COMPILE_DEFINITIONS "${_definitions}")
target_compile_options(temp PRIVATE ${Open_MP_CXX_FLAGS})
target_link_libraries(temp devel protocols.8 protocols.7 protocols_e.6 protocols_d.6 protocols_c.6 protocols_b.6 protocols_a.6 protocols.1 protocols_h.5 protocols_g.5 protocols_f.5 protocols_e.5 protocols_d.5 protocols_c.5 protocols_b.5 protocols_a.5 protocols.4  protocols_a.2 protocols_b.2 protocols.3 core.5 core.4 core.3 core.2 core.1 basic numeric utility ObjexxFCL z cppdb xml2 sqlite3 cifparse mpi boost_date_time boost_serialization boost_mpi ${Open_MP_CXX_FLAGS} ${MPI_LIBRARIES})


# My app source files
file(GLOB_RECURSE sources src/test_main.cpp src/Controller/*.cpp src/Algorithm/*.cpp src/Extra/*.cpp src/Movers/*.cpp src/MpiFiles/*.cpp src/Test/*.cpp ) 

set(_definitions MPI_ENABLED=true USE_CRYO_EM=false)
include_directories(${INCLUDE_ROSETTA_LIST})
include_directories(${MPI_INCLUDE_PATH})
include_directories(${BOOST_INCLUDE_DIRS})
link_directories(${LINK_ROSETTA_LIST})
add_executable(test_app ${sources})
target_compile_definitions(test_app PRIVATE USING_MPI=1)
set_target_properties(test_app PROPERTIES COMPILE_DEFINITIONS "${_definitions}")
target_compile_options(test_app PRIVATE ${Open_MP_CXX_FLAGS})
target_link_libraries(test_app devel protocols.8 protocols.7 protocols_e.6 protocols_d.6 protocols_c.6 protocols_b.6 protocols_a.6 protocols.1 protocols_h.5 protocols_g.5 protocols_f.5 protocols_e.5 protocols_d.5 protocols_c.5 protocols_b.5 protocols_a.5 protocols.4  protocols_a.2 protocols_b.2 protocols.3 core.5 core.4 core.3 core.2 core.1 basic numeric utility ObjexxFCL z cppdb xml2 sqlite3 cifparse mpi boost_date_time boost_serialization boost_mpi ${Open_MP_CXX_FLAGS} ${MPI_LIBRARIES})

