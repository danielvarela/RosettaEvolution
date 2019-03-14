cmake_minimum_required(VERSION 3.9)

project (EvolutionRosetta)

set (CMAKE_CXX_STANDARD 11)

# Rosetta flags for compatibility
############ Build Type, Mode, Link Rules, Compiler, Extras
SET(EXTRAS cxx11)
SET(MODE debug)
SET(LINK_RULES shared)
SET(CMAKE_BUILD_TYPE Debug)

find_package(OpenMP)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#set(BOOST_ROOT "/home/dvarela/Code/rosetta/main/source/external/boost_1_55_0/")

find_package(Boost 1.61.0 COMPONENTS serialization mpi REQUIRED )
if(Boost_FOUND)
  message("Boost found")
endif()


find_package(MPI)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin/)
SET(CMAKE_INSTALL_PREFIX ../bin/)

REMOVE_DEFINITIONS(-DPTR_MODERN)
REMOVE_DEFINITIONS(-DPTR_BOOST)

ADD_DEFINITIONS(-DCXX11)
ADD_DEFINITIONS(-DPTR_STD)

SET(COMPILE_FLAGS "${COMPILE_FLAGS} -g -Wno-error -Wno-deprecated-declarations -lz")

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
file(GLOB_RECURSE sources src/main.cpp src/*.h src/helper_apps/*.cpp src/helper_apps/*.hpp src/DifferentialEvolution/*.cpp src/DifferentialEvolution/*.hpp src/PoseFiles/*.cpp src/PoseFiles/*.hpp src/PoseFiles/moves/*.cpp src/PoseFiles/moves/*.hpp) 



include_directories(${INCLUDE_ROSETTA_LIST})
#include_directories(${MPI_INCLUDE_PATH})
link_directories(${LINK_ROSETTA_LIST})
include_directories(${BOOST_INCLUDE_DIRS})
add_executable(app ${sources})
target_compile_options(app PRIVATE -Wall ${Open_MP_CXX_FLAGS})
#target_link_libraries(app devel protocols.8 protocols.7 protocols_e.6 protocols_d.6 protocols_c.6 protocols_b.6 protocols_a.6 protocols_h.5 protocols_g.5 protocols_f.5 protocols_e.5 protocols_d.5 protocols_c.5 protocols_b.5 protocols_a.5 protocols.4 protocols.3 protocols_b.2 protocols_a.2 protocols.1 core.5 core.4 core.3 core.2 core.1 basic numeric utility ObjexxFCL z cifparse)
target_link_libraries(app devel protocols.8 protocols.7 protocols_e.6 protocols_d.6 protocols_c.6 protocols_b.6 protocols_a.6 protocols.1 protocols_h.5 protocols_g.5 protocols_f.5 protocols_e.5 protocols_d.5 protocols_c.5 protocols_b.5 protocols_a.5 protocols.4  protocols_a.2 protocols_b.2 protocols.3 core.5 core.4 core.3 core.2 core.1 basic numeric utility ObjexxFCL z cppdb xml2 sqlite3 cifparse boost_date_time ${Open_MP_CXX_FLAGS})


set(_definitions USE_CRYO_EM=true)
include_directories(${INCLUDE_ROSETTA_LIST})
#include_directories(${MPI_INCLUDE_PATH})
link_directories(${LINK_ROSETTA_LIST})
include_directories(${BOOST_INCLUDE_DIRS})
add_executable(app_dens ${sources})
set_target_properties(app_dens PROPERTIES COMPILE_DEFINITIONS "${_definitions}")
target_compile_options(app_dens PRIVATE -Wall ${Open_MP_CXX_FLAGS})
#target_link_libraries(app_dens devel protocols.8 protocols.7 protocols_e.6 protocols_d.6 protocols_c.6 protocols_b.6 protocols_a.6 protocols_h.5 protocols_g.5 protocols_f.5 protocols_e.5 protocols_d.5 protocols_c.5 protocols_b.5 protocols_a.5 protocols.4 protocols.3 protocols_b.2 protocols_a.2 protocols.1 core.5 core.4 core.3 core.2 core.1 basic numeric utility ObjexxFCL z cifparse)
target_link_libraries(app_dens devel protocols.8 protocols.7 protocols_e.6 protocols_d.6 protocols_c.6 protocols_b.6 protocols_a.6 protocols.1 protocols_h.5 protocols_g.5 protocols_f.5 protocols_e.5 protocols_d.5 protocols_c.5 protocols_b.5 protocols_a.5 protocols.4  protocols_a.2 protocols_b.2 protocols.3 core.5 core.4 core.3 core.2 core.1 basic numeric utility ObjexxFCL z cppdb xml2 sqlite3 cifparse boost_date_time ${Open_MP_CXX_FLAGS})



set(_definitions MPI_ENABLED=true)
file(GLOB_RECURSE mpi_sources src/mpi_main.cpp src/mpi_files/*.hpp src/mpi_files/*.cpp src/*.h src/DifferentialEvolution/*.cpp src/DifferentialEvolution/*.hpp src/PoseFiles/*.cpp src/PoseFiles/*.hpp src/PoseFiles/moves/*.cpp src/PoseFiles/moves/*.hpp) 
include_directories(${INCLUDE_ROSETTA_LIST})
include_directories(${MPI_INCLUDE_PATH})
include_directories(${BOOST_INCLUDE_DIRS})
link_directories(${LINK_ROSETTA_LIST})
add_executable(mpi_exec ${mpi_sources})
target_compile_definitions(mpi_exec PRIVATE USING_MPI=1)
set_target_properties(mpi_exec PROPERTIES COMPILE_DEFINITIONS "${_definitions}")
target_compile_options(mpi_exec PRIVATE ${Open_MP_CXX_FLAGS})
target_link_libraries(mpi_exec devel protocols.8 protocols.7 protocols_e.6 protocols_d.6 protocols_c.6 protocols_b.6 protocols_a.6 protocols.1 protocols_h.5 protocols_g.5 protocols_f.5 protocols_e.5 protocols_d.5 protocols_c.5 protocols_b.5 protocols_a.5 protocols.4  protocols_a.2 protocols_b.2 protocols.3 core.5 core.4 core.3 core.2 core.1 basic numeric utility ObjexxFCL z cppdb xml2 sqlite3 cifparse mpi boost_date_time boost_serialization boost_mpi ${Open_MP_CXX_FLAGS} ${MPI_LIBRARIES})


set(_definitions MPI_ENABLED=true USE_CRYO_EM=true)
include_directories(${INCLUDE_ROSETTA_LIST})
include_directories(${MPI_INCLUDE_PATH})
include_directories(${BOOST_INCLUDE_DIRS})
link_directories(${LINK_ROSETTA_LIST})
add_executable(mpi_dens ${mpi_sources})
target_compile_definitions(mpi_dens PRIVATE USING_MPI=1)
set_target_properties(mpi_dens PROPERTIES COMPILE_DEFINITIONS "${_definitions}")
target_compile_options(mpi_dens PRIVATE ${Open_MP_CXX_FLAGS})
target_link_libraries(mpi_dens devel protocols.8 protocols.7 protocols_e.6 protocols_d.6 protocols_c.6 protocols_b.6 protocols_a.6 protocols.1 protocols_h.5 protocols_g.5 protocols_f.5 protocols_e.5 protocols_d.5 protocols_c.5 protocols_b.5 protocols_a.5 protocols.4  protocols_a.2 protocols_b.2 protocols.3 core.5 core.4 core.3 core.2 core.1 basic numeric utility ObjexxFCL z cppdb xml2 sqlite3 cifparse mpi boost_date_time boost_serialization boost_mpi ${Open_MP_CXX_FLAGS} ${MPI_LIBRARIES})


file(GLOB_RECURSE app_builder_sources src/popul_builder_main.cpp src/helper_apps/*.hpp src/helper_apps/*.cpp src/*.h src/PoseFiles/*.cpp src/PoseFiles/*.hpp src/PoseFiles/moves/*.cpp src/PoseFiles/moves/*.hpp) 

include_directories(${INCLUDE_ROSETTA_LIST})
link_directories(${LINK_ROSETTA_LIST})
include_directories(${BOOST_INCLUDE_DIRS})
add_executable(popul_builder_app ${app_builder_sources})
target_compile_options(popul_builder_app PRIVATE -Wall ${Open_MP_CXX_FLAGS})
#target_link_libraries(popul_builder_app devel protocols.8 protocols.7 protocols_e.6 protocols_d.6 protocols_c.6 protocols_b.6 protocols_a.6 protocols_h.5 protocols_g.5 protocols_f.5 protocols_e.5 protocols_d.5 protocols_c.5 protocols_b.5 protocols_a.5 protocols.4 protocols.3 protocols_b.2 protocols_a.2 protocols.1 core.5 core.4 core.3 core.2 core.1 basic numeric utility ObjexxFCL z cifparse)
target_link_libraries(popul_builder_app devel protocols.8 protocols.7 protocols_e.6 protocols_d.6 protocols_c.6 protocols_b.6 protocols_a.6 protocols.1 protocols_h.5 protocols_g.5 protocols_f.5 protocols_e.5 protocols_d.5 protocols_c.5 protocols_b.5 protocols_a.5 protocols.4  protocols_a.2 protocols_b.2 protocols.3 core.5 core.4 core.3 core.2 core.1 basic numeric utility ObjexxFCL z cppdb xml2 sqlite3 cifparse boost_date_time ${Open_MP_CXX_FLAGS})


