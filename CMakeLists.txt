cmake_minimum_required(VERSION 3.9)

project (EvolutionRosetta)

# Rosetta flags for compatibility
############ Build Type, Mode, Link Rules, Compiler, Extras
SET(EXTRAS cxx11)
SET(MODE debug)
SET(LINK_RULES shared)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin/)
SET(CMAKE_INSTALL_PREFIX ../bin/)

REMOVE_DEFINITIONS(-DPTR_MODERN)
REMOVE_DEFINITIONS(-DPTR_BOOST)

ADD_DEFINITIONS(-DCXX11)
ADD_DEFINITIONS(-DPTR_STD)

SET(COMPILE_FLAGS "${COMPILE_FLAGS} -Wno-error")

# Rosetta include files
include_directories(/home/dvarela/Code/rosetta_2018/main/source/cmake/build_cxx11)
set(INCLUDE_DIRS /home/dvarela/Code/rosetta_2018/main/source/src/)

# Rosetta link library
find_library(rosetta
    NAMES basic cifparse devel core.1 core.2 core.3 core.4 core.5 utility cppdb libxml2 numeric ObjexxFCL protocols.1 protocols.3 protocols.4 protocols.7 protocols.8 protocols_a.2 protocols_a.5 protocols_a.6 protocols_b.2 protocols_b.5 protocols_b.6 protocols_c.5 protocols_c.6 protocols_d.5 protocols_d.6 protocols_e.5 protocols_e.6 protocols_f.5 protocols_g.5 protocols_h.5 sqlite3
    HINTS "/home/dvarela/Code/rosetta_2018/main/source/cmake/build_cxx11/"
)

# My app source files
file(GLOB_RECURSE sources src/*.cpp src/*.h)

add_executable(app ${sources} ${rosetta})
target_include_directories(app PRIVATE ${INCLUDE_DIRS})

