cmake_minimum_required(VERSION 3.3)
project(PumpKin)

INCLUDE_DIRECTORIES(/opt/local/include)
LINK_DIRECTORIES(/opt/local/lib)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES pk_Array.cpp pk_IO.cpp pk_Pathways.cpp PumpKin.cpp pk_Print.cpp)
add_executable(PumpKin ${SOURCE_FILES})
target_link_libraries(PumpKin glpk)
