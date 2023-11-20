cmake_minimum_required(VERSION 3.25)
set(CMAKE_CXX_COMPILER "/usr/bin/g++")
message("${CMAKE_CXX_COMPILER}")
project(GRAPH_MINOR)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fopenmp -pthread -O3 -ftree-vectorize")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories(lib/inc)

add_library(library lib/src/mmio.cpp)

file(GLOB SOURCES "src/*.cpp")
list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/pthreadsDynSched.cpp")

add_executable(output main.cpp ${SOURCES})
target_include_directories(output PUBLIC inc)
target_link_libraries(output PRIVATE library)


