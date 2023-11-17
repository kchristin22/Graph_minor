cmake_minimum_required(VERSION 3.25)
# set(CMAKE_CXX_COMPILER "usr/bin/clang++")
set(CMAKE_CXX_COMPILER "/opt/opencilk/bin/clang")
message("${CMAKE_CXX_COMPILER}")
project(CILK)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fopencilk -O3 -ftree-vectorize")

# include_directories(lib/inc)

# add_library(library lib/src/mmio.cpp)

file(GLOB SOURCES "src/*.cpp")

add_executable(cilk main.cpp ${SOURCES})
target_include_directories(cilk PUBLIC inc)