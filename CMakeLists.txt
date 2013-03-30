# Simple CMake file

project(sac-dev)

cmake_minimum_required(VERSION 2.8)

option(WITH_DMALLOC "Use dmalloc for memory profiling" OFF)

option(WITH_MALLOC_COUNT "Use malloc_count for memory profiling" ON)

# disallow in-source builds

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

# compilation options 

option(BUILD_ALL "Build different test programs" OFF)

option(BUILD_ALL_TUPLESIZE "Build different test programs" OFF)

option(BUILD_ALL_MEMLIMIT "Build different test programs" OFF)

option(BUILD_ALL_BLOCKSIZE "Build different test programs" OFF)

# Use OpenMP

include(FindOpenMP)

if(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
  message(FATAL_ERROR "OpenMP not found.")
endif()

# Enable warnings

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -ftemplate-depth-1000")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")

# Use zlib

include(FindZLIB)

include_directories(${ZLIB_INCLUDE_DIR})

# build stxxl subproject

add_subdirectory(stxxl)

# descend into source

add_subdirectory(src)

