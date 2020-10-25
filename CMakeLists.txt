# CMake entry point                                                             
cmake_minimum_required (VERSION 3.5)
project (simple_paradis)

## Add headers files                                                            
include_directories ()



SET(CMAKE_CXX_FLAGS -pthread)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

find_package(OpenMP REQUIRED)

if(OpenMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
add_compile_options(
  -std=c++17
  -O3
  -march=native
  -fopenmp
  -latomic
  -mcx16
)


 
# main.cpp
add_executable(paradis_ompf_repairgen ${HEADERS} src/paradis_ompf_repairgen.cc)
add_executable(paradis src/paradis_test.cc)
