cmake_minimum_required(VERSION 2.6)
project(omcs)

file(GLOB omcs_SRC
	"../src/*.h"
	"../src/*.cpp"
)

set (omcs_VERSION_MAJOR 2)
set (omcs_VERSION_MINOR 1)

set(CMAKE_CXX_FLAGS "-O3 -Wall -DNETCDF_CAPABLE -DFLUID_PHASE_SIMULATIONS")

include_directories("../src/")
add_executable(omcs ${omcs_SRC})
set(EXECUTABLE_OUTPUT_PATH "../bin")
