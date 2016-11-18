# This does not compile with the netCDF libraries so all output will be switched to ASCII from the executable
cmake_minimum_required(VERSION 2.6)
project(fhmc)

if(CMAKE_COMPILER_IS_GNUCXX)
    add_definitions(-std=gnu++0x)
endif()

file(GLOB fhmc_SRC
	"../src/*.h"
	"../src/*.cpp"
	"../driver/main-tmmc.cpp"
)

set (fhmc_VERSION_MAJOR 2)
set (fhmc_VERSION_MINOR 1)

set(CMAKE_CXX_FLAGS "-O3 -Wall -DFLUID_PHASE_SIMULATIONS")

include_directories("../src/")
include_directories("../driver/")
add_executable(fhmc_tmmc ${fhmc_SRC})
set(EXECUTABLE_OUTPUT_PATH "../bin")
