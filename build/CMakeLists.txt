cmake_minimum_required(VERSION 2.6)
project(fhmc)

set(CMAKE_CXX_STANDARD 11)

# Check that C++11 is supported
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    message(STATUS "${COMPILER_SUPPORTS_CXX11}")
else(COMPILER_SUPPORTS_CXX11)
    message(FATAL_ERROR "${COMPILER_SUPPORTS_CXX11}")
endif(COMPILER_SUPPORTS_CXX11)

# C++11 support for each compiler
if(CMAKE_COMPILER_IS_GNUCXX)
	# gcc
	add_definitions(-std=c++11)
else()
	# clang
	add_definitions(-std=c++11 -stdlib=libc++)
endif()

file(GLOB fhmc_SRC
	"../src/*.h"
	"../src/*.cpp"
	"../driver/main-tmmc.cpp"
)

set (fhmc_VERSION_MAJOR 2)
set (fhmc_VERSION_MINOR 1)

set(CMAKE_CXX_FLAGS "-O3 -Wno-sign-compare -Wall -DFLUID_PHASE_SIMULATIONS")

include_directories("../src/")
include_directories("../driver/")

# If you desire to compile statically, use the following 3 lines
SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
SET(BUILD_SHARED_LIBRARIES OFF)
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")

add_executable(fhmc_tmmc ${fhmc_SRC})
set(EXECUTABLE_OUTPUT_PATH "../bin")
