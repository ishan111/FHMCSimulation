cmake_minimum_required(VERSION 2.6)
project(unitTests)

################################
# GTest
################################
ADD_SUBDIRECTORY (gtest-1.7.0)
enable_testing()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

set(CMAKE_CXX_STANDARD 11)

# check that C++11 is supported
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

add_definitions(-D_GLIBCXX_DEBUG)
#add_definitions(-fsanitize=address -mllvm -asan-stack)

set(CMAKE_CXX_FLAGS "-Wall -Wno-sign-compare")# -std=c++11 -stdlib=libc++")

# do not include main.cpp
include_directories("../../src/")
file(GLOB omcs_SRC
	"../../src/*.h"
	"../../srd/atom.h"
	"../../src/barrier.cpp"
	"../../src/bias.cpp"
	"../../src/cellList.cpp"
	"../../src/crossover.cpp"
	"../../src/delete.cpp"
	"../../src/histogram.cpp"
	"../../src/input.cpp"
	"../../src/insert.cpp"
	"../../src/mover.cpp"
	"../../src/moves.cpp"
	"../../src/potentials.cpp"
	"../../src/quaternion.cpp"
	"../../src/checkpoint.cpp"
	"../../src/sanity.cpp"
	"../../src/swap.cpp"
	"../../src/system.cpp"
	"../../src/tmmc.cpp"
	"../../src/translate.cpp"
	"../../src/utilities.cpp"
	"../../src/wala.cpp"	
)

################################
# Unit Tests
################################
add_executable(runUnitTests test.cc ${omcs_SRC})
# Link test executable against gtest & gtest_main
target_link_libraries(runUnitTests gtest gtest_main)
add_test( runUnitTests runUnitTests)

set(EXECUTABLE_OUTPUT_PATH "../bin/")
