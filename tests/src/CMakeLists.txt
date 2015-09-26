cmake_minimum_required(VERSION 2.6)
project(unitTests)

################################
# GTest
################################
ADD_SUBDIRECTORY (gtest-1.7.0)
enable_testing()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})


if(CMAKE_COMPILER_IS_GNUCXX)
    add_definitions(-std=gnu++0x)
endif()

# do not include main.cpp
include_directories("../../src/")
file(GLOB omcs_SRC
	"../../src/*.h"
	"../../src/barrier.cpp"
	"../../src/bias.cpp"
	"../../src/cellList.cpp"
	"../../src/delete.cpp"
	"../../src/histogram.cpp"
	"../../src/insert.cpp"
	"../../src/moves.cpp"
	"../../src/potentials.cpp"
	"../../src/swap.cpp"
	"../../src/system.cpp"
	"../../src/translate.cpp"
	"../../src/utilities.cpp"
)

################################
# Unit Tests
################################
add_executable(runUnitTests test.cc ${omcs_SRC})
# Link test executable against gtest & gtest_main
target_link_libraries(runUnitTests gtest gtest_main)
add_test( runUnitTests runUnitTests)

set(EXECUTABLE_OUTPUT_PATH "../bin/")
