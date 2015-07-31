cmake_minimum_required(VERSION 2.6)
project(unitTests)

################################
# GTest
################################
ADD_SUBDIRECTORY (gtest-1.7.0)
enable_testing()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# do not include main.cpp
include_directories("../../../src/")
file(GLOB omcs_SRC
	"../../../src/*.h"
	"../../../src/delete.cpp"
	"../../../src/global.h"
	"../../../src/atom.h"
	"../../../src/input.h"
	"../../../src/input.cpp"
	"../../../src/insert.h"
	"../../../src/insert.cpp"
	"../../../src/moves.h"
	"../../../src/moves.cpp"
	"../../../src/potentials.h"
	"../../../src/potentials.cpp"
	"../../../src/system.h"
	"../../../src/system.cpp"
	"../../../src/utilities.h"
	"../../../src/utilities.cpp"
)

################################
# Unit Tests
################################
add_executable(runUnitTests test.cc ${omcs_SRC})
# Link test executable against gtest & gtest_main
target_link_libraries(runUnitTests gtest gtest_main)
add_test( runUnitTests runUnitTests)

set(EXECUTABLE_OUTPUT_PATH "../bin/")
