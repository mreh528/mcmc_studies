cmake_minimum_required(VERSION 3.10)

# Set project name and version
project(mcmc_studies VERSION 1.0)

# Specify C++ standard
option(CXX_WARNINGS "Enable most C++ warning flags." ON)
option(COLOR_OUTPUT "Enable colored terminal output." ON)

# Add the ROOT path to the cmake path
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

# Setup some IO paths
set(LIBRARY_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/bin")

# Compiler configs
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -std=c++11 -m64")

# Files that need to be included in all executables
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/utils)

# Files that need to be included only in specific ("target") executables
#target_include_directories()

# Find ROOT and make sure we have some necessary components
find_package(ROOT 5 REQUIRED COMPONENTS RIO Net)
if (ROOT_FOUND)
    message(STATUS "ROOT found at $ENV{ROOTSYS}")
    include_directories(${ROOT_INCLUDE_DIR})
    include(${ROOT_USE_FILE})
    link_directories(${ROOT_LIBRARY_DIR})
else()
    message(FATAL_ERROR "Correct ROOT version not found!")
endif()

# Check on things
message(STATUS "C flags      : " ${CMAKE_C_FLAGS})
message(STATUS "CXX flags    : " ${CMAKE_CXX_FLAGS})
message(STATUS "Linker flags : " ${CMAKE_EXE_LINKER_FLAGS}\n)

# Add active subdirectories
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/utils)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/macros)

