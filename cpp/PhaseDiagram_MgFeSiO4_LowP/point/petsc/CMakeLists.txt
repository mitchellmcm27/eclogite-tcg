cmake_minimum_required(VERSION 2.8)

set(PROJECT_NAME phasediagram)
project(${PROJECT_NAME})

# Set verbose output while testing CMake
#set(CMAKE_VERBOSE_MAKEFILE 1)

# Set CMake behavior
cmake_policy(SET CMP0004 NEW)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
find_package(PETSc REQUIRED)

include_directories(SYSTEM ${PETSC_INCLUDES})

# Executable
add_executable(${PROJECT_NAME} main.cpp)

# Target libraries
target_link_libraries(${PROJECT_NAME} ${PETSC_LIBRARIES})
target_link_libraries(${PROJECT_NAME} MgFeSiO4_all_slb_rx)
