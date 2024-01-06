cmake_minimum_required(VERSION 2.8)
project(Kadath_release)
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

if (EXISTS $ENV{HOME_KADATH}/Cmake/CMakeLocal.cmake)
	include($ENV{HOME_KADATH}/Cmake/CMakeLocal.cmake)
endif()

#Parallel or sequential
option(PAR_VERSION "Parallel version" ON)

if (PAR_VERSION)
	message ("Parallel version")
else(PAR_VERSION)
	message ("Sequential version")
endif(PAR_VERSION)

set(CMAKE_MODULE_PATH $ENV{HOME_KADATH}/build_release)

if(PAR_VERSION)
	file(COPY $ENV{HOME_KADATH}/src_par/do_newton.cpp DESTINATION $ENV{HOME_KADATH}/src/System_of_eqs/)
	find_package(MPI REQUIRED)
	set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
	message ("MPI Compiler " ${MPI_CXX_COMPILER})
else(PAR_VERSION)
	file(COPY $ENV{HOME_KADATH}/src_seq/do_newton.cpp DESTINATION $ENV{HOME_KADATH}/src/System_of_eqs/)
endif(PAR_VERSION)

#Get all the src files
file(GLOB_RECURSE KADATH_SOURCES $ENV{HOME_KADATH}/src/*.cpp )
#message ("Source files" : ${KADATH_SOURCES})

include_directories($ENV{HOME_KADATH}/include)
include_directories($ENV{HOME_KADATH}/include/Kadath_point_h)
#for templates in Array
include_directories($ENV{HOME_KADATH}/src/Array)

if (DEFINED FFTW_INCLUDE_DIR)
	include_directories (${FFTW_INCLUDE_DIR})
endif()

#need to use C++17
add_definitions(-std=c++17)
if ("$ENV{KAD_ARCH}" STREQUAL "")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -g -O3 -std=c++17")
else()
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -g -O3 -march=$ENV{KAD_ARCH} -std=c++17")
endif()

add_library(kadath ${KADATH_SOURCES})
set_target_properties(kadath PROPERTIES ARCHIVE_OUTPUT_DIRECTORY $ENV{HOME_KADATH}/lib COMPILE_FLAGS "-fPIC")

message("CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")
