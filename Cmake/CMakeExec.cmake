cmake_minimum_required(VERSION 3.12)  # Ensure the policy exists
cmake_policy(SET CMP0074 NEW)        # Use the new behavior
set(USE_KADATH_CMAKE_MODULES OFF)
if(USE_KADATH_CMAKE_MODULES)
  set(CMAKE_MODULE_PATH $ENV{HOME_KADATH}/Cmake)
endif()
if(EXISTS $ENV{HOME_KADATH}/Cmake/CMakeLocal.cmake)
	include ($ENV{HOME_KADATH}/Cmake/CMakeLocal.cmake)
endif()

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

option(PAR_VERSION "Parallel version" ON)
option(MKL_VERSION "Parallel version" OFF)

if (PAR_VERSION)
	message ("Parallel version")
else(PAR_VERSION)
	message ("Sequential version")
endif(PAR_VERSION)

if (MKL_VERSION)
	message ("Version with MKL ")
else(MKL_VERSION)
	message ("Version with scalapack")
endif(MKL_VERSION)

file(GLOB_RECURSE HEADERS ${CMAKE_SOURCE_DIR}/include/*.hpp)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE})

include_directories($ENV{HOME_KADATH}/include)
include_directories($ENV{HOME_KADATH}/src/Array)
include_directories($ENV{HOME_KADATH}/include/Kadath_point_h)

#Get the correct Kadath library (default == Release)
if (CMAKE_BUILD_TYPE MATCHES Debug)
	set(LIB_KADATH $ENV{HOME_KADATH}/lib/libkadath-debug.a)
else()
	set(LIB_KADATH $ENV{HOME_KADATH}/lib/libkadath.a)
endif()

#If parallel need to use the PMI wrapper
if(PAR_VERSION)
	find_package (MPI REQUIRED)
	set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
	message ("MPI CXX " ${MPI_CXX_COMPILER})
endif(PAR_VERSION)

if (NOT DEFINED GSL_LIBRARIES)
     find_package(GSL REQUIRED) #No FindGSL for cmake 2.8
endif()

if (NOT DEFINED FFTW_LIBRARIES)
     include(FindFFTW) #idem for FindFFTW
endif()

if (NOT PAR_VERSION)
if (NOT DEFINED LAPACK_LIBRARIES)
    include(FindLAPACK) #FindLAPACK does exist
endif()
endif()

#include(FindSUNDIALS)

#Need scalapack if PAR VERSION
if (PAR_VERSION)
	if (MKL_VERSION)
		include(FindMKL)
	else(MKL_VERSION)
	if (NOT DEFINED SCALAPACK_LIBRARIES)
		message ("Need to give scalapack location in CMakeLocal.cmake")
	endif()
	endif(MKL_VERSION)
endif(PAR_VERSION)

if(NOT PAR_VERSION)
#    if (NOT DEFINED PGPLOT_LIBRARIES)
#	    find_package(PGPLOT REQUIRED) #pgplot probably only used in sequential mode
#    endif()

    if (NOT DEFINED SUNDIALS_LIBRARIES)
	    include(FindSUNDIALS) #SUNDIAL probably only used in sequential mode
    endif()
endif(NOT PAR_VERSION)

option(GRHAYL_EOS "Use GRHaYL EOS library" OFF)
if (GRHAYL_EOS)
	set(CMAKE_MODULE_PATH $ENV{HOME_KADATH}/Cmake)
	include(FindGRHAYL)
	find_package(GRHAYL)
	add_definitions(-DWITH_GRHAYL_EOS)
	message("GRHAYL_INCLUDES: " ${GRHAYL_INCLUDES})
	include_directories (${GRHAYL_INCLUDES})
	unset(CMAKE_MODULE_PATH)
endif()

#need to use C++17
add_definitions(-std=c++17)
if ("$ENV{KAD_ARCH}" STREQUAL "")
  set(CMAKE_CXX_FLAGS_RELEASE "-g -O3 -std=c++17")
else()
  set(CMAKE_CXX_FLAGS_RELEASE "-g -O3 -std=c++17 -march=$ENV{KAD_ARCH}")
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DNDEBUG")
message("CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
message("CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")
