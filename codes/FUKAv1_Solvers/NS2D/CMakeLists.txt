cmake_minimum_required(VERSION 2.8)
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.0)
    cmake_policy(SET CMP0048 NEW)
endif()
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.1)
    cmake_policy(SET CMP0054 NEW)
endif()
cmake_minimum_required(VERSION 2.8)
project(myproject)

include ($ENV{HOME_KADATH}/Cmake/CMakeExec.cmake)

#Executable 1
set (PROG1 setup)
set (PROG1_SOURCES src/${PROG1}.cpp)
add_executable(${PROG1} ${PROG1_SOURCES} ${HEADERS})
target_link_libraries(${PROG1} ${LIB_KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES})

#Executable 2
set (PROG1 solve)
set (PROG1_SOURCES src/${PROG1}.cpp)
add_executable(${PROG1} ${PROG1_SOURCES} ${HEADERS})
target_link_libraries(${PROG1} ${LIB_KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES})
#
##Executable 1
set (PROG5 reader)
set (SOURCES5 src/${PROG5}.cpp)
add_executable(${PROG5} ${SOURCES5} ${HEADERS})
target_link_libraries(${PROG5} ${LIB_KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES})
##Executable 1
#set (PROG5 increase_resolution)
#set (SOURCES5 src/${PROG5}.cpp)
#add_executable(${PROG5} ${SOURCES5} ${KADATH_HEADERS})
#target_link_libraries(${PROG5} ${KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES} ${BLACS_LIBRARIES})
