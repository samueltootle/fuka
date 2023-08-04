cmake_minimum_required(VERSION 2.8)
project(myproject)

include ($ENV{HOME_KADATH}/Cmake/CMakeExec.cmake)

#Executable 1
set (PROG1 solve)
set (SOURCES1 src/${PROG1}.cpp)
add_executable(${PROG1} ${SOURCES1} ${HEADERS})
target_link_libraries(${PROG1} ${LIB_KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES} ${BLACS_LIBRARIES})

set (PROG2 reader)
set (SOURCES5 src/${PROG2}.cpp)
add_executable(${PROG2} ${SOURCES5} ${HEADERS})
target_link_libraries(${PROG2} ${LIB_KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES} ${BLACS_LIBRARIES})