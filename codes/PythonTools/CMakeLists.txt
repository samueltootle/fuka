cmake_minimum_required(VERSION 2.8)
project(myproject)

include ($ENV{HOME_KADATH}/Cmake/CMakeExec.cmake)

# overwrite Kadath linking to get shared lib working
set(LIB_KADATH -Wl,--whole-archive $ENV{HOME_KADATH}/lib/libkadath.a -Wl,--no-whole-archive)
include_directories(/usr/include)

find_package(Python3 COMPONENTS Development)
if(Python3_FOUND)
  include_directories(${Python3_INCLUDE_DIRS})
  find_package(Boost COMPONENTS python REQUIRED)
  if(Boost_FOUND)
    #Library 1
    set (LIB1 pyConfigurator)
    set (SOURCESLIB src/${LIB1}.cpp)
    add_library(${LIB1} SHARED ${SOURCESLIB} ${HEADERS})
    target_link_libraries(${LIB1} ${LIB_KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES} ${BLACS_LIBRARIES} Boost::python)
    set_target_properties(${LIB1} PROPERTIES PREFIX "_")
    set_target_properties(${LIB1} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ../lib/fukaID_readers/${LIB1}/) 
    
    #Library 1
    set (LIB1 bbh_reader)
    set (LIB1MOD bbh)
    set (SOURCESLIB src/${LIB1}.cpp)
    add_library(${LIB1} SHARED ${SOURCESLIB} ${HEADERS})
    target_link_libraries(${LIB1} ${LIB_KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES} ${BLACS_LIBRARIES} Boost::python)
    set_target_properties(${LIB1} PROPERTIES PREFIX "_")
    set_target_properties(${LIB1} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ../lib/fukaID_readers/${LIB1MOD}/) 
  
    #Library 1
    set (LIB1 bh_reader)
    set (LIB1MOD bh)
    set (SOURCESLIB src/${LIB1}.cpp)
    add_library(${LIB1} SHARED ${SOURCESLIB} ${HEADERS})
    target_link_libraries(${LIB1} ${LIB_KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES} ${BLACS_LIBRARIES} Boost::python)
    set_target_properties(${LIB1} PROPERTIES PREFIX "_")
    set_target_properties(${LIB1} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ../lib/fukaID_readers/${LIB1MOD}/) 

    #Library 1
    set (LIB1 bhns_reader)
    set (LIB1MOD bhns)    
    set (SOURCESLIB src/${LIB1}.cpp)
    add_library(${LIB1} SHARED ${SOURCESLIB} ${HEADERS})
    target_link_libraries(${LIB1} ${LIB_KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES} ${BLACS_LIBRARIES} Boost::python)
    set_target_properties(${LIB1} PROPERTIES PREFIX "_")
    set_target_properties(${LIB1} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ../lib/fukaID_readers/${LIB1MOD}/) 

    #Library 1
    set (LIB1 bns_reader)
    set (LIB1MOD bns)    
    set (SOURCESLIB src/${LIB1}.cpp)
    add_library(${LIB1} SHARED ${SOURCESLIB} ${HEADERS})
    target_link_libraries(${LIB1} ${LIB_KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES} ${BLACS_LIBRARIES} Boost::python)
    set_target_properties(${LIB1} PROPERTIES PREFIX "_")
    set_target_properties(${LIB1} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ../lib/fukaID_readers/${LIB1MOD}/) 
   
    #Library 1
    set (LIB1 ns_reader)
    set (LIB1MOD ns)    
    set (SOURCESLIB src/${LIB1}.cpp)
    add_library(${LIB1} SHARED ${SOURCESLIB} ${HEADERS})
    target_link_libraries(${LIB1} ${LIB_KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES} ${BLACS_LIBRARIES} Boost::python)
    set_target_properties(${LIB1} PROPERTIES PREFIX "_")
    set_target_properties(${LIB1} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ../lib/fukaID_readers/${LIB1MOD}/) 
    
    #Library 1
    set (LIB1 ns_isotropic_diffrot_reader)
    set (LIB1MOD ns)    
    set (SOURCESLIB src/${LIB1}.cpp)
    add_library(${LIB1} SHARED ${SOURCESLIB} ${HEADERS})
    target_link_libraries(${LIB1} ${LIB_KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES} ${BLACS_LIBRARIES} Boost::python)
    set_target_properties(${LIB1} PROPERTIES PREFIX "_")
    set_target_properties(${LIB1} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ../lib/fukaID_readers/${LIB1MOD}/) 

    #Library 1
    set (LIB1 ns_isotropic_uniformrot_reader)
    set (LIB1MOD ns)    
    set (SOURCESLIB src/${LIB1}.cpp)
    add_library(${LIB1} SHARED ${SOURCESLIB} ${HEADERS})
    target_link_libraries(${LIB1} ${LIB_KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES} ${BLACS_LIBRARIES} Boost::python)
    set_target_properties(${LIB1} PROPERTIES PREFIX "_")
    set_target_properties(${LIB1} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ../lib/fukaID_readers/${LIB1MOD}/) 
    
    #Library 1
    set (LIB1 ns_isotropic_norot_reader)
    set (LIB1MOD ns)    
    set (SOURCESLIB src/${LIB1}.cpp)
    add_library(${LIB1} SHARED ${SOURCESLIB} ${HEADERS})
    target_link_libraries(${LIB1} ${LIB_KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES} ${BLACS_LIBRARIES} Boost::python)
    set_target_properties(${LIB1} PROPERTIES PREFIX "_")
    set_target_properties(${LIB1} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ../lib/fukaID_readers/${LIB1MOD}/) 

    #Library 1
    set (LIB1 ns_diffrot_reader)
    set (LIB1MOD ns)    
    set (SOURCESLIB src/${LIB1}.cpp)
    add_library(${LIB1} SHARED ${SOURCESLIB} ${HEADERS})
    target_link_libraries(${LIB1} ${LIB_KADATH} ${PGPLOT_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES} ${BLACS_LIBRARIES} Boost::python)
    set_target_properties(${LIB1} PROPERTIES PREFIX "_")
    set_target_properties(${LIB1} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ../lib/fukaID_readers/${LIB1MOD}/) 


  else()
    message("Boost/Boost Python libraries not found. pyconfig cannot be built.")
  endif()
else()
  message("Python3 developmnent components not found. pyconfig cannot be built.")
endif()


