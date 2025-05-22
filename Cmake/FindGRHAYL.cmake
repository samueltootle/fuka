# - Find the GRHAYL library
#
# Usage:
#   find_package(GRHAYL [REQUIRED] [QUIET] )
#
# It sets the following variables:
#   GRHAYL_FOUND               ... true if grhayl is found on the system
#   GRHAYL_LIBRARIES           ... full path to grhayl library
#   GRHAYL_INCLUDES            ... grhayl include directory
#
# The following variables will be checked by the function
#   GRHAYL_USE_STATIC_LIBS    ... if true, only static libraries are found
#   GRHAYL_ROOT               ... if set, the libraries are exclusively searched
#                               under this path
#   GRHAYL_LIBRARY            ... grhayl library to use
#   GRHAYL_INCLUDE_DIR        ... grhayl include directory
#

#If environment variable GRHAYLDIR is specified, it has same effect as GRHAYL_ROOT
if( NOT GRHAYL_ROOT AND DEFINED ENV{GRHAYLDIR} )
  set( GRHAYL_ROOT $ENV{GRHAYLDIR} )
elseif( NOT GRHAYL_ROOT AND DEFINED ENV{GRHAYL_ROOT} )
  set( GRHAYL_ROOT $ENV{GRHAYL_ROOT} )
endif()

# Check if we can use PkgConfig
find_package(PkgConfig)

# #Determine from PKG
if( PKG_CONFIG_FOUND AND NOT GRHAYL_ROOT )
  pkg_check_modules( PKG_GRHAYL QUIET "grhayl" )
endif()

#Check whether to search static or dynamic libs
set( CMAKE_FIND_LIBRARY_SUFFIXES_SAV ${CMAKE_FIND_LIBRARY_SUFFIXES} )

if( ${GRHAYL_USE_STATIC_LIBS} )
  set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX} )
else()
  set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_SHARED_LIBRARY_SUFFIX} )
endif()

if( GRHAYL_ROOT )
  #find libs
  find_library(
    GRHAYL_LIB
    NAMES "grhayl"
    PATHS ${GRHAYL_ROOT}
    PATH_SUFFIXES "lib" "lib64"
    NO_DEFAULT_PATH
  )
else()

  find_library(
    GRHAYL_LIB
    NAMES "grhayl"
    PATHS ${PKG_GRHAYL_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
  )
endif( GRHAYL_ROOT )
set(GRHAYL_INCLUDES ${GRHAYL_ROOT}/include)

set(GRHAYL_LIBRARIES ${GRHAYL_LIB})

set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_SAV} )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GRHAYL DEFAULT_MSG
                                  GRHAYL_INCLUDES GRHAYL_LIBRARIES)

mark_as_advanced(GRHAYL_INCLUDES GRHAYL_LIBRARIES GRHAYL_LIB)
