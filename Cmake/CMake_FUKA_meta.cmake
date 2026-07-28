
project(FUKA VERSION 2.4.0)
find_package(Git QUIET)

set(GIT_HASH "unknown")
set(GIT_DESCRIBE "unknown")
set(GIT_DIRTY 0)


if(GIT_FOUND)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
        WORKING_DIRECTORY $ENV{HOME_KADATH}
        OUTPUT_VARIABLE GIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --always --tags --dirty
        WORKING_DIRECTORY $ENV{HOME_KADATH}
        OUTPUT_VARIABLE GIT_DESCRIBE
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    execute_process(
        COMMAND ${GIT_EXECUTABLE} diff-index --quiet HEAD --
        WORKING_DIRECTORY $ENV{HOME_KADATH}
        RESULT_VARIABLE GIT_DIFF_RESULT
    )

    if(GIT_DIFF_RESULT EQUAL 0)
        set(GIT_DIRTY false)
    else()
        set(GIT_DIRTY true)
    endif()

    string(TIMESTAMP BUILD_DATE "%Y-%m-%dT%H:%M:%SZ" UTC)
endif()

configure_file(
    $ENV{HOME_KADATH}/Cmake/fuka_version.hpp.in
    $ENV{HOME_KADATH}/include/FUKA_Solvers/utilities/fuka_version.hpp
    @ONLY
)