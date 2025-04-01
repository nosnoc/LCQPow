##
##	This file is part of LCQPow.
##
##	LCQPow -- A Solver for Quadratic Programs with Commplementarity Constraints.
##	Copyright (C) 2020 - 2022 by Jonas Hall et al.
##
##	LCQPow is free software; you can redistribute it and/or
##	modify it under the terms of the GNU Lesser General Public
##	License as published by the Free Software Foundation; either
##	version 2.1 of the License, or (at your option) any later version.
##
##	LCQPow is distributed in the hope that it will be useful,
##	but WITHOUT ANY WARRANTY; without even the implied warranty of
##	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##	See the GNU Lesser General Public License for more details.
##
##	You should have received a copy of the GNU Lesser General Public
##	License along with LCQPow; if not, write to the Free Software
##	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
##

# cmake minimum required version
cmake_minimum_required(VERSION 3.13)

# Project name
project(lcqpow CXX)

set(CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# set some default flags
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# set install prefix
set(
    CMAKE_INSTALL_PREFIX
    "${CMAKE_BINARY_DIR}"
)

## Options ------------------------------------------------------------------------------
option(
    BUILD_EXAMPLES
    "Option to build examples"
    ON
)

option(
    UNIT_TESTS
    "Option to build and run unit tests"
    ON
)

option(
    PROFILING
    "Option to enable profiling (gprof)."
    OFF
)

option(
    QPOASES_SCHUR
    "Use the Schur Complement Method with MA57 Solver"
    OFF
)

option(
    BUILD_MATLAB_INTERFACE
    "Option to build Matlab interface"
    ON
)

option(
    BUILD_PYTHON_INTERFACE
    "Option to build Python interface"
    ON
)

# print options
message(
    "CMAKE_BUILD_TYPE           ${CMAKE_BUILD_TYPE}\n"
    "EXAMPLES                   ${BUILD_EXAMPLES}\n"
    "BUILD_MATLAB_INTERFACE     ${BUILD_MATLAB_INTERFACE}\n"
    "BUILD_PYTHON_INTERFACE     ${BUILD_PYTHON_INTERFACE}\n"
    "UNIT_TESTS                 ${UNIT_TESTS}\n"
    "PROFILING                  ${PROFILING}\n"
    "QPOASES_SCHUR              ${QPOASES_SCHUR}\n"
)

## ADD ALL EXTERNAL PROJECTS ------------------------------------------------------------
# external modules
include(ExternalProject)

# 1) qpOASES
set(
    DEF_SOLVER
    "SOLVER_NONE"
)

if (${QPOASES_SCHUR})

    find_package(Matlab)
    if (NOT ${Matlab_FOUND})
        # throw error if schur is enabled and matlab not found
        message(
            FATAL_ERROR
            "qpOASES with Schur Complement method is currently only supported in combination with matlab, which can not be located."
        )
    endif()
    
    message(
        "qpOASES Matlab found"
    )

    # if found, then set some variables
    get_filename_component(
        MATLAB_LIBDIR
        ${Matlab_MEX_LIBRARY}
        DIRECTORY
    )

    set(
        LIB_SOLVER
        "${MATLAB_LIBDIR}/libmwma57.so"
    )

    set(
        DEF_SOLVER
        "SOLVER_MA57"
    )

    set(
        USE_SOLVER
        "MA57"
    )

    set(
        LINKHSL
        "-Wl,-rpath=${MATLAB_LIBDIR}"
    )

    set(
        LIB_LAPACK
        "${MATLAB_LIBDIR}/libmwlapack.so"
    )

    set(
        LIB_BLAS
        "${MATLAB_LIBDIR}/libmwblas.so"
    )

    set(
        LIB_METIS
        "${MATLAB_LIBDIR}/libmwmetis.so"
    )

    set(
        sparse_libs
        "${LIB_SOLVER} ${LIB_BLAS} ${LIB_LAPACK} ${LIB_METIS}"
    )

    set(
        qpOASES_MAKE_ARGS
        ${qpOASES_MAKE_ARGS} LIB_SOLVER=${LIB_SOLVER} LIB_BLAS=${LIB_BLAS} LIB_LAPACK=${LIB_LAPACK} LINKHSL=${LINKHSL}
    )

    if(NOT EXISTS ${LIB_SOLVER} OR NOT EXISTS ${LIB_BLAS} OR NOT EXISTS ${LIB_LAPACK} OR NOT EXISTS ${LIB_METIS})
        message(
            FATAL_ERROR
            "Failed to locate one of the following dependencies for Schur Complement Method:\n"
            "MA57   ${LIB_SOLVER}\n"
            "BLAS   ${LIB_BLAS}\n"
            "LAPACK ${LIB_LAPACK}\n"
            "METIS  ${LIB_METIS}"
        )
    else()
        message(
            "Located the following dependencies for Schur Comeplement Method:\n"
            "MA57   ${LIB_SOLVER}\n"
            "BLAS   ${LIB_BLAS}\n"
            "LAPACK ${LIB_LAPACK}\n"
            "METIS  ${LIB_METIS}"
        ) 
    endif()
endif()

set(
    qpOASES_CPP_FLAGS
    "-w -O3 -fPIC -DLINUX -D${DEF_SOLVER} -D__USE_LONG_FINTS__ -D__NO_COPYRIGHT__"
)

set(
    qpOASES_MAKE_ARGS
    ${qpOASES_MAKE_ARGS} CPPFLAGS=${qpOASES_CPP_FLAGS}
)

ExternalProject_Add(
    qpoases
    DOWNLOAD_COMMAND cp -a ${CMAKE_SOURCE_DIR}/external/qpOASES/. ${CMAKE_BINARY_DIR}/external/src/qpoases
    PREFIX external
    CONFIGURE_COMMAND mkdir -p bin;
    BUILD_IN_SOURCE ON
    BUILD_COMMAND make ${qpOASES_MAKE_ARGS} src;
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    INSTALL_COMMAND cp ${CMAKE_BINARY_DIR}/external/src/qpoases/bin/libqpOASES.so ${CMAKE_INSTALL_PREFIX}/lib/libqpOASES.so
)

add_library(qpoases_lib SHARED IMPORTED)
set_target_properties(qpoases_lib PROPERTIES
  IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}qpOASES${CMAKE_SHARED_LIBRARY_SUFFIX}
  IMPORTED_NO_SONAME TRUE
)
set(qpoases_include "${CMAKE_BINARY_DIR}/external/src/qpoases/include")

# 2) OSQP
ExternalProject_Add(
    osqp
    DOWNLOAD_COMMAND cp -a ${CMAKE_SOURCE_DIR}/external/osqp/. ${CMAKE_BINARY_DIR}/external/src/osqp
    PREFIX external
    CMAKE_ARGS
        -DDLONG=OFF
        -DBUILD_SHARED_LIBS=ON
    BUILD_COMMAND cmake --build .
    INSTALL_COMMAND cp ${CMAKE_BINARY_DIR}/external/src/osqp-build/out/libosqp.so ${CMAKE_INSTALL_PREFIX}/lib/libosqp.so
)

add_library(osqp_lib SHARED IMPORTED)
set_target_properties(osqp_lib PROPERTIES
  IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}osqp${CMAKE_SHARED_LIBRARY_SUFFIX}
)
set(osqp_include "${CMAKE_BINARY_DIR}/external/src/osqp/include")

# 3) googletest
ExternalProject_Add(
    gtest
    DOWNLOAD_COMMAND cp -a ${CMAKE_SOURCE_DIR}/external/googletest/. ${CMAKE_BINARY_DIR}/external/src/gtest
    PREFIX external
    INSTALL_COMMAND cp ${CMAKE_BINARY_DIR}/external/src/gtest-build/lib/libgtest.a ${CMAKE_INSTALL_PREFIX}/lib/
)

set(gtest_lib "-L${CMAKE_BINARY_DIR}/external/src/gtest-build/out -lgtest")
set(gtest_include "${CMAKE_BINARY_DIR}/external/src/gtest/googletest/include")

## CREATE LCQPow libraries ----------------------------------------------------------
# Include project headers
include_directories(${PROJECT_NAME} PRIVATE include)

# Add compile options
add_compile_options(
    -D__USE_LONG_FINTS__
    -D${DEF_SOLVER}
)

if (${PROFILING})
    add_compile_options(-pg)
endif()

# Save auxiliar source files to variable
aux_source_directory(src SRC_FILES)

# create static lib
add_library(${PROJECT_NAME}-static STATIC ${SRC_FILES})
set_target_properties(
    ${PROJECT_NAME}-static
    PROPERTIES
    OUTPUT_NAME ${PROJECT_NAME}
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

if (${QPOASES_SCHUR})
    target_link_libraries(
        ${PROJECT_NAME}-static
        PRIVATE ${Matlab_LIBRARIES}
    )
endif()

# create shared lib
add_library(${PROJECT_NAME}-shared SHARED ${SRC_FILES})
set_target_properties(
    ${PROJECT_NAME}-shared
    PROPERTIES
    OUTPUT_NAME ${PROJECT_NAME}
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

target_include_directories(${PROJECT_NAME}-shared SYSTEM
    PUBLIC ${qpoases_include}
    PUBLIC ${osqp_include}
)
target_link_libraries(${PROJECT_NAME}-shared
    PUBLIC qpoases_lib
    PUBLIC osqp_lib
)

if (${QPOASES_SCHUR})
    target_link_libraries(${PROJECT_NAME}-shared
        PRIVATE ${Matlab_LIBRARIES}
    )

    get_filename_component(LIB_SOLVER_NAME ${LIB_SOLVER} NAME)
    get_filename_component(LIB_BLAS_NAME ${LIB_BLAS} NAME)
    get_filename_component(LIB_LAPACK_NAME ${LIB_LAPACK} NAME)
    get_filename_component(LIB_METIS_NAME ${LIB_METIS} NAME)


    ADD_CUSTOM_TARGET(link_schur ALL 
        COMMAND ${CMAKE_COMMAND} -E create_symlink ${LIB_SOLVER} lib/${LIB_SOLVER_NAME}
        COMMAND ${CMAKE_COMMAND} -E create_symlink ${LIB_BLAS} lib/${LIB_BLAS_NAME}
        COMMAND ${CMAKE_COMMAND} -E create_symlink ${LIB_LAPACK} lib/${LIB_LAPACK_NAME}
        COMMAND ${CMAKE_COMMAND} -E create_symlink ${LIB_METIS} lib/${LIB_METIS_NAME}
    )

endif()

## Build examples -----------------------------------------------------------------------
if (${BUILD_EXAMPLES})
    aux_source_directory(examples EXAMPLE_FILES)

    FOREACH(ELEMENT ${EXAMPLE_FILES})
        # get filename w/o dir and extension
        get_filename_component(EXAMPLE_NAME ${ELEMENT} NAME_WE)

        # generate executable target
        add_executable(${EXAMPLE_NAME} ${ELEMENT})

        # link libraries
        target_link_libraries(
            ${EXAMPLE_NAME}
            PUBLIC ${PROJECT_NAME}-shared
            #PUBLIC ${PROJECT_NAME}-static
            #PRIVATE qpoases_lib osqp_lib
        )

        if (${QPOASES_SCHUR})
            target_link_libraries(
                ${EXAMPLE_NAME}
                PRIVATE ${Matlab_LIBRARIES}
            )
        endif()

        # specify output directory
        set_target_properties(
            ${EXAMPLE_NAME}
            PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/examples"
        )
    endforeach()
endif()

## Build Matlab interface ---------------------------------------------------------------
if (${BUILD_MATLAB_INTERFACE})

    # Ensure cmake can locate matlab and set matlab paths
    find_package(Matlab)
    if (NOT ${Matlab_FOUND})
        # send warning that the matlab interface could not be built
        message(
            WARNING
            "Unable to build Matlab interface (Matlab was not found)."
        )
    
    # if found, then build interface
    else()
        get_filename_component(
            MATLAB_LIBDIR
            ${Matlab_MEX_LIBRARY}
            DIRECTORY
        )

        if (${QPOASES_SCHUR})
            set(
                MA57_lib
                -lmwma57
            )
        endif()

        set(
            MATLAB_INTERFACE_DIR
            ${CMAKE_SOURCE_DIR}/interfaces/matlab
        )

        set(
            MATLAB_INTERFACE_NAME
            LCQPow
        )

        set(
            MATLAB_INTERFACE_DESTINATION
            ${CMAKE_BINARY_DIR}/lib
        )

        matlab_add_mex(
            NAME matlab_interface
            SRC ${MATLAB_INTERFACE_DIR}/${MATLAB_INTERFACE_NAME}.cpp
            OUTPUT_NAME ${MATLAB_INTERFACE_DESTINATION}/${MATLAB_INTERFACE_NAME}
            LINK_TO "-llcqpow -lqpOASES -losqp -lmwblas -lmwlapack -lmwma57"
        )

        target_link_directories(
            matlab_interface
            PUBLIC ${CMAKE_BINARY_DIR}/lib
        )

        target_compile_options(
            matlab_interface
            PUBLIC -D${DEF_SOLVER}
            PUBLIC -D__USE_LONG_FINTS__
            PUBLIC -D__MATLAB__
            PUBLIC -D__NO_COPYRIGHT__
            PUBLIC -D__cpluplus
            PUBLIC -O
            PUBLIC -largeArrayDims
            PUBLIC -lmwblas
            PUBLIC -lmwlapack
            PUBLIC -lmwma57
            PUBLIC -llcqpow
            PUBLIC -lqpOASES
            PUBLIC -losqp
        )

        target_include_directories(
            matlab_interface
            PRIVATE include
            SYSTEM ${osqp_include}
            SYSTEM ${qpoases_include}
        )

        file(
            COPY ${MATLAB_INTERFACE_DIR}/${MATLAB_INTERFACE_NAME}.m
            DESTINATION ${MATLAB_INTERFACE_DESTINATION}
        )
    
    endif()

endif()


## Build Python interface ---------------------------------------------------------------
if (${BUILD_PYTHON_INTERFACE})
    
    add_subdirectory(external/pybind11)
    add_subdirectory(interfaces/python/lcqpow)

    set(
        PYTHON_INTERFACE_DIR
        ${CMAKE_SOURCE_DIR}/interfaces/python
    )

    set(
        PYTHON_INTERFACE_NAME
        lcqpow
    )

    set(
        PYTHON_INTERFACE_DESTINATION
        ${CMAKE_BINARY_DIR}/interfaces/python
    )

    file(
        COPY ${PYTHON_INTERFACE_DIR}/${PYTHON_INTERFACE_NAME}/__init__.py
        DESTINATION ${PYTHON_INTERFACE_DESTINATION}/${PYTHON_INTERFACE_NAME}
    )

    file(
        COPY ${PYTHON_INTERFACE_DIR}/examples/
        DESTINATION ${PYTHON_INTERFACE_DESTINATION}
    )

    file(
        COPY ${PROJECT_SOURCE_DIR}/examples/example_data
        DESTINATION ${PYTHON_INTERFACE_DESTINATION}
    )
    
endif()

## Build unit tests ---------------------------------------------------------------------
if (${UNIT_TESTS})
    # Add include directories of dependencies: qpOASES, OSQP
    include_directories(${PROJECT_NAME}
        SYSTEM ${gtest_include}
        SYSTEM ${qpoases_include}
        SYSTEM ${osqp_include}
    )

    # Add Unit testing source file
    add_executable(RunUnitTests test/RunUnitTests.cpp)

    set_target_properties(
        RunUnitTests
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/tests"
    )

    # Add libraries: gtest, qpOASES, OSQP
    target_link_libraries(
        RunUnitTests
        PUBLIC ${PROJECT_NAME}-shared
        PRIVATE ${gtest_lib} ${LIB_SOLVER}
        PRIVATE -lpthread
    )

    # Build examples for testing
    aux_source_directory(test/examples TEST_FILES)

    FOREACH(ELEMENT ${TEST_FILES})
        # get filename w/o dir and extension
        get_filename_component(EXAMPLE_NAME ${ELEMENT} NAME_WE)

        # save to example test name
        set(EXAMPLE_TEST_NAME "${EXAMPLE_NAME}_TEST")

        # generate executable target
        add_executable(${EXAMPLE_TEST_NAME} ${ELEMENT})

        # link libraries
        target_link_libraries(
            ${EXAMPLE_TEST_NAME}
            PUBLIC ${PROJECT_NAME}-shared
        )

        if (${QPOASES_SCHUR})
            target_link_libraries(
                ${EXAMPLE_TEST_NAME}
                PRIVATE ${Matlab_LIBRARIES}
            )
        endif()

        # specify output directory
        set_target_properties(
            ${EXAMPLE_TEST_NAME}
            PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/tests"
        )
    endforeach()
endif()
