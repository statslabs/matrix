#
# Copyright 2018-2019 The Statslabs Authors.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
cmake_minimum_required(VERSION 3.6.0)

set(CMAKE_CXX_STANDARD 11)

set(PROJECT_NAME "Statslabs::Matrix")
set(PROJECT_FULL_NAME "Statslabs::Matrix C++ Linear Algebra Library")

set(MATRIX_VERSION_MAJOR 0)
set(MATRIX_VERSION_MINOR 4)
set(MATRIX_VERSION_PATCH 0)
set(PROJECT_VERSION "${MATRIX_VERSION_MAJOR}.${MATRIX_VERSION_MINOR}.${MATRIX_VERSION_PATCH}")

project(${PROJECT_NAME} VERSION ${PROJECT_VERSION} LANGUAGES CXX)
message(STATUS "${PROJECT_NAME} ${PROJECT_VERSION}")

if("${CMAKE_BUILD_TYPE}" STREQUAL "")
    message(STATUS "CMAKE_BUILD_TYPE is unset, defaulting to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

if (CMAKE_HOST_WIN32)
    message(STATUS "Windows detected")
elseif(CMAKE_HOST_APPLE)
    message(STATUS "macOS detected")
    #include("cmake/vecLib.cmake")
    include("cmake/MKL.cmake")
elseif(CMAKE_HOST_SOLARIS)
    message(STATUS, "Solaris detected")
    #include("cmake/SunPerf.cmake")
else()
    include("cmake/MKL.cmake")
endif()

if (HAVE_MKL)
    list(APPEND SLAB_LINKER_LIBS ${MKL_LINKER_LIBS})
else ()
    remove_definitions(-DUSE_MKL)
    remove_definitions(-DUSE_CBLAS)

    find_package(BLAS REQUIRED)
    message(STATUS "Found BLAS library at: ${BLAS_LIBRARIES}")
    include_directories(AFTER ${BLAS_INCLUDE_DIRS})
    list(APPEND SLAB_LINKER_LIBS ${BLAS_LIBRARIES})

    find_package(LAPACK REQUIRED)
    message(STATUS "Found LAPACK library at: ${LAPACK_LIBRARIES}")
    list(APPEND SLAB_LINKER_LIBS ${LAPACK_LIBRARIES})
endif ()


##################################################
# Create target and set properties
##################################################

add_executable(main src/main.cc)

target_link_libraries(main PUBLIC ${SLAB_LINKER_LIBS})
target_compile_features(main PRIVATE cxx_alias_templates)

if(MSVC)
    target_compile_options(main PRIVATE /W4 /WX)
else()
    target_compile_options(main PRIVATE -Wall -Wextra -pedantic)
endif()

# Set target properties
target_include_directories(main
        PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        )

add_library(matrix src/library.cc)

if (HAVE_MKL)
    target_compile_definitions(matrix PUBLIC "USE_MKL")
endif()

target_compile_features(matrix PRIVATE cxx_alias_templates)

if(MSVC)
    target_compile_options(matrix PRIVATE /W4 /WX)
else()
    target_compile_options(matrix PRIVATE -Wall -Wextra -pedantic)
endif()

target_link_libraries(matrix PUBLIC ${SLAB_LINKER_LIBS})

# Add an alias so that library can be used inside the build tree, e.g. when testing
add_library(Statslabs::matrix ALIAS matrix)

# Set target properties
target_include_directories(matrix
        PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        )

##################################################
# Installation instructions
##################################################

include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/Matrix)

install(TARGETS matrix
        EXPORT matrix-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )

#This is required so that the exported target has the name Matrix and not matrix
set_target_properties(matrix PROPERTIES EXPORT_NAME Matrix)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Export the targets to a script
install(EXPORT matrix-targets
        FILE
        MatrixTargets.cmake
        NAMESPACE
        Statslabs::
        DESTINATION
        ${INSTALL_CONFIGDIR}
        )

#Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/MatrixConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/MatrixConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/MatrixConfig.cmake
        INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
        )

#Install the config, configversion and custom find modules
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/MatrixConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/MatrixConfigVersion.cmake
        DESTINATION ${INSTALL_CONFIGDIR}
        )

export(EXPORT matrix-targets FILE ${CMAKE_CURRENT_BINARY_DIR}/MatrixTargets.cmake NAMESPACE Statslabs::)

# Register package in user's package registry
export(PACKAGE Matrix)

find_package(Doxygen)
if(DOXYGEN_FOUND AND SLAB_INCLUDE_DOCS)
    set(BUILD_DOC_DIR ${CMAKE_SOURCE_DIR}/build/doc)
    if(NOT EXISTS ${BUILD_DOC_DIR})
        file(MAKE_DIRECTORY ${BUILD_DOC_DIR})
    endif()

    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    message("Doxygen build started")
    add_custom_target(Doxygen ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM)
endif()

if(SLAB_INCLUDE_EXAMPLES)
    add_subdirectory(examples)
endif()
if(SLAB_INCLUDE_TESTS)
    add_subdirectory(test)
endif()