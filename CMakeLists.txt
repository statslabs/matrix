#
# Copyright 2018 The StatsLabs Authors.
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
cmake_minimum_required(VERSION 3.0)
project(matrix VERSION 0.3.0 LANGUAGES CXX)

    set(CMAKE_CXX_STANDARD 11)

if("${CMAKE_BUILD_TYPE}" STREQUAL "")
    message(STATUS "CMAKE_BUILD_TYPE is unset, defaulting to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

include("cmake/platform.cmake")
include("cmake/OpenMP.cmake")
include("cmake/SDL.cmake")
include("cmake/MKL.cmake")
include("cmake/Doxygen.cmake")
include("cmake/profiling.cmake")

##################################################
# Create target and set properties
##################################################

add_executable(main
        src/main.cpp)

target_link_libraries(main PUBLIC ${MKL_LINKER_LIBS})

# Set target properties
target_include_directories(main
        PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        )

add_library(matrix
        src/library.cpp)

target_link_libraries(matrix PUBLIC ${MKL_LINKER_LIBS})

# Add an alias so that library can be used inside the build tree, e.g. when testing
add_library(Matrix::matrix ALIAS matrix)

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
        Matrix::
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

export(EXPORT matrix-targets FILE ${CMAKE_CURRENT_BINARY_DIR}/MatrixTargets.cmake NAMESPACE Matrix::)

# Register package in user's package registry
export(PACKAGE Matrix)

add_subdirectory(examples)
add_subdirectory(test)
