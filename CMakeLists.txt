# Copyright 2019, 2020, Collabora, Ltd.
# Copyright 2019, 2020, Visual Computing Lab, ISTI - Italian National Research Council
# SPDX-License-Identifier: BSL-1.0

cmake_minimum_required(VERSION 3.9)
project(MeshLabExtraPlugins)

set (MESHLAB_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/meshlab/src)
include(${MESHLAB_SRC_DIR}/cmake/meshlab_global_settings.cmake)
include(${MESHLAB_SRC_DIR}/cmake/meshlab_tools.cmake)

### Build options
option(BUILD_STRICT "Strictly enforce resolution of all symbols" ON)
option(BUILD_WITH_DOUBLE_SCALAR "Use double type instead of float type for scalars" OFF)
option(BUILD_PLUGIN_EXAMPLES "Build the plugin examples" OFF)

### Dependencies
#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
message(STATUS "Searching for required components")
find_package(OpenGL REQUIRED)

find_package(
	Qt5
	COMPONENTS OpenGL Xml Network
	REQUIRED)
find_package(OpenMP)
find_package(Qt5Qml)

message(STATUS "Searching for required components with bundled fallback")
find_package(GLEW)
find_package(Eigen3)

### Build directories
set(MESHLAB_BUILD_DISTRIB_DIR ${CMAKE_CURRENT_BINARY_DIR}/distrib)
set(MESHLAB_LIB_OUTPUT_DIR ${MESHLAB_BUILD_DISTRIB_DIR})
set(MESHLAB_PLUGIN_OUTPUT_DIR ${MESHLAB_BUILD_DISTRIB_DIR}/plugins)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${MESHLAB_BUILD_DISTRIB_DIR})

### Install directories
if(WIN32 OR APPLE)
	set(INSTALL_TO_UNIX_LAYOUT OFF)
else()
	set(INSTALL_TO_UNIX_LAYOUT ON)
endif()
include(GNUInstallDirs)
if(INSTALL_TO_UNIX_LAYOUT)
	set(MESHLAB_LIB_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR}/meshlab)
	set(MESHLAB_PLUGIN_INSTALL_DIR ${MESHLAB_LIB_INSTALL_DIR}/plugins)
elseif(APPLE)
	set(MESHLAB_LIB_INSTALL_DIR meshlab.app/Contents/Frameworks)
	set(MESHLAB_PLUGIN_INSTALL_DIR meshlab.app/Contents/PlugIns ${MESHLAB_PARENT_SCOPE})
else()
	set(MESHLAB_LIB_INSTALL_DIR .)
	set(MESHLAB_PLUGIN_INSTALL_DIR plugins ${MESHLAB_PARENT_SCOPE})
endif()

### Install Settings
if (NOT APPLE)
	set(CMAKE_INSTALL_RPATH $ORIGIN/../${MESHLAB_LIB_INSTALL_DIR};$ORIGIN/../${CMAKE_INSTALL_LIBDIR})
else()
	SET(CMAKE_INSTALL_RPATH $ORIGIN/../Frameworks)
endif()

### Enter subdirectories

# VCGLib -- required
if(NOT VCGDIR)
	get_filename_component(VCGDIR "${MESHLAB_SRC_DIR}/vcglib" ABSOLUTE)
	if(NOT EXISTS ${VCGDIR})
		set(VCGDIR NOTFOUND)
	endif()
endif()
set(VCGDIR "${VCGDIR}")

if(NOT VCGDIR)
	message(FATAL_ERROR "VCGLib not found. Please clone recursively the MeshLab repo.")
endif()


# External
set(EXTERNAL_DIR ${MESHLAB_SRC_DIR}/external)
include(${EXTERNAL_DIR}/external_common.cmake)
set(EXTRA_EXTERNAL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external)
add_subdirectory(${EXTRA_EXTERNAL_DIR})

add_subdirectory(${MESHLAB_SRC_DIR}/common)

# Clear RPATH before adding plugins
set(CMAKE_INSTALL_RPATH)

### Plugin subdirectories
set (MESHLAB_PLUGINS
	plugins/filter_globalregistration
	plugins/filter_poisson
	plugins/filter_virtualgoniometer
)

set (MESHLAB_PLUGIN_EXAMPLES
	plugin_examples/iomesh_example
	plugin_examples/filter_example
)

message(STATUS "\nConfiguring plugins")
foreach(PLUGIN ${MESHLAB_PLUGINS})
	if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${PLUGIN}/CMakeLists.txt)
		message(STATUS "- ${PLUGIN}")
		add_subdirectory(${PLUGIN})
	else()
		message(STATUS " - ${PLUGIN} - Skipping, plugin or build system not found.")
	endif()
endforeach()

if (BUILD_PLUGIN_EXAMPLES)
	message(STATUS "\nConfiguring plugin examples")
	foreach(PLUGIN ${MESHLAB_PLUGIN_EXAMPLES})
		if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${PLUGIN}/CMakeLists.txt)
			message(STATUS "- ${PLUGIN}")
			add_subdirectory(${PLUGIN})
		else()
			message(STATUS " - ${PLUGIN} - Skipping, plugin or build system not found.")
		endif()
	endforeach()
endif()
