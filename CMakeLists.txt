# Copyright 2019, 2020, Collabora, Ltd.
# Copyright 2019, 2020, Visual Computing Lab, ISTI - Italian National Research Council
# SPDX-License-Identifier: BSL-1.0

cmake_minimum_required(VERSION 3.18)
project(MeshLabExtraPlugins)

set (MESHLAB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/meshlab)
set (MESHLAB_SRC_DIR ${MESHLAB_DIR}/src)
include(${MESHLAB_SRC_DIR}/cmake/meshlab_global_settings.cmake)

### Build options
option(BUILD_PLUGIN_EXAMPLES "Build the plugin examples" OFF)

### Dependencies
list(APPEND CMAKE_MODULE_PATH "${MESHLAB_SRC_DIR}/cmake")
message(STATUS "Searching for required components")
find_package(OpenGL REQUIRED)

find_package(
	Qt5
	COMPONENTS OpenGL Xml Network
	REQUIRED)
find_package(OpenMP)

message(STATUS "Searching for required components with bundled fallback")
find_package(GLEW)
find_package(Eigen3)

### Enter subdirectories
set (MESHLAB_PLUGINS "") #no meshlab plugins needed
set (MESHLAB_USE_DEFAULT_BUILD_AND_INSTALL_DIRS ON)
set (MESHLAB_BUILD_ONLY_LIBRARIES ON)
set (MESHLAB_ALLOW_OPTIONAL_EXTERNAL_LIBRARIES OFF)
add_subdirectory(${MESHLAB_DIR})
set (MESHLAB_PLUGIN_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/meshlab/plugins)
message(STATUS "Plugins output dir: " ${MESHLAB_PLUGIN_OUTPUT_DIR})


# External
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external)

# Clear RPATH before adding plugins
set(CMAKE_INSTALL_RPATH)

### Plugin subdirectories
set (MESHLAB_PLUGINS
	plugins/filter_csg
	plugins/filter_globalregistration
	plugins/filter_orientedbbox
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
