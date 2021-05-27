# Copyright 2019, 2020, Collabora, Ltd.
# Copyright 2019, 2020, Visual Computing Lab, ISTI - Italian National Research Council
# SPDX-License-Identifier: BSL-1.0

cmake_minimum_required(VERSION 3.13)
project(MeshLabExtraPlugins)

set (MESHLAB_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/meshlab/src)
include(${MESHLAB_SRC_DIR}/cmake/meshlab_global_settings.cmake)
include(${MESHLAB_SRC_DIR}/cmake/meshlab_tools.cmake)

### Build options
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

message(STATUS "Searching for required components with bundled fallback")
find_package(GLEW)
find_package(Eigen3)

### Enter subdirectories
set (MESHLAB_PLUGINS "") #no meshlab plugins needed
set (BUILD_ONLY_MESHLAB_LIBRARIES ON)
set (ALLOW_OPTIONAL_EXTERNAL_MESHLAB_LIBRARIES OFF)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
add_subdirectory(${MESHLAB_SRC_DIR})
message(STATUS "Plugins output dir: " ${MESHLAB_PLUGIN_OUTPUT_DIR})


# External
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external)

# Clear RPATH before adding plugins
set(CMAKE_INSTALL_RPATH)

### Plugin subdirectories
set (MESHLAB_PLUGINS
	plugins/filter_csg
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
