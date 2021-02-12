# Copyright 2019 Collabora, Ltd.
# SPDX-License-Identifier: BSL-1.0

### Generated file! Edit the templates in src/templates then re-run ./make-cmake.py

set(SOURCES virtual_goniometer.cpp)

set(HEADERS virtual_goniometer.h)

set(RESOURCES meshlab.qrc)

add_library(filter_virtualgoniometer MODULE ${SOURCES} ${HEADERS} ${RESOURCES})

target_include_directories(filter_virtualgoniometer PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(filter_virtualgoniometer PUBLIC meshlab-common)

set_property(TARGET filter_virtualgoniometer PROPERTY FOLDER Plugins)

set_property(TARGET filter_virtualgoniometer PROPERTY RUNTIME_OUTPUT_DIRECTORY
                                           ${MESHLAB_PLUGIN_OUTPUT_DIR})

set_property(TARGET filter_virtualgoniometer PROPERTY LIBRARY_OUTPUT_DIRECTORY
                                           ${MESHLAB_PLUGIN_OUTPUT_DIR})

install(TARGETS filter_virtualgoniometer DESTINATION ${MESHLAB_PLUGIN_INSTALL_DIR}
                COMPONENT Plugins)
