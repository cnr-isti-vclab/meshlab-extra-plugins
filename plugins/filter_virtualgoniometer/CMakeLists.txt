# Copyright 2019 Collabora, Ltd.
# SPDX-License-Identifier: BSL-1.0

### Generated file! Edit the templates in src/templates then re-run ./make-cmake.py

set(SOURCES virtual_goniometer.cpp)

set(HEADERS virtual_goniometer.h)

set(RESOURCES meshlab.qrc)

add_meshlab_plugin(filter_virtualgoniometer ${SOURCES} ${HEADERS} ${RESOURCES})
