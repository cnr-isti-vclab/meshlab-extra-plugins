# Copyright 2019-2020, Collabora, Ltd.
# SPDX-License-Identifier: BSL-1.0


set(SOURCES sampleio.cpp)

set(HEADERS sampleio.h)

add_meshlab_plugin(iomesh_example MODULE ${SOURCES} ${HEADERS})
