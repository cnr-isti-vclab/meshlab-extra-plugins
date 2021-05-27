# Copyright 2019, 2020, Collabora, Ltd.
# Copyright 2019, 2020, Visual Computing Lab, ISTI - Italian National Research Council
# SPDX-License-Identifier: BSL-1.0

message(STATUS "Searching for optional components")

# gmp or mpir - optional, for filter_csg
include(${CMAKE_CURRENT_SOURCE_DIR}/gmp-mpir.cmake)

# opengr - optional and header-only, for filter_globalregistration
include(${CMAKE_CURRENT_SOURCE_DIR}/opengr.cmake)
