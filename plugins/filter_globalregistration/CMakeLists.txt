# Copyright 2019-2020, Collabora, Ltd.
# SPDX-License-Identifier: BSL-1.0


# Only build if we have opengr
if(TARGET external-opengr)

	set(SOURCES globalregistration.cpp)

	set(HEADERS globalregistration.h)

	add_meshlab_plugin(filter_globalregistration ${SOURCES} ${HEADERS})

	target_link_libraries(filter_globalregistration PRIVATE external-opengr)
	if(OpenMP_CXX_FOUND)
		target_link_libraries(filter_globalregistration PRIVATE OpenMP::OpenMP_CXX)
	endif()

else()
	message(STATUS "Skipping filter_globalregistration - missing opengr in external directory.")
endif()
