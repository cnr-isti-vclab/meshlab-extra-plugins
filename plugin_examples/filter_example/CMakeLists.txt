# Copyright 2019, 2020, Visual Computing Lab, ISTI - Italian National Research Council

# Simple CMake example of a MeshLab Plugin

# You just need to call the function "add_meshlab_plugin", passing
# as arguments the name of the plugin and the list of headers and sources.

# The name of the plugin will be the name of the target, therefore you can
# use that name if you need to link the plugin to some other libraries.

set(SOURCES filter_example.cpp)

set(HEADERS filter_example.h)

add_meshlab_plugin(filter_example ${SOURCES} ${HEADERS})
