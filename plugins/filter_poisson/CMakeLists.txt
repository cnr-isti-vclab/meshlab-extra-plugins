# Copyright 2019-2020, Collabora, Ltd.
# SPDX-License-Identifier: BSL-1.0


set(SOURCES filter_poisson.cpp src/Factor.cpp src/Geometry.cpp
            src/MarchingCubes.cpp src/MultiGridOctest.cpp)

set(HEADERS
    filter_poisson.h
    src/Allocator.h
    src/BinaryNode.h
    src/Factor.h
    src/FunctionData.h
    src/Geometry.h
    src/Hash.h
    src/MarchingCubes.h
    src/MultiGridOctreeData.h
    src/Octree.h
    src/PPolynomial.h
    src/PoissonParam.h
    src/Polynomial.h
    src/SparseMatrix.h
    src/Vector.h)

add_meshlab_plugin(filter_poisson ${SOURCES} ${HEADERS})
