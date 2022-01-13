# Copyright(C) 2022, Alfonso Sánchez-Beato
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

option(BUILD_TESTS "Set me to enable building test_hybbrid binary")

set(SOURCES filter_orientedbbox.cpp genetic.cpp quick_hull.cpp rotating_calipers.cpp)

set(HEADERS filter_orientedbbox.h genetic.h quick_hull.h rotating_calipers.h)

add_meshlab_plugin(filter_orientedbbox ${SOURCES} ${HEADERS})

if(BUILD_TESTS)
  add_executable(test_hybbrid test_hybbrid.cpp genetic.cpp)
  target_include_directories(test_hybbrid PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(test_hybbrid PUBLIC meshlab-common)
endif()
