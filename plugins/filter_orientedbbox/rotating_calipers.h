/*
 * Copyright (c) 2022 Alfonso Sánchez-Beato
 * Based on matlab implementation by Chia-Tche Chang, Bastien Gorissen,
 * and Samuel Melchior (https://github.com/chadogome/OptimalOBB).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ROTATING_CALIPERS_H
#define ROTATING_CALIPERS_H

#include <vector>
#include <Eigen/Dense>

// Angle and resulting area for the set of 2d vertices passed to the
// rotatingCalipers function.
struct RotAreaPair {
	float angle;
	float area;
};

RotAreaPair rotatingCalipers(const std::vector<Eigen::Vector2f> &vert);

#endif // ROTATING_CALIPERS_H
