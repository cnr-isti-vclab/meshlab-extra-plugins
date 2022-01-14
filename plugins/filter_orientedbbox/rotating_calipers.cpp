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

#include "quick_hull.h"
#include "rotating_calipers.h"

using namespace std;
using namespace Eigen;

template <int C>
static bool cmp_v_coord(const Vector2f& v1, const Vector2f& v2)
{
	return v1[C] < v2[C];
}

// Returns 2D lines intersection. Lines are represented as coefficients:
// li[0]*x + li[1]*y + li[2] = 0
static Vector2f intersectLines(const Vector3f& l1, const Vector3f& l2)
{
	Vector2f pt;
	float denom = l1[1]*l2[0] - l1[0]*l2[1];
	pt[0] = (l1[2]*l2[1] - l1[1]*l2[2])/denom;
	pt[1] = (l1[0]*l2[2] - l1[2]*l2[0])/denom;
	return pt;
}

static inline size_t getPrevCircularQ(size_t i, size_t size)
{
	return i == 0? size - 1: i - 1;
}

template <int C>
static size_t getFirstDifferentCounterClock(const vector<Vector2f> &vert, size_t start)
{
	size_t vert_sz = vert.size();
	size_t idx;
	for (idx = getPrevCircularQ(start, vert_sz);
	     vert[idx][C] == vert[start][C];
	     idx = getPrevCircularQ(idx, vert_sz));

        return (idx + 1) % vert_sz;
}

// Implements the 2D rotating calipers algorithm [1].
//
// [1] Toussaint, G. T. (1983, May). Solving geometric problems with
// the rotating calipers. In Proc. IEEE Melecon (Vol. 83, p. A10).
RotAreaPair rotatingCalipers(const std::vector<Eigen::Vector2f> &vert)
{
	// 4 quadrants: top-right, top-left, bottom-left, bottom-right
	static constexpr size_t NSIDES = 4;
	vector<Vector2f> hull;
	double rot = 0.;
	double vol = INFINITY, minVolAng = NAN;
	// Angles order: right, up, left, down
	vector<double> angle(NSIDES);
	// Index to current hull vertex for each side
	vector<size_t> idx(NSIDES), idxNext(NSIDES);
	vector<Vector2f> vecNext(NSIDES);

	// We get convex hull with Vector2fs in counterclockwise order
	hull = quickHull(vert);
	size_t hull_sz = hull.size();
	auto x_minmax = minmax_element(hull.begin(), hull.end(), cmp_v_coord<0>);
	auto y_minmax = minmax_element(hull.begin(), hull.end(), cmp_v_coord<1>);
	// More than one max or min? which one to take? (repeated indexes)
        // -> take the one that appears first when running counterclockwise
	idx[2] = x_minmax.first - hull.begin();
	idx[2] = getFirstDifferentCounterClock<0>(hull, idx[2]);
	idx[0] = x_minmax.second - hull.begin();
	idx[0] = getFirstDifferentCounterClock<0>(hull, idx[0]);
	idx[3] = y_minmax.first - hull.begin();
	idx[3] = getFirstDifferentCounterClock<1>(hull, idx[3]);
	idx[1] = y_minmax.second - hull.begin();
	idx[1] = getFirstDifferentCounterClock<1>(hull, idx[1]);
	for (size_t i = 0; i < NSIDES; ++i) {
		idxNext[i] = (idx[i] + 1) % hull_sz;
		vecNext[i] = hull[idxNext[i]] - hull[idx[i]];
	}
	angle[0] = acos(vecNext[0][1]/vecNext[0].norm());
	angle[1] = acos(-vecNext[1][0]/vecNext[1].norm());
	angle[2] = acos(-vecNext[2][1]/vecNext[2].norm());
	angle[3] = acos(vecNext[3][0]/vecNext[3].norm());

	// TODO We could also do a for with the number of sides + 1
        while (2*rot < M_PI) {
		// Unit vectors for right and up side of square
		Vector2f dir{-sin(rot), cos(rot)};
		Vector2f norm{cos(rot), sin(rot)};
		// Taking dir as normal and R as point we get a line from right to left
		Vector3f lineWithRHorz{dir[0], dir[1],
				       - dir[0]*hull[idx[0]][0] - dir[1]*hull[idx[0]][1]};
		// Now up from down line passing through L
		Vector3f lineWithLVert{norm[0], norm[1],
				       - norm[0]*hull[idx[2]][0] - norm[1]*hull[idx[2]][1]};
		Vector2f orthL = intersectLines(lineWithRHorz, lineWithLVert);
		Vector3f lineWithUHorz{dir[0], dir[1],
				       - dir[0]*hull[idx[1]][0] - dir[1]*hull[idx[1]][1]};
		Vector3f lineWithDVert{norm[0], norm[1],
				       - norm[0]*hull[idx[3]][0] - norm[1]*hull[idx[3]][1]};
		Vector2f orthD = intersectLines(lineWithUHorz, lineWithDVert);
		double curVol = (hull[idx[0]] - orthL).norm()*(hull[idx[3]] - orthD).norm();
		if (curVol < vol) {
			vol = curVol;
			minVolAng = rot;
		}

		auto minAngIt = min_element(angle.begin(), angle.end());
		size_t idxMinAng = minAngIt - angle.begin();
		rot = angle[idxMinAng];
		idx[idxMinAng] = idxNext[idxMinAng];
		idxNext[idxMinAng] = (idx[idxMinAng] + 1) % hull_sz;
		vecNext[idxMinAng] = hull[idxNext[idxMinAng]] - hull[idx[idxMinAng]];
		switch (idxMinAng) {
		case 0:
			angle[0] = acos(vecNext[0][1]/vecNext[0].norm());
			break;
		case 1:
			angle[1] = acos(-vecNext[1][0]/vecNext[1].norm());
			break;
		case 2:
			angle[2] = acos(-vecNext[2][1]/vecNext[2].norm());
			break;
		case 3:
			angle[3] = acos(vecNext[3][0]/vecNext[3].norm());
			break;
		};
	}

	return RotAreaPair(minVolAng, vol);
}
