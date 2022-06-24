// The MIT License (MIT)
//
// Copyright (c) 2014 MiguelVieira
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
// Original code can be found in
// https://github.com/MiguelVieira/ConvexHull2D/blob/master/ConvexHull.cpp

#include "quick_hull.h"

using namespace std;
using namespace Eigen;

// The length of segment (a, b).
static float len(const Vector2f& a, const Vector2f& b) {
	return sqrt((b[0] - a[0]) * (b[0] - a[0]) + (b[1] - a[1]) * (b[1] - a[1]));
}

// The unsigned distance of p from segment (a, b).
static float dist(const Vector2f& a, const Vector2f& b, const Vector2f& p) {
	return fabs((b[0] - a[0]) * (a[1] - p[1]) - (b[1] - a[1]) * (a[0] - p[0]))
		/ len(a, b);
}

// Returns the index of the farthest Vector2f from segment (a, b).
static size_t getFarthest(const Vector2f& a, const Vector2f& b, const vector<Vector2f>& v) {
	size_t idxMax = 0;
	float distMax = dist(a, b, v[idxMax]);

	for (size_t i = 1; i < v.size(); ++i) {
		float distCurr = dist(a, b, v[i]);
		if (distCurr > distMax) {
			idxMax = i;
			distMax = distCurr;
		}
	}

	return idxMax;
}

// The z-value of the cross product of segments
// (a, b) and (a, c). Positive means c is ccw
// from (a, b), negative cw. Zero means its collinear.
static float ccw(const Vector2f& a, const Vector2f& b, const Vector2f& c) {
	return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
}

//Check if 3 Vector2fs are in positive (counter clockwise) orientation
// Recursive call of the quickhull algorithm.
static void quickHullRec(const vector<Vector2f>& v, const Vector2f& a, const Vector2f& b,
		  vector<Vector2f>& hull) {
	if (v.empty()) {
		return;
	}

	Vector2f f = v[getFarthest(a, b, v)];

	// Collect Vector2fs to the left of segment (a, f)
	vector<Vector2f> left;
	for (auto p : v) {
		if (ccw(a, f, p) > 0) {
			left.push_back(p);
		}
	}
	quickHullRec(left, a, f, hull);

	// Add f to the hull
	hull.push_back(f);

	// Collect Vector2fs to the left of segment (f, b)
	vector<Vector2f> right;
	for (auto p : v) {
		if (ccw(f, b, p) > 0) {
			right.push_back(p);
		}
	}
	quickHullRec(right, f, b, hull);
}
// Returns true if a is lexicographically before b.
static bool isLeftOf(const Vector2f& a, const Vector2f& b) {
	return (a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]));
}

// QuickHull algorithm.
// https://en.wikipedia.org/wiki/QuickHull
vector<Vector2f> quickHull(const vector<Vector2f>& v) {
	vector<Vector2f> hull;

	// Start with the leftmost and rightmost vertex
	Vector2f a = *min_element(v.begin(), v.end(), isLeftOf);
	Vector2f b = *max_element(v.begin(), v.end(), isLeftOf);

	// Split the vertex on either side of segment (a, b)
	vector<Vector2f> left, right;
	for (auto p: v) {
		if (p != a && p != b)
			ccw(a, b, p) > 0 ? left.push_back(p) : right.push_back(p);
	}

	// Be careful to add Vector2fs to the hull
	// in the correct order. Add our leftmost Vector2f.
	hull.push_back(a);

	// Add hull Vector2fs from the left (top)
	quickHullRec(left, a, b, hull);

	// Add our rightmost Vector2f
	hull.push_back(b);

	// Add hull Vector2fs from the right (bottom)
	quickHullRec(right, b, a, hull);

	// Result is clockwise, we need it counter-clockwise
	reverse(hull.begin(), hull.end());
	return hull;
}
