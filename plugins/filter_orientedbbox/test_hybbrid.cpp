/*
 * Copyright (c) 2022 Alfonso Sánchez-Beato
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

#include <cassert>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <cstddef>
#include <iostream>
#include <vector>

// Include all algorithm so we can test static functions
#include "Eigen/src/Core/Matrix.h"
#include "quick_hull.cpp"
#include "rotating_calipers.cpp"
#include "genetic.h"

using namespace std;
using namespace Eigen;

static bool cmp_float(float a, float b, float eps=0.001f)
{
	return fabsf(a - b) < eps;
}

template <typename M>
static bool cmp_matrix(const M &a, const M &b, float eps=0.001f)
{
	assert(a.size() == b.size());
	for (size_t i = 0, sz = a.size(); i < sz; ++i) {
		if (fabsf(a(i) - b(i)) >= eps)
			return false;
	}

        return true;
}

struct TestGeneticAlgorithm {
        TestGeneticAlgorithm(GeneticAlgorithm &ga): ga_(ga) {}
	GeneticAlgorithm &ga_;

        float bbVolume(const Eigen::Matrix3f &rot) const {
		return ga_.bbVolume(rot);
	}
	Matrix3f rotMatricesMean(Matrix3f *rot_v, size_t rot_n) const {
		return ga_.rotMatricesMean(rot_v, rot_n);
	}
	struct GeneticAlgorithm::RotVolPair localOptimum(Matrix3f rot) const {
		return ga_.localOptimum(rot);
	}
};

static void test_bbVolume(void)
{
	Matrix3f rot;
	//rot << 1., 0., 0., 0., 1., 0., 0., 0., 1.;
	rot = AngleAxisf(0.25*M_PI, Vector3f::UnitX())
		* AngleAxisf(0.5*M_PI,  Vector3f::UnitY())
		* AngleAxisf(0.33*M_PI, Vector3f::UnitZ());
	vector<Vector3f> vert{Vector3f{1., 0., 0.},
			      Vector3f{0., 0., 1.},
			      Vector3f{0., 2., 0.}};
	vector<Vector3f> vertRot;
	for (auto &pt: vert)
		vertRot.push_back(rot*pt);

	GeneticAlgorithm ga(vertRot);
	TestGeneticAlgorithm tga(ga);

	assert(cmp_float(tga.bbVolume(Matrix3f::Identity()), 2.8398f));
	assert(cmp_float(tga.bbVolume(rot), 2.f));
}

float frobeniusNorm(const Matrix3f &m)
{
	float norm = 0;
	size_t sz = m.size();
	for (size_t i = 0; i < sz; ++i) {
		float v = m(i);
		norm += v*v;
	}
	return norm;
}

Matrix3f calcProjectionToSO3(const Matrix3f &mat)
{
	std::vector<Eigen::Matrix3f> diag_rot_mats_(3);
	std::vector<Eigen::Matrix3f> diag_reflect_mats_(4);
	diag_rot_mats_[0] <<  1, 0, 0, 0, -1, 0, 0, 0, -1;
	diag_rot_mats_[1] << -1, 0, 0, 0,  1, 0, 0, 0, -1;
	diag_rot_mats_[2] << -1, 0, 0, 0, -1, 0, 0, 0,  1;
	diag_reflect_mats_[0] <<  1, 0, 0, 0,  1, 0, 0, 0, -1;
	diag_reflect_mats_[1] <<  1, 0, 0, 0, -1, 0, 0, 0,  1;
	diag_reflect_mats_[2] << -1, 0, 0, 0,  1, 0, 0, 0,  1;
	diag_reflect_mats_[3] << -1, 0, 0, 0, -1, 0, 0, 0, -1;

        HouseholderQR<Matrix3f> qr(mat);
	Matrix3f q = qr.householderQ();
	Matrix3f r = qr.matrixQR().triangularView<Upper>();
	cout << "Project:\n" << mat << endl;
	cout << "det(mat):\t" << mat.determinant() << endl;
	const vector<Matrix3f> *test_mats;
	Matrix3f min_q, min_r;
	float min_dist = INFINITY;

	if (q.determinant() > 0) {
		// q is a rotation and a possible solution already, we'll
		// check also by multiplying by diagonal rotation matrices.
		test_mats = &diag_rot_mats_;
		min_q = q;
		min_r = r;
		min_dist = frobeniusNorm(mat - q);
	} else {
		// As q is a reflection, multiplying by other reflections will
		// create rotation matrices, and we will choose the nearest to
		// the input matrix.
		test_mats = &diag_reflect_mats_;
	}

	// Find best Q of all possible decompositions. The test matrices fulfill
	// q*d*d^T*r=q*r, while q*d and d^T*r are still orthogonal and
	// upper-triangular respectively.
	size_t num_diag = test_mats->size();
	for (size_t i = 0; i < num_diag; ++i) {
		Matrix3f q_i = q*(*test_mats)[i];
		float dist = frobeniusNorm(mat - q_i);
		if (dist < min_dist) {
			min_dist = dist;
			min_q = q_i;
			min_r = (*test_mats)[i].transpose()*r;
		}
	}

	cout << "Q:\n" << min_q << endl;
	cout << "R:\n" << min_r << endl;
	return min_q;
}

// Mean of an array of rotation matrices.
Matrix3f rotMatricesMean(const Matrix3f *rot_v, size_t rot_n)
{
	Matrix3f mean(Matrix3f::Zero());

	for (size_t r_i = 0; r_i < rot_n; ++r_i) {
		mean += rot_v[r_i];
	}
	mean /= rot_n;

	return calcProjectionToSO3(mean);
}

static void test_rotMatricesMean(void)
{
        Matrix3f rot_v[3];
	// identity, no rotation
	rot_v[0] = Matrix3f::Identity();
	rot_v[1] = Matrix3f();
	// pi/2 with x as rotation axis
	rot_v[1] << 1, 0, 0, 0, 0, -1, 0, 1, 0;
	Matrix3f mean = rotMatricesMean(rot_v, 2);
	Matrix3f mean_res;
	// Mean is pi/4 around x axis
	mean_res << 1,         0,         0,
		    0,  0.707107, -0.707107,
		    0,  0.707107,  0.707107;
	assert(cmp_matrix<Matrix3f>(mean, mean_res));
	assert(cmp_matrix<Matrix3f>(mean*mean.transpose(), Matrix3f::Identity()));

	// Including the own mean needs to produce the mean again
	rot_v[2] = mean_res;
	mean = rotMatricesMean(rot_v, 3);
	assert(cmp_matrix<Matrix3f>(mean, mean_res));
	assert(cmp_matrix<Matrix3f>(mean*mean.transpose(), Matrix3f::Identity()));

	rot_v[0] = Matrix3f::Identity();
	// pi/6 around x axis
	rot_v[1] << 1, 0, 0, 0, 0.8660254, -0.5, 0, 0.5, 0.8660254;
	mean = rotMatricesMean(rot_v, 2);
	// Mean is pi/12 around x axis
	mean_res << 1,         0,         0,
		    0,  0.965926, -0.258819,
		    0,  0.258819,  0.965926;
	assert(cmp_matrix<Matrix3f>(mean, mean_res));
	assert(cmp_matrix<Matrix3f>(mean*mean.transpose(), Matrix3f::Identity()));

	rot_v[0] = Matrix3f::Identity();
	// pi/2 with y as rotation axis
	rot_v[1] <<  0, 0, 1,
		     0, 1, 0,
		    -1, 0, 0;
	mean = rotMatricesMean(rot_v, 2);
	mean_res <<  0.707107,         0,  0.707107,
		     0,                1,         0,
		    -0.707107,         0,  0.707107;
	assert(cmp_matrix<Matrix3f>(mean, mean_res));
	assert(cmp_matrix<Matrix3f>(mean*mean.transpose(), Matrix3f::Identity()));

	rot_v[0] = Matrix3f::Identity();
	// pi/2 with z as rotation axis
	rot_v[1] <<  0, -1, 0,
		     1,  0, 0,
		     0,  0, 0;
	mean = rotMatricesMean(rot_v, 2);
	mean_res <<  0.707107,        -0.707107,  0,
		     0.707107,         0.707107,  0,
		     0,                0,         1;
	assert(cmp_matrix<Matrix3f>(mean, mean_res));
	assert(cmp_matrix<Matrix3f>(mean*mean.transpose(), Matrix3f::Identity()));

	rot_v[0] = Matrix3f::Identity();
	// Two rotations - this is ecpected to fail if in the future the
	// projection is calculated by a polar decomposition.
	rot_v[1] = AngleAxisf(0.25*M_PI, Vector3f::UnitX())
		* AngleAxisf(0.25*M_PI, Vector3f::UnitY());
	cout << rot_v[1] << endl;
	mean = rotMatricesMean(rot_v, 2);
	mean_res <<  0.92388,  -0.136774,  0.357407,
		     0.270598,  0.893889, -0.357407,
		    -0.270598,  0.426914,  0.862856;
	assert(cmp_matrix<Matrix3f>(mean, mean_res));
	assert(cmp_matrix<Matrix3f>(mean*mean.transpose(), Matrix3f::Identity()));
}

static void checkVec2f(const Vector2f v1, const Vector2f v2)
{
	for (long i = 0; i < v1.size(); ++i)
		assert(cmp_float(v1[i], v2[i]));
}

static void checkLineCross(Vector3f l1, Vector3f l2, Vector2f cross)
{
	Vector2f pt = intersectLines(l1, l2);
	assert(cmp_float(pt[0], cross[0]));
	assert(cmp_float(pt[1], cross[1]));
}

static void test_intersectLines(void)
{
	checkLineCross({1.f, 1.f, -1.f}, {-1.f, 1.f, 0.f}, {.5f, .5f});
	checkLineCross({0.f, 1.f, 0.f}, {1.f, 0.f, 0.f}, {.0f, .0f});
	checkLineCross({1.f, 1.f, -1.f}, {0.f, 1.f, 0.f}, {1.f, .0f});
	// Parallel lines - should not happen in this algorithm
	//checkLineCross({0.f, 1.f, 0.f}, {0.f, 1.f, 1.f}, {.0f, .0f});
}

static void test_quickHull(void)
{
	// Try a grid
	{
		vector<Vector2f> vert;
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				vert.push_back({i, j});
		vector<Vector2f> hull = quickHull(vert);
		checkVec2f(hull[0], {2., 0.});
		checkVec2f(hull[1], {2., 2.});
		checkVec2f(hull[2], {0., 2.});
		checkVec2f(hull[3], {0., 0.});
	}

	// Check return order is the same with different input order
	{
		vector<Vector2f> vert;
		// Try a grid
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				vert.push_back({j, i});
		vector<Vector2f> hull = quickHull(vert);
		checkVec2f(hull[0], {2., 0.});
		checkVec2f(hull[1], {2., 2.});
		checkVec2f(hull[2], {0., 2.});
		checkVec2f(hull[3], {0., 0.});
	}

	{
		vector<Vector2f> vert {
			{0, 0},
			{-0.647477388, -0.762049496},
			{0.547499537, -0.458458036},
			{-0.530110657, 0.457271039}
		};
		vector<Vector2f> hull = quickHull(vert);
		assert(hull.size() == 3);
		checkVec2f(hull[0], {0.547499537, -0.458458036});
		checkVec2f(hull[1], {-0.530110657, 0.457271039});
		checkVec2f(hull[2], {-0.647477388, -0.762049496});
	}
}

static void test_rotatingCalipers(void)
{
	vector<Vector3f> vert3d;
	float d = 0.2;
	/* First test: 6 sides
	 *
	 *    |-\
	 *    \  \
	 *     \  \
	 *      \_|
	 */
	vector<vector<Vector2f>> vertTests {
		{{1.,0.}, {1.0-d,0.}, {1.,d},
		 {0.,1.}, {d,1.}, {0.,1.-d}},
		{{.25, 0.}, {1.0, 0.5},
	         {0.25, 1.}, {0., 0.}},
		{{0., 1.}, {.5, 0.}, {1., 0.}}
	};
	vector<RotAreaPair> raTests {
		{.angle = 0.785398f, .area = 0.400000f},
		{.angle = 0.982794f, .area = 0.875000f},
		{.angle = 0.785398f, .area = 0.500000f},
	};
	size_t numTsts = vertTests.size();
	for (size_t i = 0; i < numTsts; ++i) {
		RotAreaPair ra = rotatingCalipers(vertTests[i]);
		printf("angle %f, area %f\n", ra.angle, ra.area);
		assert(cmp_float(ra.angle, raTests[i].angle));
		assert(cmp_float(ra.area, raTests[i].area));
	}
}

static void test_localOptimum1(void)
{
	// Cube of size 1
    	vector<Vector3f> vert {
		{0,        0,        0},
		{0.866025, 0.5,      0},
		{    -0.5, 0.866025, 0},
		{0.258819, 0.965926, 0},
		{0,        0,        1},
		{0.866025, 0.5,      1},
		{    -0.5, 0.866025, 1},
		{0.258819, 0.965926, 1},
	};
	GeneticAlgorithm ga(vert);
	TestGeneticAlgorithm tga(ga);
	// Non-optimal box, wrong by 30º
	Matrix3f notBestRot;
	notBestRot << 1, 0, 0,
		      0, 1, 0,
		      0, 0, 1;
	GeneticAlgorithm::RotVolPair rvp = tga.localOptimum(notBestRot);
	cout << "test local optimum:\n";
	cout << rvp.om << endl;
	cout << rvp.vol << endl;
	Matrix3f bestRot;
	bestRot << 0.866025, -0.5,      0,
		   0.5,       0.866025, 0,
		   0,         0,        1;
	assert(cmp_matrix<Matrix3f>(bestRot, rvp.om));
	assert(cmp_float(rvp.vol, 1));
}

static void test_localOptimum2(void)
{
	// Cube of size 1
    	vector<Vector3f> vert {
		{0,0,0},
		{1,0,0},
		{0,1,0},
		{1,1,0},
		{0,0,1},
		{1,0,1},
		{0,1,1},
		{1,1,1},
	};
	GeneticAlgorithm ga(vert);
	TestGeneticAlgorithm tga(ga);
	// Non-optimal box, wrong by 30º
	Matrix3f notBestRot;
	notBestRot << 1,  0,         0,
		      0,  0.866025,  0.5,
		      0, -0.5,       0.866025;
	GeneticAlgorithm::RotVolPair rvp = tga.localOptimum(notBestRot);
	cout << "test local optimum:\n";
	cout << rvp.om << endl;
	cout << rvp.vol << endl;
	Matrix3f bestRot;
	bestRot << 1, 0, 0,
		   0, 1, 0,
		   0, 0, 1;
	assert(cmp_matrix<Matrix3f>(bestRot, rvp.om));
	assert(cmp_float(rvp.vol, 1));
}

static void test_rotatingCalipers2(void)
{
    	vector<Vector2f> vert {
		{0,        0,      },
		{0.866025, 0.5,    },
		{    -0.5, 0.866025},
		{0.258819, 0.965926},
		{0,        0,      },
		{0.866025, 0.5,    },
		{    -0.5, 0.866025},
		{0.258819, 0.965926},
	};
	RotAreaPair ra = rotatingCalipers(vert);
	cout << "rotatingCalipers2:" << endl;
	cout << ra.angle << endl;
	cout << ra.area << endl;
	assert(cmp_float(ra.angle, M_PI/6));
	assert(cmp_float(ra.area, 1));
}

static void test_genetic(void)
{
	vector<Vector3f> vert {
		{0,0,0},
		{1,0,0},
		{0,1,0},
		{0,0,1},
	};
	GeneticAlgorithm ga(vert);

	GeneticAlgorithm::RotVolPair rvp = ga.run(NULL, NULL);
	cout << "Optimal box:\n";
	cout << rvp.om << endl;
	cout << rvp.vol << endl;
}

static void test_genetic2(void)
{
	vector<Vector3f> vert {
		{0,0,0},
		{1,0,0},
		{0,1,0},
		{1,1,0},
		{0,0,1},
		{1,0,1},
		{0,1,1},
		{1,1,1},
	};
	GeneticAlgorithm ga(vert);

	GeneticAlgorithm::RotVolPair rvp = ga.run(NULL, NULL);
	cout << "Optimal box:\n";
	cout << rvp.om << endl;
	cout << rvp.vol << endl;
}

int main(void)
{
	test_intersectLines();
	test_quickHull();
	test_rotatingCalipers();
	test_rotatingCalipers2();
	test_bbVolume();
	test_rotMatricesMean();
	test_localOptimum1();
	test_localOptimum2();
	test_genetic();
	test_genetic2();
	return 0;
}
