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

#include "genetic.h"
#include "Eigen/src/Core/Matrix.h"
#include "Eigen/src/Core/util/Constants.h"
#include "rotating_calipers.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// Population size
static int M_ = 30;

GeneticAlgorithm::GeneticAlgorithm(const std::vector<Eigen::Vector3f> &vert) :
	vert_(vert), diag_rot_mats_(3), diag_reflect_mats_(4)
{
	// We do not include the identity here to avoid an unnecessary
	// multiplication while finding the best QR decomposition.
	diag_rot_mats_[0] <<  1, 0, 0, 0, -1, 0, 0, 0, -1;
	diag_rot_mats_[1] << -1, 0, 0, 0,  1, 0, 0, 0, -1;
	diag_rot_mats_[2] << -1, 0, 0, 0, -1, 0, 0, 0,  1;

	diag_reflect_mats_[0] <<  1, 0, 0, 0,  1, 0, 0, 0, -1;
	diag_reflect_mats_[1] <<  1, 0, 0, 0, -1, 0, 0, 0,  1;
	diag_reflect_mats_[2] << -1, 0, 0, 0,  1, 0, 0, 0,  1;
	diag_reflect_mats_[3] << -1, 0, 0, 0, -1, 0, 0, 0, -1;
}

GeneticAlgorithm::RotVolPair
GeneticAlgorithm::run(CallBackIter cb, void *cbData, int numIter) const
{
	// Initial population via qr factorization of random rotation matrices
	vector<Simplex> pop(M_);
	for (int i = 0; i < M_; ++i) {
		Simplex sx;
		for (int j = 0; j < DIM + 1; ++j) {
			Matrix3f rot;
			rot.setRandom();
			HouseholderQR<Matrix3f> qr(rot);
			Matrix3f q = qr.householderQ();
			// det(Q) == -1 for reflections, we convert
			// it here to a rotation in that case.
			q = q*q.determinant();
			sx.rv[j] = rotVol(q);
		}
		pop[i] = sx;
	}

	RotVolPair oldBest, bestRv(INFINITY);
	int iterNoProgress = 0;
	// Main loop
	for (int it_i = 0; it_i < numIter; ++it_i) {
		if (cb != NULL) {
			if (cb(cbData, it_i + 1) == false)
				break;
		}
		// Order according to fitness function
		sort(pop.begin(), pop.end());
		// Local optimization of current best, by running 2d rotating
		// calipers in the 3 axis
		// TODO Consider if we can skip local optimum if rough estimate
		// is "bad enough"
		RotVolPair locOpt = localOptimum(pop[0].rv[pop[0].minVolIdx()].om);
		if (locOpt.vol < bestRv.vol) {
			oldBest = bestRv;
			bestRv = locOpt;
		}

		// Stop if good enough (less than 1% improvement over 5 generations)
		if (fabsf(bestRv.vol - oldBest.vol) < 0.01f*oldBest.vol) {
			if (++iterNoProgress >= 5)
				break;
		} else {
			iterNoProgress = 0;
		}

		// Get next generation
		vector<Simplex> nextPop;
		nelderMeadBreed(pop, nextPop);
		pop = nextPop;
	}

	return bestRv;
}

// Local optimum for a given rotation matrix
struct GeneticAlgorithm::RotVolPair GeneticAlgorithm::localOptimum(const Matrix3f &rot) const
{
	// Rotate vectors, get axis box dimensions
	vector<Vector3f> rotPts;
	rotateVertex(rot, rotPts);
	Vector3f size = getAxisBoxSize(rotPts);

	// Projections to different axis. Order is important as we multiply
	// by size[i], so first we project on y,z and multiply for size in x
	// direction, etc.
	int idx[3][2] = {{1, 2}, {0, 2}, {0, 1}};
	float volume = INFINITY;
	size_t bestDim = 0;
	float projAngle = 0.f;
	for (size_t i = 0; i < sizeof idx/sizeof idx[0]; ++i) {
		vector<Eigen::Vector2f> proj;
		for (const auto &vec_it: rotPts)
			proj.push_back(Vector2f(vec_it[idx[i][0]], vec_it[idx[i][1]]));
		// Get best angle for this plane by using rotating calipers algorithm
		RotAreaPair rap = rotatingCalipers(proj);
		float dimVol = rap.area*size[i];
		if (dimVol < volume) {
			volume = dimVol;
			projAngle = -rap.angle;
			bestDim = i;
		}
	}

	Matrix3f projR;
	switch (bestDim) {
	case 0:
		projR << 1, 0, 0,
			0, cos(projAngle), sin(projAngle),
			0, -sin(projAngle), cos(projAngle);
		break;
	case 1:
		projR << cos(projAngle), 0, -sin(projAngle),
			0, 1, 0,
			sin(projAngle), 0, cos(projAngle);
		break;
	case 2:
		projR << cos(projAngle), sin(projAngle), 0,
			-sin(projAngle), cos(projAngle), 0,
			0, 0, 1;
		break;
	}

	return RotVolPair(projR*rot, volume);
}

float GeneticAlgorithm::Simplex::minVolume(void) const
{
	return rv[minVolIdx()].vol;
}

size_t GeneticAlgorithm::Simplex::minVolIdx(void) const
{
	float min = INFINITY;
	size_t minIdx = 0;
	for (size_t i = 0; i < sizeof rv/sizeof rv[0]; ++i) {
		if (rv[i].vol < min)
			minIdx = i, min = rv[i].vol;
	}

	return minIdx;
}

// Input is the current population (we will use the first M_/2
// samples), which is expected to be ordered, being the best sample
// the first one (smallest is best, as it is volume)
// pop: out, new population
void GeneticAlgorithm::nelderMeadBreed(const vector<Simplex> &currentPop,
				       vector<Simplex> &pop) const
{
	constexpr int numGroups = 4;
	vector<int> idx[numGroups];
	srand(time(NULL));
	float factor = (M_/2.f)*(float)RAND_MAX;
	int h1 = M_/2, h2 = (M_ + 1)/2;

	// Create the 4 groups (sets of random indexes)
	for (int g_i = 0; g_i < numGroups; ++g_i) {
		for (int i = 0; i < M_/2; ++i)
			idx[g_i].push_back((int) floorf(rand()/factor));
	}

	// Crossover first two
	for (int sx_i = 0; sx_i < h1; ++sx_i) {
		float fit1 = currentPop[idx[0][sx_i]].minVolume();
		float fit2 = currentPop[idx[1][sx_i]].minVolume();
		float cutoff = 0.5 + 0.1*(fit1 <= fit2) - 0.1*(fit1 >= fit2);
		Simplex sx;
		for (int om_i = 0; om_i < DIM + 1; ++om_i) {
			if (rand() < cutoff)
				sx.rv[om_i] = currentPop[idx[0][sx_i]].rv[om_i];
			else
				sx.rv[om_i] = currentPop[idx[1][sx_i]].rv[om_i];
		}
		pop.push_back(sx);
	}
	// Crossover last two
	for (int sx_i = 0; sx_i < h2; ++sx_i) {
		float fit1 = currentPop[idx[2][sx_i]].minVolume();
		float fit2 = currentPop[idx[3][sx_i]].minVolume();
		float prob = 0.5 + 0.1*(fit1 <= fit2) - 0.1*(fit1 >= fit2);
		Simplex sx;
		for (int om_i = 0; om_i < DIM + 1; ++om_i) {
			Matrix3f avg = averageMatrices(
				currentPop[idx[2][sx_i]].rv[om_i].om, prob,
				currentPop[idx[3][sx_i]].rv[om_i].om, 1 - prob);
			sx.rv[om_i] = rotVol(avg);
		}
		pop.push_back(sx);
	}

	// Mutation
	for (Simplex& sx: pop)
		sx = nelderMead(sx);
}

// Rotates the vertex by the rotation matrix rot, returning those points
// in rotPts
void GeneticAlgorithm::rotateVertex(const Matrix3f &rot, vector<Vector3f> &rotPts) const
{
	for (const auto &pt: vert_)
		rotPts.push_back(rot.transpose()*pt);
}

// Returns vertex indices of maximum and minimum coordinates in axis d
void GeneticAlgorithm::getIdxMinMax(const vector<Vector3f> &points, int d,
				    size_t &minIdx, size_t &maxIdx) const
{
	float min = numeric_limits<float>::max();
	float max = -numeric_limits<float>::max();
	size_t num_v = points.size();
	for (size_t i = 0; i < num_v; ++i) {
		if (points[i][d] < min) {
			min = points[i][d];
			minIdx = i;
		}
		if (points[i][d] > max) {
			max = points[i][d];
			maxIdx = i;
		}
	}
}

// Returns the 8 points that define the bounding box, first the four in
// one side, then the 4 in the opposite side. It returns also the
// dimensions of the bounding box.
void GeneticAlgorithm::getBBoxPoints(const Matrix3f &rot,
				     vector<Vector3f> &bbPts, vector<float> &dim) const
{
	// Rotate points
	vector<Vector3f> rotPts;
	rotateVertex(rot, rotPts);
	// Get limits of rotated points. Orden in index array:
	// 0,1: x min, x max
	// 2,3: y min, y max
	// 4,5: z min, z max
	size_t idx[6];
		for (int d = 0; d < DIM; ++d) {
		getIdxMinMax(rotPts, d, idx[2*d], idx[2*d + 1]);
		dim.push_back(rotPts[idx[2*d + 1]][d] - rotPts[idx[2*d]][d]);
	}
	vector<Vector3f> axisBox;
	// z min plane
	// (xmin, ymin), (xmax, ymin), (xmin, ymax), (xmax, ymax),
	axisBox.push_back({rotPts[idx[0]][0], rotPts[idx[2]][1], rotPts[idx[4]][2]});
	axisBox.push_back({rotPts[idx[1]][0], rotPts[idx[2]][1], rotPts[idx[4]][2]});
	axisBox.push_back({rotPts[idx[0]][0], rotPts[idx[3]][1], rotPts[idx[4]][2]});
	axisBox.push_back({rotPts[idx[1]][0], rotPts[idx[3]][1], rotPts[idx[4]][2]});
	// z max plane
	// (xmin, ymax), (xmax, ymax), (xmin, ymin), (xmax, ymin)
	axisBox.push_back({rotPts[idx[0]][0], rotPts[idx[3]][1], rotPts[idx[5]][2]});
	axisBox.push_back({rotPts[idx[1]][0], rotPts[idx[3]][1], rotPts[idx[5]][2]});
	axisBox.push_back({rotPts[idx[0]][0], rotPts[idx[2]][1], rotPts[idx[5]][2]});
	axisBox.push_back({rotPts[idx[1]][0], rotPts[idx[2]][1], rotPts[idx[5]][2]});
	// Rotate back the bounding box
	for (const auto &pt: axisBox)
		bbPts.push_back(rot*pt);
}

// Returns x-y-z dimensions of the axis oriented box enclosing a set
// of points
Vector3f GeneticAlgorithm::getAxisBoxSize(const vector<Vector3f> &points) const
{
	Vector3f size;

	for (int d = 0; d < DIM; ++d) {
		size_t minIdx, maxIdx;
		getIdxMinMax(points, d, minIdx, maxIdx);
		size[d] = points[maxIdx][d] - points[minIdx][d];
	}

	return size;
}

// Volume when bounding points inside a box with angles defined by rot
float GeneticAlgorithm::bbVolume(const Matrix3f &rot) const
{
	// Rotate points
	vector<Vector3f> rotPts;
	rotateVertex(rot, rotPts);
	// Find limit in each axis
	Vector3f size = getAxisBoxSize(rotPts);
	return size[0]*size[1]*size[2];
}

float GeneticAlgorithm::frobeniusNorm(const Matrix3f &m) const
{
	float norm = 0;
	size_t sz = m.size();
	for (size_t i = 0; i < sz; ++i) {
		float v = m(i);
		norm += v*v;
	}
	return norm;
}

// Calculate projection of matrix mat to the space of orthogonal matrices, by
// obtaining the nearest QR decomposition. We restrict the projection to
// rotations.
//
// TODO Seems to work well enough, but this is not the orthogonal projection.
// We should instead use the polar decomposition of mat [1].  However, this
// seems to actually produce the orthogonal projection when mat rotates around
// one axis, as apparently in that case the QR decomposition is the same as a
// polar decomposition. This is not true though when rotating around more than
// one axis. However, the approximation seems to be good enough for the
// Nelder-Mead method, and calculating the QR decomposition is probably faster
// than calculating the polar decomposition. This needs more research though.
//
// [1] Moakher, M. (2002). Means and averaging in the group of rotations.
// SIAM journal on matrix analysis and applications, 24(1), 1-16.
Matrix3f GeneticAlgorithm::calcProjectionToSO3(const Matrix3f &mat) const
{
	HouseholderQR<Matrix3f> qr(mat);
	Matrix3f q = qr.householderQ();
	const vector<Matrix3f> *test_mats;
	Matrix3f min_q;
	float min_dist = INFINITY;

	if (q.determinant() > 0) {
		// q is a rotation and a possible solution already, we'll
		// check also by multiplying by diagonal rotation matrices.
		test_mats = &diag_rot_mats_;
		min_q = q;
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
		}
	}

	return min_q;
}

// Mean of an array of rotation matrices.
Matrix3f GeneticAlgorithm::rotMatricesMean(const Matrix3f *rot_v, size_t rot_n) const
{
	Matrix3f mean(Matrix3f::Zero());

	for (size_t r_i = 0; r_i < rot_n; ++r_i) {
		mean += rot_v[r_i];
	}
	mean /= rot_n;

	return calcProjectionToSO3(mean);
}

Eigen::Matrix3f GeneticAlgorithm::averageMatrices(const Eigen::Matrix3f& m1, float w1,
						  const Eigen::Matrix3f& m2, float w2) const
{
	Eigen::Matrix3f avg = w1*m1 + w2*m2;
	return calcProjectionToSO3(avg);
}

GeneticAlgorithm::Simplex GeneticAlgorithm::nelderMead(const Simplex &sx) const
{
	Simplex out = sx;
	for (int i = 0; i < nelMeadIter; ++i) {
		// Order by volume
		sort(begin(out.rv), end(out.rv));

		// Centroid, excluding worst value
		Eigen::Matrix3f matv[DIM];
		for (size_t i = 0; i < DIM; ++i)
			matv[i] = out.rv[i].om;
		Eigen::Matrix3f cent = rotMatricesMean(matv, DIM);
		RotVolPair refl = rotVol(cent*out.rv[DIM].om.transpose()*cent);
		if (refl.vol < out.rv[DIM - 1].vol) {
			if (refl.vol >= out.rv[0].vol) {
				// Reflection
				out.rv[DIM] = refl;
			} else {
				// Expansion
				RotVolPair exp =
					rotVol(cent*out.rv[DIM].om.transpose()*cent*refl.om);
				if (exp.vol < refl.vol)
					out.rv[DIM] = exp;
				else
					out.rv[DIM] = refl;
			}
		} else {
			RotVolPair avg = rotVol(averageMatrices(cent, 0.5,
								out.rv[DIM].om, 0.5));
			if (avg.vol <= out.rv[DIM].vol) {
				out.rv[DIM] = avg;
			} else {
				for (size_t i = 1; i < DIM + 1; ++i)
					out.rv[i] = rotVol(averageMatrices(
								   out.rv[0].om, 0.5,
								   out.rv[i].om, 0.5));
			}
		}
	}

	return out;
}
