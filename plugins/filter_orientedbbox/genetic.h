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

#ifndef ORIENTEDBBOX_GENETIC_H
#define ORIENTEDBBOX_GENETIC_H

#include <vector>
#include <Eigen/Dense>

/**
 * @brief The GeneticAlgorithm class
 * Uses a genetic algorithm to optimize a target function.
 */
class GeneticAlgorithm
{
	// Just for unit tests, no run-time effects
	friend struct TestGeneticAlgorithm;

public:
	// Callback type, used to return the current iteration. If it returns
	// false, the algorithm will stop.
	typedef bool CallBackIter(void *cbData, int iter);

	GeneticAlgorithm(const std::vector<Eigen::Vector3f> &vert);

	// Encodes a rotation matrix and the associated volume for our vertices
	struct RotVolPair {
		RotVolPair() = default;
		RotVolPair(float vol) : vol(vol) {}
		RotVolPair(const Eigen::Matrix3f& om, float vol) : om(om), vol(vol) {}
		// Orthogonal matrix
		Eigen::Matrix3f om;
		// Volume of minbb when using the rotation matrix
		float vol;
		bool operator<(const RotVolPair& rv) const
		{
			return vol < rv.vol;
		}
	};

        RotVolPair run(CallBackIter cb, void *cbData, int numIter = maxIter) const;
	void getBBoxPoints(const Eigen::Matrix3f &rot,
			   std::vector<Eigen::Vector3f> &bbPts,
			   std::vector<float> &dim) const;

private:
	static const int DIM = 3;
	static const int maxIter = 100;
	static const int nelMeadIter = 20;

	const std::vector<Eigen::Vector3f> &vert_;
	// Set of diagonal rotation matrices
	std::vector<Eigen::Matrix3f> diag_rot_mats_;
	// Set of diagonal reflection matrices
	std::vector<Eigen::Matrix3f> diag_reflect_mats_;

        struct RotVolPair rotVol(const Eigen::Matrix3f& om) const {
		return RotVolPair{om, bbVolume(om)};
	}
	// Simplex used in Nelder-Mead step.
	// It is composed of DIM+1 orthogonal matrices
	struct Simplex {
		RotVolPair rv[DIM + 1];
		bool operator<(const Simplex& sx) const
		{
			return minVolume() < sx.minVolume();
		}
		float minVolume(void) const;
		size_t minVolIdx(void) const;
	};
	struct RotVolPair localOptimum(const Eigen::Matrix3f &om) const;
	void rotateVertex(const Eigen::Matrix3f &rot,
			  std::vector<Eigen::Vector3f> &rotPts) const;
	void getIdxMinMax(const std::vector<Eigen::Vector3f> &points, int d,
			  size_t &minIdx, size_t &maxIdx) const;
	Eigen::Vector3f getAxisBoxSize(const std::vector<Eigen::Vector3f> &points) const;
	float bbVolume(const Eigen::Matrix3f &rot) const;
	float frobeniusNorm(const Eigen::Matrix3f &m) const;
	Eigen::Matrix3f calcProjectionToSO3(const Eigen::Matrix3f &mat) const;
	Eigen::Matrix3f rotMatricesMean(const Eigen::Matrix3f *rot_v, size_t rot_n) const;
	Eigen::Matrix3f averageMatrices(const Eigen::Matrix3f& m1, float w1,
					const Eigen::Matrix3f& m2, float w2) const;
	Simplex nelderMead(const Simplex &sx) const;
	void nelderMeadBreed(const std::vector<Simplex> &best_half,
			     std::vector<Simplex> &pop) const;
};

#endif //ORIENTEDBBOX_GENETIC_H
