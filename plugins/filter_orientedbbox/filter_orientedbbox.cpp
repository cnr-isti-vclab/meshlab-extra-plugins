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

#include "filter_orientedbbox.h"
#include "genetic.h"

#include <vcg/complex/algorithms/convex_hull.h>

using namespace vcg;

/**
 * @brief
 * Constructor usually performs only two simple tasks of filling the two lists
 *  - typeList: with all the possible id of the filtering actions
 *  - actionList with the corresponding actions.
 * If you want to add icons to your filtering actions you can do here by
 * construction the QActions accordingly
 */
FilterOrientedbboxPlugin::FilterOrientedbboxPlugin()
{
	typeList = { FP_CALC_ORIENTEDBBOX };

	for(const ActionIDType& tt : typeList)
		actionList.push_back(new QAction(filterName(tt), this));
}

QString FilterOrientedbboxPlugin::pluginName() const
{
	return "FilterOrientedbbox";
}

QString FilterOrientedbboxPlugin::vendor() const
{
	return "Alfonso Sanchez-Beato";
}

/**
 * @brief filterName() must return the very short string describing each
 * filtering action (this string is used also to define the menu entry)
 * @param filterId: the id of the filter
 * @return the name of the filter
 */
QString FilterOrientedbboxPlugin::filterName(ActionIDType filterId) const
{
	switch(filterId) {
	case FP_CALC_ORIENTEDBBOX :
		return "Oriented bounding box";
	default :
		assert(0);
		return "";
	}
}


/**
 * @brief filterInfo() must return the longer string describing each filtering action
 * (this string is used in the About plugin dialog)
 * @param filterId: the id of the filter
 * @return an info string of the filter
 */
QString FilterOrientedbboxPlugin::filterInfo(ActionIDType filterId) const
{
	switch(filterId) {
	case FP_CALC_ORIENTEDBBOX :
        return "Calculate oriented bounding box of a mesh using the HYBBRID method. The "
               "algorithm this filter implements is based on the paper: <br>"
               "<i>Chang, C. T., Gorissen, B., & Melchior, S.</i><br>"
               "<b>Fast oriented bounding box optimization on the rotation group SO(3,R)</b><br>"
               "ACM Transactions on Graphics (TOG), 30(5), 1-16.";
	default :
		assert(0);
		return "Unknown Filter";
	}
}

/**
 * @brief The FilterClass describes in which generic class of filters it fits.
 * This choice affects the submenu in which each filter will be placed
 * More than a single class can be chosen.
 * @param a: the action of the filter
 * @return the class of the filter
 */
FilterOrientedbboxPlugin::FilterClass
FilterOrientedbboxPlugin::getClass(const QAction *a) const
{
	switch(ID(a)) {
	case FP_CALC_ORIENTEDBBOX :
		return FilterPlugin::Remeshing;
	default :
		assert(0);
		return FilterPlugin::Generic;
	}
}

/**
 * @brief FilterSamplePlugin::filterArity
 * @return
 */
FilterPlugin::FilterArity FilterOrientedbboxPlugin::filterArity(const QAction*) const
{
	return SINGLE_MESH;
}

/**
 * @brief FilterSamplePlugin::getPreConditions
 * @return
 */
int FilterOrientedbboxPlugin::getPreConditions(const QAction*) const
{
	return MeshModel::MM_NONE;
}

/**
 * @brief This function returns a list of parameters needed by each filter.
 * For each parameter you need to define,
 * - the name of the parameter,
 * - the default value
 * - the string shown in the dialog
 * - a possibly long string describing the meaning of that parameter (shown
 *   as a popup help in the dialog)
 * @param action
 * @param m
 */
RichParameterList FilterOrientedbboxPlugin::initParameterList(const QAction *action,
							      const MeshModel &/*m*/)
{
	RichParameterList parlst;
	switch(ID(action)) {
	case FP_CALC_ORIENTEDBBOX:
		parlst.addParam(RichInt("MaxIter",
					100,
					"Max iterations",
					"Maximum number of iterations. Note that the "
					"algorithm will finish before reaching these "
					" iterations if no further progress is "
					"happening. If the solution is not good enough"
					" you could experiment with changing this"
					" value.\n"));
		break;
	default:
		assert(0);
	}
	return parlst;
}

/**
 * @brief The Real Core Function doing the actual mesh processing.
 * @param action
 * @param md: an object containing all the meshes and rasters of MeshLab
 * @param par: the set of parameters of each filter, with the values set by the user
 * @param cb: callback object to tell MeshLab the percentage of execution of the filter
 * @return true if the filter has been applied correctly, false otherwise
 */
std::map<std::string, QVariant> FilterOrientedbboxPlugin::applyFilter(
		const QAction * action,
		const RichParameterList &par,
		MeshDocument &md,
		unsigned int &/*postConditionMask*/,
		vcg::CallBackPos *cb)
{
	switch(ID(action)) {
	case FP_CALC_ORIENTEDBBOX:
		calcOrientedbbox(md, cb, par);
		break;
	default:
		wrongActionCalled(action);
	}
	return std::map<std::string, QVariant>();
}

bool FilterOrientedbboxPlugin::callBackIter(void *cbData, int iter)
{
	FilterOrientedbboxPlugin *filt = static_cast<FilterOrientedbboxPlugin *>(cbData);
	filt->currentIt = iter;
	return filt->cbPos(100*filt->currentIt/filt->numIter, "Genetic algorithm loop");
}

void FilterOrientedbboxPlugin::calcOrientedbbox(
		MeshDocument &md,
		vcg::CallBackPos *cb,
		const RichParameterList &par)
{
	CMeshO &m = md.mm()->cm;

	log("Calculating bounding box for %d vertices", m.vn);

	// Use convex hull to minimize number of points
	MeshModel &mm = *md.mm();
	MeshModel &convexm = *md.addNewMesh("", "Convex Hull");
	convexm.updateDataMask(MeshModel::MM_FACEFACETOPO);
	bool result = vcg::tri::ConvexHull<CMeshO, CMeshO>::ComputeConvexHull(
		mm.cm, convexm.cm);
	convexm.clearDataMask(MeshModel::MM_FACEFACETOPO);
	convexm.updateBoxAndNormals();
	if (!result)
		throw MLException("Failed computing convex hull.");

	// Convert CMeshO to Eigen::Vector3f so we can use the vertex in the algorithm
	std::vector<Eigen::Vector3f> vert;
	auto end_it = convexm.cm.vert.end();
	for (auto it = convexm.cm.vert.begin(); it != end_it; ++it) {
		Eigen::Vector3f ev;
		it->cP().ToEigenVector(ev);
		vert.push_back(ev);
	}

	// Calculate now rotation matrix and bounding box mesh
	numIter = par.getInt("MaxIter");
	currentIt = 0;
	cbPos = cb;
	GeneticAlgorithm ga(vert);
	GeneticAlgorithm::RotVolPair rvp = ga.run(callBackIter, this, numIter);
	std::vector<Eigen::Vector3f> bbox;
	std::vector<float> dim;
	ga.getBBoxPoints(rvp.om, bbox, dim);
	log("Finished after %d iterations", currentIt);

	log("Rotation matrix:");
	log("%f\t%f\t%f", rvp.om(0), rvp.om(1), rvp.om(2));
	log("%f\t%f\t%f", rvp.om(3), rvp.om(4), rvp.om(5));
	log("%f\t%f\t%f", rvp.om(6), rvp.om(7), rvp.om(8));
	log("Box is %f x %f x %f", dim[0], dim[1], dim[2]);
	log("Volume is %f", rvp.vol);

	MeshModel &bboxm = *md.addNewMesh("", "Bounding box");
	size_t num_vert_box = bbox.size();
	static constexpr size_t num_cube_faces = 12;
	tri::Allocator<CMeshO>::AddVertices(bboxm.cm, num_vert_box);
	tri::Allocator<CMeshO>::AddFaces(bboxm.cm, num_cube_faces);
	for (size_t i = 0; i < num_vert_box; ++i) {
		for (size_t c = 0; c <3 ; ++c)
			bboxm.cm.vert[i].P()[c] = bbox[i][c];
	}
	// 12 faces of the cube (must be triangles)
	size_t vertIdx[num_cube_faces][3] = {
		{0,2,1}, {1,2,3}, {2,4,3}, {3,4,5}, {4,6,5}, {5,6,7},
		{4,2,6}, {6,2,0}, {6,0,7}, {7,0,1}, {7,1,5}, {5,1,3}};
	for (size_t fc_i = 0; fc_i < num_cube_faces; ++fc_i) {
		for (size_t i = 0; i < 3; ++i)
			bboxm.cm.face[fc_i].V(i) = &bboxm.cm.vert[vertIdx[fc_i][i]];
	}

	bboxm.updateBoxAndNormals();
}

MESHLAB_PLUGIN_NAME_EXPORTER(FilterSamplePlugin)
