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

#ifndef MESHLAB_FILTER_ORIENTEDBBOX_PLUGIN_H
#define MESHLAB_FILTER_ORIENTEDBBOX_PLUGIN_H

#include <common/plugins/interfaces/filter_plugin.h>

/**
 * @brief The FilterOrientedbboxPlugin class
 * Calculates optimal bounding box for a mesh, using the HYBBRID method [1].
 *
 * [1] Chang, C. T., Gorissen, B., & Melchior, S. (2011). Fast
 * oriented bounding box optimization on the rotation group SO(3,R).
 * ACM Transactions on Graphics (TOG), 30(5), 1-16.
 */
class FilterOrientedbboxPlugin : public QObject, public FilterPlugin
{
	// keep these three lines unchanged
	Q_OBJECT
	MESHLAB_PLUGIN_IID_EXPORTER(FILTER_PLUGIN_IID)
	Q_INTERFACES(FilterPlugin)

public:
	// enum used to give an ID to every filter implemented in the plugin
	enum FileterIds { FP_CALC_ORIENTEDBBOX };

	FilterOrientedbboxPlugin();

	QString pluginName() const;
	QString vendor() const;

	QString filterName(ActionIDType filter) const;
	QString filterInfo(ActionIDType filter) const;
	FilterClass getClass(const QAction* a) const;
	FilterArity filterArity(const QAction*) const;
	int getPreConditions(const QAction *) const;
	RichParameterList initParameterList(const QAction*, const MeshModel &/*m*/);
	std::map<std::string, QVariant> applyFilter(
			const QAction* action,
			const RichParameterList & params,
			MeshDocument &md,
			unsigned int& postConditionMask,
			vcg::CallBackPos *cb);

private:
	std::map<std::string, QVariant> calcOrientedbbox(
			MeshDocument &md,
			vcg::CallBackPos *cb,
			const RichParameterList &par);

	// Method and variables used to do callbacks from algorithm
	static bool callBackIter(void *cbData, int iter);
	int numIter, currentIt;
	vcg::CallBackPos *cbPos;
};

#endif //MESHLAB_FILTER_ORIENTEDBBOX_PLUGIN_H
