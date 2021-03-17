/****************************************************************************
* MeshLab                                                           o o     *
* A versatile mesh processing toolbox                             o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005-2021                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#include <common/ml_document/mesh_model.h>
#include "sampleio.h"

#include <wrap/io_trimesh/export_smf.h>
#include <wrap/io_trimesh/import_smf.h>
#include <wrap/io_trimesh/export.h>

QString SampleIOPlugin::pluginName() const
{
	return "IOMeshExample";
}

/*
	returns the list of the file's type which can be imported
*/
std::list<FileFormat> SampleIOPlugin::importFormats() const
{
	return { FileFormat("Simple Model Format", tr("SMF")) };
}

/*
	returns the list of the file's type which can be exported
*/
std::list<FileFormat> SampleIOPlugin::exportFormats() const
{
	return { FileFormat("Simple Model Format", tr("SMF")) };
}

/*
	returns the mask on the basis of the file's type.
	otherwise it returns 0 if the file format is unknown
*/
void SampleIOPlugin::exportMaskCapability(
		const QString&,
		int &capability,
		int &defaultBits) const
{
	capability=defaultBits=0;
	return;
}


bool SampleIOPlugin::open(
		const QString &,
		const QString &fileName,
		MeshModel &m,
		int& ,
		const RichParameterList & ,
		vcg::CallBackPos *,
		QWidget *)
{
	int result = vcg::tri::io::ImporterSMF<CMeshO>::Open(m.cm, qPrintable(fileName));
	if (result != vcg::tri::io::ImporterSMF<CMeshO>::E_NOERROR)
	{
		return false;
	}
	return true;
}

bool SampleIOPlugin::save(
		const QString &,
		const QString &fileName,
		MeshModel &m,
		const int mask,
		const RichParameterList &,
		vcg::CallBackPos *,
		QWidget *)
{
	QString errorMsgFormat = "Error encountered while exportering file %1:\n%2";

	int result = vcg::tri::io::ExporterSMF<CMeshO>::Save(m.cm,qPrintable(fileName),mask);
	if(result!=0) {
		errorMessage = "Saving Error: " + errorMsgFormat.arg(fileName, vcg::tri::io::Exporter<CMeshO>::ErrorMsg(result));
		return false;
	}
	return true;
}

MESHLAB_PLUGIN_NAME_EXPORTER(SampleIOPlugin)
