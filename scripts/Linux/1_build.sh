#!/bin/bash

#default paths wrt the script folder
SCRIPTS_PATH="$(dirname "$(realpath "$0")")"
MESHLAB_SCRIPTS_PATH=$SCRIPTS_PATH/../../meshlab/scripts/Linux
SOURCE_PATH=$SCRIPTS_PATH/../..
INSTALL_PATH=$SOURCE_PATH/install/

bash $MESHLAB_SCRIPTS_PATH/1_build.sh -s=$SOURCE_PATH -i=$INSTALL_PATH $@ 

mkdir -p $INSTALL_PATH/meshlab/plugins/
mv $INSTALL_PATH/usr/meshlab/plugins/* $INSTALL_PATH/meshlab/plugins/
