#!/bin/bash

#default paths wrt the script folder
SCRIPTS_PATH="$(dirname "$(realpath "$0")")"
MESHLAB_SCRIPTS_PATH=$SCRIPTS_PATH/../../meshlab/scripts/Windows
SOURCE_PATH=$SCRIPTS_PATH/../..

bash $MESHLAB_SCRIPTS_PATH/1_build.sh -s=$SOURCE_PATH $@
