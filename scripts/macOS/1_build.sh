#!/bin/bash

#default paths wrt the script folder
SCRIPTS_PATH=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
MESHLAB_SCRIPTS_PATH=$SCRIPTS_PATH/../../meshlab/scripts/macOS
SOURCE_PATH=$SCRIPTS_PATH/../..
INSTALL_PATH=$SOURCE_PATH/install/

bash $MESHLAB_SCRIPTS_PATH/1_build.sh -s=$SOURCE_PATH -i=$INSTALL_PATH $@
