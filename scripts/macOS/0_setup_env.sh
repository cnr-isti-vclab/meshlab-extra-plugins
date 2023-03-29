#!/bin/bash
# this is a script shell sets up an ubuntu (18.04, 20.04 and 22.04) environment where
# MeshLab can be compiled.
#
# Run this script if you never installed any of the MeshLab dependencies.

SCRIPTS_PATH=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
MESHLAB_SCRIPTS_PATH=$SCRIPTS_PATH/../../meshlab/scripts/macOS

bash $MESHLAB_SCRIPTS_PATH/0_setup_env.sh $@

