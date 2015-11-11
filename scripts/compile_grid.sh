#!/bin/bash

##
## @file compile_grid.sh
##
## Description of the script file
##
## @author cuto <Jakub.Cuth@cern.ch>
## @date 2015-11-11


# setup ENV
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
# LHAPDF and ROOT
lsetup "sft --cmtConfig=x86_64-slc6-gcc47-opt MCGenerators_lcgcmt67c/lhapdf/6.1.5 "  "root --skipConfirm"

mkdir m4 || exit 1
echo install Cuba
./install-cuba || exit 5
echo install DYTURBO
autoreconf -i || exit 2
./configure --enable-Ofast --enable-root || exit 3
make install || exit 4

exit 0
