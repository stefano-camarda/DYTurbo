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
lsetup "root --skipConfirm"
lsetup "sft MCGenerators_lcgcmt67c/lhapdf/6.1.5 "
lsetup "root --skipConfirm"
export LHAPDF_DATA_PATH=./

mkdir -p m4 || exit 1
echo install Cuba
./install-cuba || exit 5
echo install DYTURBO
autoreconf -i || exit 2
./configure --enable-Ofast --enable-root || exit 3
make install || exit 4

./bin/dytests || exit 6

exit 0
