#!/bin/bash

##
## @file run_grid.sh
##
## Description of the script file
##
## @author cuto <Jakub.Cuth@cern.ch>
## @date 2015-10-27

jobname=$1
ln -sf infiles/${jobname}.in input.in
random_seed=$2
sed -i "s|= seed|= 1$random_seed|g        " input.in


# setup ENV
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
# LHAPDF and ROOT
lsetup "root 6.04.14-x86_64-slc6-gcc49-opt"
export PATH=/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/lhapdf/6.1.5-6570e/x86_64-slc6-gcc49-opt/bin:$PATH
lhapdf-config --version || exit 2
# own lhapdf set
#export LHAPDF_DATA_PATH=./
# official sets
export LHAPDF_DATA_PATH=/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current/
export LHAPATH=$LHAPDF_DATA_PATH
# dyturbo libs
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:./lib

echo
echo LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH
echo

# compile just for sure
#make install

echo
echo input.in :
cat input.in
echo
rm -f results*.root

/usr/bin/time -v ./bin/dyturbo input.in
hadd -f results_merge.root results*.root


