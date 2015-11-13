#!/bin/bash

##
## @file run_grid.sh
##
## Description of the script file
##
## @author cuto <Jakub.Cuth@cern.ch>
## @date 2015-10-27

# setup ENV
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
# LHAPDF and ROOT
lsetup "root --skipConfirm"
#lsetup "sft MCGenerators_lcgcmt67c/lhapdf/6.1.5 " 
# own lhapdf set
#export LHAPDF_DATA_PATH=/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current/
export LHAPDF_DATA_PATH=./
export LHAPATH=$LHAPDF_DATA_PATH
# dyturbo libs
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:./lib

echo
echo LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH
echo

# compile just for sure
#make install

jobname=$1
ln -sf scripts/infiles/$jobname.in input.in
random_seed=$2
sed -i "s|= seed|= 1$random_seed|g        " input.in

echo
echo input.in :
cat input.in
echo
rm -f results*.root

/usr/bin/time -v ./bin/dyturbo input.in
hadd -f results_merge.root results*.root
