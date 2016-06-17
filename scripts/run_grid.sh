#!/bin/bash

##
## @file run_grid.sh
##
## Description of the script file
##
## @author cuto <Jakub.Cuth@cern.ch>
## @date 2015-10-27

arguments=$*

# setup ENV
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
# ROOT
lsetup "root 6.04.14-x86_64-slc6-gcc49-opt"

# LHAPDF
if [[ SETTARGET =~ compile ]]
then
    # setup sft
    lsetup "sft releases/MCGenerators/lhapdf/6.1.5-2f446"
    lhapdf-config --version || exit 2
    # official sets
    export LHAPDF_DATA_PATH=/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current/
else
    # own lhapdf set
    export LHAPDF_DATA_PATH=./LHAPDF
fi
export LHAPATH=$LHAPDF_DATA_PATH
# dyturbo libs
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:./lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/cvmfs/sft.cern.ch/lcg/external/gcc/4.8.1/x86_64-slc6-gcc48-opt/lib/../lib64

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

./bin/dyturbo $arguments

hadd -f results_merge.root results*.root


