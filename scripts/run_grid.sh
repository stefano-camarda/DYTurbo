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

# random seed (or PDF variations)
if [[ $jobname =~ _array_ ]]
then
    sed -i "s|= seed|= 123|g           " input.in
    sed -i "s|= array|= $random_seed|g " input.in
else
    sed -i "s|= seed|= 1$random_seed|g " input.in
fi

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

/usr/bin/time -v ./bin/dyturbo input.in
hadd -f results_merge.root results*.root


