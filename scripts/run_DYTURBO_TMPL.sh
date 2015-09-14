#!/bin/bash

##
## @file run_DYTURBO_TMPL.sh
##
## Batch script template for mogon, process by submit_DYTURBO.sh
##
## @author cuto <Jakub.Cuth@cern.ch>
## @date 2015-08-06

# Setup of BSUB
#BSUB -J JOBNAME
#BSUB -L /bin/bash
#BSUB -o OUTDIR/JOBNAME.out
#BSUB -e OUTDIR/JOBNAME.err
# #BSUB -q atlaslong
#BSUB -q atlasshort
#BSUB -W 5:00
#BSUB -app Reserve5G
#BSUB -n SETNPROCESSORS
#BSUB -R 'rusage[atlasio=0]'

shopt -s expand_aliases

alias CP='rsync -avPL'

#source ATLAS + ROOT
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
localSetupROOT # 5.34.25-x86_64-slc6-gcc48-opt
# setup lhapdf
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/etapfs03/atlashpc/cuth/resbos/lhapdf/LHAPDF-6.1.4/install/lib/



# run the code
printenv
cd /jobdir/$LSB_JOBID || echo local run: staying in `pwd`
rm -rf run_dir
mkdir run_dir
cd run_dir
date
# cp all resbos stuff
CP DYTURBOROOTDIR/bin/dyturbo .
CP DYTURBOROOTDIR/input/default.in .
CP DYTURBOINPUTFILE input.in
/usr/bin/time -v ./dyturbo input.in | tee OUTDIR/JOBNAME.log
hadd -f results_merge.root results*.root
CP results_merge.root OUTDIR/JOBNAME.root


exit 0


