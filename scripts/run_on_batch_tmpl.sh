#!/bin/bash

##
## @file run_on_batch_tmpl.sh
##
## Description of the script file
##
## @author Jakub Cuth <Jakub.Cuth@cern.ch>
## @date 2016-06-16



# Setup of BSUB
#BSUB -J JOBNAME[SEEDLIST]
#BSUB -L /bin/bash
#BSUB -oo OUTDIR/JOBNAME_%I.out
#BSUB -q SETQUEUE
#BSUB -W SETWALLTIME
#BSUB -app Reserve5G
#BSUB -n SETNPROCESSORS
#BSUB -R "rusage[atlasio=0] select[hname!='a0135' && hname!='a0356' && hname!='a0545']"



shopt -s expand_aliases

alias CP='rsync -avPL'

#source ATLAS + ROOT
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
lsetup root # 5.34.25-x86_64-slc6-gcc48-opt
# setup lhapdf
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/etapfs03/atlashpc/cuth/resbos/lhapdf/LHAPDF-6.1.4/install/lib/
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:DYTURBOROOTDIR/../RESBOS/lhapdf/lhapdf-5.6.0/install/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:SETLHAPDFLIB
export LHAPDF_DATA_PATH=SETLHAPDFDATA
export LHAPATH=$LHAPDF_DATA_PATH
#mkdir PDFsets/


# file exists
if [[ -n "`find OUTDIR -name "JOBNAME_${LSB_JOBINDEX}.root" `" ]]
then
    echo "Output file 'OUTDIR/JOBNAME_${LSB_JOBINDEX}.root' already exists. Bye."
    exit 0
fi;

# go to pool directory
cd /jobdir/${LSB_JOBID}-${LSB_JOBINDEX} || echo local run: staying in `pwd`
rm -rf run_dir_$LSB_JOBINDEX
mkdir run_dir_$LSB_JOBINDEX
cd run_dir_$LSB_JOBINDEX
date

# cp all stuff
CP DYTURBOROOTDIR/bin/dyturbo .
CP DYTURBOROOTDIR/input/default.in .
CP DYTURBOINPUTFILE input.in

# run
/usr/bin/time -v ./dyturbo SETPROGARGUMETS
hadd -f results_merge.root results*.root || exit 3
CP results_merge.root OUTDIR/JOBNAME_${LSB_JOBINDEX}.root

exit 0
