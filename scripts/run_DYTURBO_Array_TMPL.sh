#!/bin/bash

##
## @file run_DYTURBO_TMPL.sh
##
## Batch script template for mogon, process by submit_DYTURBO.sh
##
## @author cuto <Jakub.Cuth@cern.ch>
## @date 2015-08-06

# Setup of BSUB
#BSUB -J JOBNAME[SEEDLIST]
#BSUB -L /bin/bash
#BSUB -o OUTDIR/JOBNAME_%I.out
#BSUB -e OUTDIR/JOBNAME_%I.err
# #BSUB -q atlaslong
#BSUB -q atlasshort
#BSUB -W 5:00
#BSUB -app Reserve5G
#BSUB -n SETNPROCESSORS
#BSUB -R "rusage[atlasio=0],select[hname!='a0135']"

shopt -s expand_aliases

alias CP='rsync -avPL'

#source ATLAS + ROOT
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
localSetupROOT # 5.34.25-x86_64-slc6-gcc48-opt
# setup lhapdf
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/etapfs03/atlashpc/cuth/resbos/lhapdf/LHAPDF-6.1.4/install/lib/
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:DYTURBOROOTDIR/../RESBOS/lhapdf/lhapdf-5.6.0/install/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:DYTURBOROOTDIR/lhapdf6/lib/
#mkdir PDFsets/




# run the code
printenv
if [[ JOBNAME_${LSB_JOBINDEX} =~ ^dyturbo_ ]]
then
    if [[ -n "`find OUTDIR -name "JOBNAME_${LSB_JOBINDEX}.root" `" ]]; then exit 0; fi;
    cd /tmp/${LSB_JOBID}_${LSB_JOBINDEX}.tmpdir || echo local run: staying in `pwd`
else
    cd /jobdir/$LSB_JOBID || echo local run: staying in `pwd`
fi
rm -rf run_dir_$LSB_JOBINDEX
mkdir run_dir_$LSB_JOBINDEX
cd run_dir_$LSB_JOBINDEX
date
# cp all stuff
# bins
CP DYTURBOROOTDIR/bin/dyturbo .
CP DYTURBOROOTDIR/../DYRES/bin/dyres .
CP DYTURBOROOTDIR/../MCFM/Bin/mcfm .

# config files
CP DYTURBOROOTDIR/input/default.in .
CP DYTURBOROOTDIR/../MCFM/Bin/process.DAT .
CP DYTURBOINPUTFILE input.in
sed -i "s|= seed|= ${LSB_JOBINDEX}|g        " input.in
ln -sf input.in infile
ln -sf input.in input.DAT

# run
if [[ JOBNAME_${LSB_JOBINDEX} =~ ^dyturbo_ ]]
then
    /usr/bin/time -v ./dyturbo input.in #| tee OUTDIR/JOBNAME_${LSB_JOBINDEX}.log
    hadd -f results_merge.root results*.root
    CP results_merge.root OUTDIR/JOBNAME_${LSB_JOBINDEX}.root
elif [[ JOBNAME_${LSB_JOBINDEX} =~ ^dyres_ ]]
then
    sed -i "s|seed|1${LSB_JOBINDEX}|g        " input.in
    /usr/bin/time -v ./dyres < infile #| tee OUTDIR/JOBNAME_${LSB_JOBINDEX}.log
    #hadd -f results_merge.root results*.root
    #CP result.top OUTDIR/JOBNAME_${LSB_JOBINDEX}.top
    hadd -f results_merge.root results*.root
    CP results_merge.root OUTDIR/JOBNAME_${LSB_JOBINDEX}.root
elif [[ JOBNAME_${LSB_JOBINDEX} =~ ^mcfm_ ]]
then
    /usr/bin/time -v ./mcfm #| tee OUTDIR/JOBNAME_${LSB_JOBINDEX}.log
    CP *.C OUTDIR/JOBNAME_${LSB_JOBINDEX}.C
    hadd -f results_merge.root results*.root
    CP results_merge.root OUTDIR/JOBNAME_${LSB_JOBINDEX}.root
fi


exit 0


