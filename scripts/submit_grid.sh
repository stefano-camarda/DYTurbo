#!/bin/bash

##
## @file submit_grid.sh
##
## Description of the script file
##
## @author cuto <Jakub.Cuth@cern.ch>
## @date 2015-10-27



DRYRUN=echo 
DRYRUN=

input=dyturbo_test_t`date +%s`


# lsetup rucio panda
# voms-proxy-init atlas

# make tarbal
tarbalfile=scripts/infiles/${input}.tar.gz
tar -cvzf $tarbalfile --transform='s|.*/||g' bin/dyturbo grid/infile.in scripts/run_grid.sh


$DRYRUN prun --exec ". run_grid.sh ${input} " \
     --outDS user.jcuth.${input}.out \
     --outputs=HIST:results_merge.root \
     --nJobs 1 \
     --rootVer=6.02/12 --cmtConfig=x86_64-slc6-gcc48-opt \
     --inTarBall=$tarbalfile \
     #--site ANALY_CERN_SHORT
     #--excludeFile="out_*"


