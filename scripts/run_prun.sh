#!/bin/bash

##
## @file run_prun.sh
##
## Description of the script file
##
## @author Jakub Cuth <Jakub.Cuth@cern.ch>
## @date 2016-04-26


PRUN(){

    jobname=$1
    njobs=$2

    if [[ $target =~ compile ]]
    then

prun \
--bexec "chmod 777 compile_grid.sh; ./compile_grid.sh  $DYTURBOVERSION;" \
--exec "chmod 777 run_grid.sh; ./run_grid.sh ${jobname} %RNDM:800000 ;" \
--extFile dyturbo-${DYTURBOVERSION}.tar.gz \
--nJobs $njobs \
--cloud CERN \
--outDS user.${CERNUSER}.${jobname}_${gridv}/ \
--outputs results_merge.root \
--tmpDir /tmp/${CERNUSER} \
--nGBPerJob=MAX \
--rootVer=$ROOTVERSION --cmtConfig=$CMTVERSION

    else

prun \
--exec "chmod 777 run_grid.sh; ./run_grid.sh ${jobname} %RNDM:800000 ;" \
--nJobs $njobs \
--maxFileSize=20000000 \
--cloud CERN \
--outDS user.${CERNUSER}.${jobname}_${gridv}/ \
--outputs results_merge.root \
--noCompile \
--tmpDir /tmp/${CERNUSER} \
--nGBPerJob=MAX \
--rootVer=$ROOTVERSION --cmtConfig=$CMTVERSION

    fi
}


