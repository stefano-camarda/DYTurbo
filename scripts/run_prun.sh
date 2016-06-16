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
    shift
    njobs=$1
    shift
    arguments=$*

    if [[ $target =~ compile ]]
    then

        echo NOT TESTED $target && exit 6
 prun \
 --bexec "chmod 777 compile_grid.sh; ./compile_grid.sh  $DYTURBOVERSION;" \
 --exec "chmod 777 run_grid.sh; ./run_grid.sh ${arguments} --seed %RNDM ;" \
 --extFile dyturbo-${DYTURBOVERSION}.tar.gz \
 --nJobs $njobs \
 --outDS ${CERNGROUP}.${CERNUSER}.${jobname}_${gridv}/ \
 --outputs results_merge.root \
 --tmpDir /tmp/${CERNUSER} \
 --nGBPerJob=MAX \
 --rootVer=$ROOTVERSION --cmtConfig=$CMTVERSION 

    else

prun \
--exec "chmod 777 run_grid.sh; ./run_grid.sh ${arguments} --seed %RNDM ;" \
--nJobs $njobs \
--maxFileSize=25000000 \
--outDS ${CERNGROUP}.${CERNUSER}.${jobname}_${gridv}/ \
--outputs results_merge.root \
--noCompile \
--tmpDir /tmp/${CERNUSER} \
--nGBPerJob=MAX \
--rootVer=$ROOTVERSION --cmtConfig=$CMTVERSION $OFFICIAL

    fi
}


