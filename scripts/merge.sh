#!/bin/bash

##
## @file merge.sh
##
## Description of the script file
##
## @author cuto <Jakub.Cuth@cern.ch>
## @date 2015-08-13

DRYRUN=echo 
DRYRUN=


for variation in `seq 0 50` g_05 g_15 as_0117 as_0119
do
    qtregion=qt*y*
    wild=results/dyturbo_z0_lhc7_CT10nnlo_${variation}_${qtregion}_100101
    inlist=`echo $wild.root `
    for i in results/dyturbo_z0_lhc7_CT10nnlo_0_${qtregion}_100101.out
    do
        ls `echo $i | sed "s/_0_/_${variation}_/; s/.out/.root/"` > /dev/null
    done
    lobin=$ibin
    echo
    outfile=results_merge/dyturbo_z0_lhc7_CT10nnlo_${variation}_qtyMergeRESUM_100101.root
    $DRYRUN hadd -f $outfile $inlist
done

exit 0
