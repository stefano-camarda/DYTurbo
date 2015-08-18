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
    inlist=
    lobin=unset
    hibin=unset
    for ibin in 0 2 4 6 8 10 12 14 16 18 22 26 30 34 38 42 46 50 54 60 70 80 100 150 200 300 800
    do
        if [[ $lobin == unset ]]
        then
            lobin=$ibin
            continue
        fi
        hibin=$ibin
        qtregion=qt${lobin}${hibin}
        inlist="$inlist results/dyturbo_z0_lhc7_CT10nnlo_${variation}_${qtregion}_100101.root"
        lobin=$ibin
    done
    echo
    outfile=results_merge/dyturbo_z0_lhc7_CT10nnlo_${variation}_qtMerge_100101.root
    $DRYRUN hadd -f $outfile $inlist
done

exit 0
