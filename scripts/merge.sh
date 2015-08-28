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


mergre_z_run(){
    for variation in `seq 0 50` g_05 g_15 as_0117 as_0119
    do
        qtregion=qt*y*
        wild=results/dyturbo_z0_lhc7_CT10nnlo_${variation}_${qtregion}_100101
        inlist=`echo $wild.root `
        for i in results/dyturbo_z0_lhc7_CT10nnlo_0_${qtregion}_100101.out
        do
            ls `echo $i | sed "s/_0_/_${variation}_/; s/.out/.root/"` > /dev/null
        done
        echo
        outfile=results_merge/dyturbo_z0_lhc7_CT10nnlo_${variation}_qtyMergeCT_100101.root
        $DRYRUN hadd -f $outfile $inlist
    done
}

merge_w_pt(){
    #resumct
    outfile=results_merge/dyturbo_wp_lhc7_CTZPT2_0_qtyMergetRESCT_100101.root
    inlist=`echo results/dyturbo_wp_lhc7_CTZPT2_0_qt*y*tRESCT_100101.root`
    #$DRYRUN hadd -f $outfile $inlist
    #real
    outfile=results_merge/dyturbo_wp_lhc7_CTZPT2_0_qtyMergetREAL_100101.root
    inlist=`echo results/dyturbo_wp_lhc7_CTZPT2_0_qt0020y05tREAL_*.root`
    $DRYRUN hadd -f $outfile $inlist
    #virt
    outfile=results_merge/dyturbo_wp_lhc7_CTZPT2_0_qtyMergetVIRT_100101.root
    inlist=`echo results/dyturbo_wp_lhc7_CTZPT2_0_qt0020y05tVIRT_*.root`
    $DRYRUN hadd -f $outfile $inlist
}


merge_w_pt

exit 0
