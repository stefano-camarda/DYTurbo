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
    proc=z0
    #resumct
    echo
    outfile=results_merge/dyturbo_${proc}_lhc7_CTZPT2_0_qtyMergetRES3DCT_100101.root
    #inlist=`echo results/dyturbo_${proc}_lhc7_CTZPT2_0_qt*y05tRESCT_100101.root`
    inlist=`echo results/dyturbo_${proc}_lhc7_CTZPT2_0_qt*y*tRES3DCT_100101.root`
    $DRYRUN hadd -f $outfile $inlist
    exit 0
    #real
    echo
    outfile=results_merge/dyturbo_${proc}_lhc7_CTZPT2_0_qtyMergetREAL_100101.root
    inlist=`echo results/dyturbo_${proc}_lhc7_CTZPT2_0_qt0020y05tREAL_*.root`
    inlist=`echo results/dyturbo_${proc}_lhc7_CTZPT2_0_qt0020y05tREAL_100101.root`
    $DRYRUN hadd -f $outfile $inlist
    #virt
    echo
    outfile=results_merge/dyturbo_${proc}_lhc7_CTZPT2_0_qtyMergetVIRT_100101.root
    inlist=`echo results/dyturbo_${proc}_lhc7_CTZPT2_0_qt0020y05tVIRT_*.root`
    inlist=`echo results/dyturbo_${proc}_lhc7_CTZPT2_0_qt0020y05tVIRT_100101.root`
    $DRYRUN hadd -f $outfile $inlist
}


merge_w_pt

exit 0
