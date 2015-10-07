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


merge_pt_y(){
    DRYRUN=echo 
    DRYRUN=

    defaults_seed=100101
    search="results/*_100101.root"
    search="results/dyturbo_wp_lhc7_ZPT-CT10_0_qt0100y05t*_$defaults_seed.root" # wp 
    for inrootfile in `ls $search`
    do
        root_base=`echo $inrootfile | sed "s|$defaults_seed|*|g"`
        root_merge=`echo $inrootfile | sed "s|$defaults_seed|merge|g; s|results|results_merge|g"`
        $DRYRUN hadd -f $root_merge $root_base
        echo
    done

    # fin and tot
    root_base=`echo "$search" | sed  "s|$defaults_seed|merge|g; s|results|results_merge|g"  | rev  | cut -d"_" -f2- | cut -dt -f2- | rev `t
    $DRYRUN hadd -f ${root_base}FIN_merge.root ${root_base}{CT,REAL,VIRT}_merge.root
    echo
    $DRYRUN hadd -f ${root_base}TOT_merge.root ${root_base}{RES,FIN}_merge.root
    echo
}


#merge_w_pt

merge_pt_y

exit 0
