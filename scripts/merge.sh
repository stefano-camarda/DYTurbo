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

merge_wwidth(){
    DRYRUN=echo 
    DRYRUN=
    defaultsTERM=RES
    #outdir=results_merge/wwidth_D0_zmumu
    search="results/dyturbo_*_*_CT10nnlo_0_f*qt01000y-55t${defaultsTERM}_100101.root"
    search="results_merge/wwidth_D0_zmumu/dyturbo_*_*_CT10nnlo_0_f*qt01000y-55t${defaultsTERM}_seed_10100.root"
    for inrootfile in `ls $search`
    do
        root_base=`echo $inrootfile | sed "s|${defaultsTERM}_|*|g"`
        #root_merge=`echo $inrootfile | sed "s|${defaultsTERM}_|TOT_|g; s|results|$outdir|g"`
        $DRYRUN hadd -f $root_merge $root_base
        echo
    done
}


merge_pt_y(){
    DRYRUN=echo 
    DRYRUN=

    outdir=results_merge/wm_real_1000seeds_151029
    mkdir -p $outdir

    defaults_seed=seed_1000
    search="results/*_100101.root"
    search="results/dyturbo_wm_lhc7_ZPT-CT10_0_qt0100y05t*_$defaults_seed.root" # DYTURBO wpm,z0
    #search="results/dyturbo_z0_lhc7_ZPT-CT10_0_qt0100y05t*_$defaults_seed.root" # z0 zpt
    #search="results/dyturbo_z0_lhc7_WZZPT-CT10_0_qt0100y05t*_$defaults_seed.root"  # z0 wz zpwt
    #search="results/dyres_w{p,m}_lhc7_ZPT-CT10_0_qt0100y05t*_seed_$defaults_seed.root" # DYRES seed
    #search="results/dyres_wp_lhc7_ZPT-CT10_0_qt0100y05t*_$defaults_seed.root"
    for inrootfile in `ls $search`
    do
        root_base=`echo $inrootfile | sed "s|$defaults_seed|*|g"`
        first_file=`echo $inrootfile | sed "s|$defaults_seed|seed_1015|g"`
        root_merge=`echo $inrootfile | sed "s|$defaults_seed|merge|g; s|results|$outdir|g"`
        $DRYRUN hadd -f $root_merge $root_base
        echo
        #if [[ $root_merge =~ REAL ]]
        #then # outlier removal
        root_merge=`echo $inrootfile | sed "s|$defaults_seed|outliers|g; s|results|$outdir|g"`
        $DRYRUN ./../merger/bin/merger $root_merge $first_file  $root_base
        #fi
        echo
    done

    # no hadd terms
    return


    # fin and tot
    root_base=`echo "$search" | sed  "s|$defaults_seed|merge|g; s|results|results_merge|g"  | rev  | cut -d"_" -f2- | cut -dt -f2- | rev `t
    $DRYRUN hadd -f ${root_base}FIN_merge.root ${root_base}{CT,REAL,VIRT}_merge.root
    echo
    $DRYRUN hadd -f ${root_base}TOT_merge.root ${root_base}{RES,FIN}_merge.root
    echo
    #
    root_base=`echo "$search" | sed  "s|$defaults_seed|outliers|g; s|results|results_merge|g"  | rev  | cut -d"_" -f2- | cut -dt -f2- | rev `t
    $DRYRUN hadd -f ${root_base}FIN_outliers.root ${root_base}{CT,REAL,VIRT}_outliers.root
    echo
    $DRYRUN hadd -f ${root_base}TOT_outliers.root ${root_base}{RES,FIN}_outliers.root
    echo
}

mergedir=results_merge/grid_151201
MERGER="./bin/merger -X"
merge_grid(){
    DRYRUN=echo 
    DRYRUN=
    #griddownloads='results_grid/user.jcuth.*_0_* results_grid/group.phys-sm.*_5{1,2,3,4}_*'
    griddownloads=`echo results_grid/user.jcuth.*_0_* `
    griddownloads="$griddownloads `echo results_grid/group.phys-sm.dyturbo_{wp,wm,z0}_*_4_*_seed.out_HIST`"
    #griddownloads="$griddownloads `echo results_grid/group.phys-sm.*_wm_*_5{1,2}_*`"
    #griddownloads="$griddownloads `echo results_grid/group.phys-sm.*_wp_*_5{1,2}_*`"
    #griddownloads="$griddownloads `echo results_grid/group.phys-sm.*_z0_*_5{1,2}_*`"
    #griddownloads="$griddownloads `echo results_grid/group.phys-sm.*_wm_*_5{3,4}_*`"
    #griddownloads="$griddownloads `echo results_grid/group.phys-sm.*_wp_*_5{3,4}_*`"
    #griddownloads="$griddownloads `echo results_grid/group.phys-sm.*_z0_*_5{3,4}_*`"
    mkdir -p $mergedir
    for griddir in `echo $griddownloads`
    do
        echo $griddir
        base=`echo $griddir | cut -d. -f3 | sed "s|_seed|_merge|g"`
        outfile=$mergedir/$base.root
        echo skip $outfile
        #$DRYRUN hadd -f $mergedir/$base.root $griddir/*results_merge.root*
        base=`echo $base | sed "s|_merge|_outliers|g"`
        outfile=$mergedir/$base.root
        if [ -f $outfile  ]
        then
            echo skip $outfile
        else
            $DRYRUN $MERGER $outfile $griddir/*results_merge.root*
        fi
    done
}

merge_grid_TOT(){
    DRYRUN=echo 
    DRYRUN=
    for fres in `ls $mergedir/*RES*outliers.root`
    do
        ffin=`echo  $fres | sed "s|RES_|FIN_|g"`
        ftot=`echo  $fres | sed "s|RES_|TOT_|g"`
        fct=`echo   $fres | sed "s|RES_|CT_|g"`
        freal=`echo $fres | sed "s|RES_|REAL_|g"`
        fvirt=`echo $fres | sed "s|RES_|VIRT_|g"`
        $DRYRUN hadd -f $ffin $fct $freal $fvirt
        $DRYRUN hadd -f $ftot $fres $ffin
    done
}


MERGER="hadd -f "
merge_cubatures(){
    DRYRUN=echo 
    DRYRUN=
    tmp_phase=qt010y01
    tot_phase=qt0100ym55
    #outdir=results_merge/quad_160118
    #resdir=/etapfs03/atlashpc/cuth/DYTURBO_PROD/dyturbo-0.9.6.2
    #resdir=/etapfs03/atlashpc/cuth/DYTURBO_PROD/results_wrongMURMUF_for_WpWm
    resdir=/etapfs03/atlashpc/cuth/DYTURBO_PROD/dyturbo-0.9.6.2/results
    prodname=benchmark_v2_160204_WZ
    outdir=results_merge/$prodname
    for proc in wp wm # wp wm z0
    do
        for pdfmem in 0 # `seq 0 54`
        do
            mem=$(( $pdfmem + 100 ))
            #name="results/dyturbo_wm_lhc7_WZZPT-CT10_array_${tmp_phase}t*3D_seed_$mem.root"
            #name="results/dyturbo_wm_lhc7_WZZPT-CT10_array_${tmp_phase}tRES3D_seed_$mem.root"
            name=$resdir/"dyturbo_${proc}_lhc7_CT10nnlo_0_bm2${tmp_phase}t*3D_seed_10100.root"
            mkdir -p $outdir
            for f in `ls $name`
            do
                infiles=`echo  $f | sed "s|$tmp_phase|*|g" `
                #outfile=`echo  $f | sed "s|$tmp_phase|$tot_phase|g; s|seed_[0-9]*|seed_merge|g; s|results|$outdir|g; s|array|$pdfmem|g;" `
                outfile=$outdir/`basename  $f | sed "s|$tmp_phase|$tot_phase|g; s|seed_[0-9]*|seed_outliers|g;" `
                $DRYRUN $MERGER $outfile $infiles
            done
        done
    done
}

merge_stefano(){
    DRYRUN=echo 
    DRYRUN=
    #
    indir=results_Stefano/
    #outdir=results_merge/Stefano_dyturbo_v1_160201_o2
    outdir=results_merge/Stefano_dyturbo_v1_160204_WZ
    mkdir -p $outdir

    for proc in Wp Wm Z
    do
        for term in v r
        do
            # adding '0' for 'z0' and then cut on first two characters for 'wp' and 'wm'
            procOut=`echo "${proc}0" | tr '[:upper:]' '[:lower:]'`
            outf=$outdir/${procOut:0:2}${term}.root
            inf=$indir/${proc}${term}nnlo/benchmark/CT10nnlo/0/integr/*/AiMoments.root
            $DRYRUN ./bin/merger $outf $inf 
            #exit 0
        done
    done
}

merge_benchmark(){
    DRYRUN=echo 
    DRYRUN=
    #
    resdir=/etapfs03/atlashpc/cuth/DYTURBO_PROD
    #
    #resdir=$resdir/dyturbo-0.9.6/results_benchmark0
    #prodname=benchmark_v0_160125
    #resdir=$resdir/dyturbo-0.9.6/results_benchmark1
    #prodname=benchmark_v1_160125
    #resdir=$resdir/dyturbo-0.9.6.1/results
    #prodname=benchmark_v0.1_160125
    #testsample="*{REAL,VIRT}"
    #seednum=10105
    #
    resdir=$resdir/results_wrongMURMUF_for_WpWm
    #testsample="w*bm0*{tCT}"
    #resdir=$resdir/dyturbo-0.9.6.2/results/
    #testsample="w*bm0*{RES,CT}"
    seednum=11105
    #
    #prodname=benchmark_v1_160202
    outdir=results_merge/$prodname
    mkdir -p $outdir
    #
    for fres in `ls $resdir/dyturbo_z0*bm1*{RES,CT}*$seednum*root`
    do
        #echo $fres
        infiles=`echo $fres | sed "s|$seednum|*|g"`
        outfilebase=`basename $fres | sed "s|$seednum.*||g;s|y-55|ym55|g;s|bm1|bm0|g"`
        #$DRYRUN hadd -f $outdir/${outfilebase}merge.root $infiles
        $DRYRUN ./bin/merger -X $outdir/${outfilebase}outliers.root $infiles #&& return 0 || return 3
        echo
    done
}

merge_general(){
    DRYRUN=echo 
    DRYRUN=
    #
    resdir=/home/cuth/work/unimainz/ATLAS/workarea/DYTURBO/
    qtymerge=qt010y01
    #
    rm -f mergeconf
    # z0
    #echo  " z0 bm1          {REAL,VIRT}  $resdir/dyturbo-0.9.6/results_benchmark1 bm0 10105 benchmark_v0_160204_WZ " >> mergeconf
    #echo  " z0 bm1          {REAL,VIRT}  $resdir/dyturbo-0.9.6/results_benchmark1 bm1 10105 benchmark_v1_160204_WZ " >> mergeconf
    #echo  " z0 bm1          {REAL,VIRT}  $resdir/dyturbo-0.9.6/results_benchmark1 bm2 10105 benchmark_v2_160204_WZ " >> mergeconf
    #echo " z0 bm0          {REAL,VIRT}  $resdir/dyturbo-0.9.6.2/results_bm0_WZREALVIRT          bm0 11105 benchmark_v0_160204_WZ " >> mergeconf
    #echo " z0 bm0          {REAL,VIRT}  $resdir/dyturbo-0.9.6.2/results_bm0_WZREALVIRT          bm1 11105 benchmark_v1_160204_WZ " >> mergeconf
    #echo " z0 bm0          {REAL,VIRT}  $resdir/dyturbo-0.9.6.2/results_bm0_WZREALVIRT          bm2 11105 benchmark_v2_160204_WZ " >> mergeconf
    #echo  " z0 bm0          {RES,CT}     $resdir/results_wrongMURMUF_for_WpWm     bm0 11105 benchmark_v0_160204_WZ " >> mergeconf
    #echo  " z0 bm1          {RES,CT}     $resdir/results_wrongMURMUF_for_WpWm     bm1 11105 benchmark_v1_160204_WZ " >> mergeconf
    #echo  " z0 bm2$qtymerge {RES3D,CT3D} $resdir/results_wrongMURMUF_for_WpWm     bm2 10100 benchmark_v2_160204_WZ " >> mergeconf

    # wpm
    #echo  " w{p,m} bm0          VIRT         $resdir/dyturbo-0.9.6.2/results_bm012_WpWm     bm0 11105 benchmark_v0_160204_WZ " >> mergeconf
    #echo  " w{p,m} bm0          VIRT         $resdir/dyturbo-0.9.6.2/results_bm012_WpWm     bm1 11105 benchmark_v1_160204_WZ " >> mergeconf
    #echo  " w{p,m} bm0          VIRT         $resdir/dyturbo-0.9.6.2/results_bm012_WpWm     bm2 11105 benchmark_v2_160204_WZ " >> mergeconf
    #echo " w{p,m} bm0          REAL         $resdir/dyturbo-0.9.6.2/results_bm0_WZREALVIRT bm0 11105 benchmark_v0_160204_WZ " >> mergeconf
    #echo " w{p,m} bm0          REAL         $resdir/dyturbo-0.9.6.2/results_bm0_WZREALVIRT bm1 11105 benchmark_v1_160204_WZ " >> mergeconf
    #echo " w{p,m} bm0          REAL         $resdir/dyturbo-0.9.6.2/results_bm0_WZREALVIRT bm2 11105 benchmark_v2_160204_WZ " >> mergeconf
    #echo  " w{p,m} bm0          {RES,CT}     $resdir/dyturbo-0.9.6.2/results_bm0_WpmRESCT   bm0 11105 benchmark_v0_160204_WZ " >> mergeconf
    #echo  " w{p,m} bm1          {RES,CT}     $resdir/dyturbo-0.9.6.2/results_bm012_WpWm     bm1 11105 benchmark_v1_160204_WZ " >> mergeconf
    #echo  " w{p,m} bm2$qtymerge {RES3D,CT3D} $resdir/dyturbo-0.9.6.2/results_bm012_WpWm     bm2 10100 benchmark_v2_160204_WZ " >> mergeconf

    #
    echo  " z0 o1 {RES,CT}  results o1 1010 aimoments_test_160211 " >> mergeconf
    echo  " z0 o2 {RES,CT}  results o2 1010 aimoments_test_160211 " >> mergeconf




    while read -r mline
    do
        read proc inbm term indir outbm seednum prodname <<< $(echo $mline)
        outdir=results_merge/$prodname
        mkdir -p $outdir
        #
        for fres in ` eval "ls $indir/dyturbo_${proc}_*_${inbm}*t${term}_*$seednum*root"`
        do
            #MERGER="./bin/merger -X "
            MERGER="./bin/merger "
            infiles=`echo $fres | sed "s|$seednum|*|g"`
            outfilebase=`basename $fres | sed "s|$seednum.*||g;s|y-55|ym55|g;s|$inbm|$outbm|g"`
            ##
            [[ $fres =~ t*3D_ ]] &&
                MERGER="hadd -f " &&
                infiles=`echo $fres | sed "s|$seednum|*|g;s|$qtymerge|*|g"` &&
                outfilebase=`basename $fres | sed "s|$seednum.*||g;s|$qtymerge|qt0100ym55|g;"`
            #echo infiles: `ls $infiles | wc -l`
            #ls -1 $infiles
            $DRYRUN $MERGER $outdir/${outfilebase}outliers.root $infiles #&& return 0 || return 3
            #echo
        done
    done < mergeconf > mergeout 2>&1
    echo " output written in mergeout"
}



#merge_w_pt
#merge_pt_y
#merge_wwidth
#merge_grid
#merge_grid_TOT

#merge_cubatures
#merge_benchmark
merge_general
#merge_stefano

exit 0
